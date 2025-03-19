function compute_atmos_profile_fields(T::AbstractArray{FT,1}, p_half::AbstractArray{FT,1}, q, vmr; g₀=9.807) where FT
    
    # convert q from g/kg to kg/kg
    q = q ./ FT(1000)
    #FT = eltype(T)
    Nₐ = FT(6.02214179e+23)
    R  = FT(8.3144598)
    # Calculate full pressure levels
    p_full = (p_half[2:end] + p_half[1:end-1]) / 2

    # Dry and wet mass
    dry_mass = FT(28.9644e-3)    # in kg/molec, weighted average for N2 and O2
    wet_mass = FT(18.01534e-3)   # just H2O
    ratio = dry_mass / wet_mass
    n_layers = length(T)

    # Also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )
    Δz      = zeros(FT, n_layers)
    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i]/(1-q[i]) * ratio # dry_mass/(dry_mass-wet_mass*(1-1/q[i]))
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd = Nₐ * Δp / (M  * g₀ * 100^2) * 100
        vcd_dry[i] = vmr_dry    * vcd   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * vcd
        Δz[i] =  (log(p_half[i + 1]) - log(p_half[i])) / (g₀ * M  / (R * T[i]) )
        #@show Δz, T[i], M, Δp
    end

    # TODO: This is still a bit clumsy:
    new_vmr = Dict{String, Union{Real, Vector}}()


    for molec_i in keys(vmr)
        if vmr[molec_i] isa AbstractArray
            if length(vmr[molec_i]) == length(p_full)
                new_vmr[molec_i] = vmr[molec_i]
            else
                @info "Warning, make sure that the VMR is interpolated correctly! Right now, it might be tricky"
                pressure_grid = collect(range(minimum(p_full), maximum(p_full), length=length(vmr[molec_i])))
                interp_linear = LinearInterpolation(pressure_grid, vmr[molec_i])
                new_vmr[molec_i] = [interp_linear(x) for x in p_full]
            end
        else
            new_vmr[molec_i] = vmr[molec_i]
        end
    end
    
    return p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz

end

function get_profiles_from_geoschem(geos::GeosData, iLon::Int, iLat::Int, iTime::Int, gasStrings)
    # Make sure to go from TOA to Bottom
    q = reverse(geos.data["q"][iLon, iLat, :, iTime]);
    T = reverse(geos.data["T"][iLon, iLat, :, iTime]);
    ilev = geos.data["ilev"]
    lev = geos.data["lev"]
    psurf = geos.data["psurf_start"][iLon, iLat, iTime]
    p = reverse(ilev * psurf)
    p_full = reverse(lev * psurf)
    new_vmr = Dict{String, Union{Real, Vector}}()
    for gas in gasStrings
        new_vmr[gas] = reverse(geos.data[gas][iLon, iLat, :, iTime])
    end
    return q, T, p, p_full, new_vmr
end

function get_aerosols_from_geoschem(geos::GeosData, iLon::Int, iLat::Int, iTime::Int, aerosolStrings)

    new_aod = Dict{String, Union{Real, Vector}}()
    for aero in aerosolStrings
        new_aod[aero] = reverse(geos.data[aero][iLon, iLat, :, iTime])
    end
    return new_aod
end