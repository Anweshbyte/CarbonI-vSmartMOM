using Unitful 

#ensure_artifact_installed("cross_sections", find_artifacts_toml(@__DIR__), quiet_download = true)
function get_artifacts_path(name)
    # This returns the full path to the current Project.toml or parent environment file
    proj_toml = Base.current_project()
    # Extract just the directory part
    proj_dir = dirname(proj_toml)
    return ensure_artifact_installed(name, find_artifacts_toml(proj_dir), quiet_download = true)
end
#ensure_artifact_installed("cross_sections", find_artifacts_toml(@__DIR__), quiet_download = true)
# Load cross section tables (can toggle to also use ABSCO, doesn't cover full range though)
function loadXSModels(path=CarbonI.xs_folder)
    co2 = load_interpolation_model(path*"/co2_model.jld2")
    co2_iso2 = load_interpolation_model(path*"/co2_model_iso2.jld2")
    # Commenting out the ABSCO CO2 table for now, as I want 13CO2 and CO2 in the same table
    #co2_ = load_interpolation_model(artifact"cross_sections"*"/sco2_v52.jld2")
    #etp0 = extrapolate(co2_.itp, 0.0)
    #co2__ = InterpolationModel( etp0, co2_.mol, co2_.iso, co2_.ν_grid, co2_.p_grid, co2_.t_grid);
    ch4 = load_interpolation_model(path*"/ch4_model.jld2")
    h2o = load_interpolation_model(path*"/h2o_model_tccon.jld2")
    hdo = load_interpolation_model(path*"/hdo_model.jld2")
    n2o = load_interpolation_model(path*"/n2o_model.jld2")
    co  = load_interpolation_model(path*"/co_model.jld2")
    c2h6  = load_interpolation_model(path*"/c2h6_model.jld2")
    #return co2__, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6
    return co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6
end

#using Pkg.Artifacts


abstract type InstrumentOperator end

"Struct for an atmospheric profile"
struct AtmosphericProfile{FT}
    lat::FT
    lon::FT
    psurf::FT
    T::Array{FT,1}
    q::Array{FT,1}
    p::Array{FT,1}
    p_levels::Array{FT,1}
    vmr_h2o::Array{FT,1}
    vcd_dry::Array{FT,1}
    vcd_h2o::Array{FT,1}
end;

"Struct for Kernel Instrument Function"
@with_kw struct KernelInstrument{FT} <: InstrumentOperator 
    kernel::FT
    ν_out::Array
end;

"Read atmospheric profile (just works for our file, can be generalized"
function read_atmos_profile(file::String, lat::Real, lon::Real, timeIndex; g₀=9.807)
    @assert 1 <= timeIndex <= 4

    ds = Dataset(file)

    # See how easy it is to actually extract data? Note the [:] in the end reads in ALL the data in one step
    lat_   = ds["YDim"][:]
    lon_   = ds["XDim"][:]
    
    FT = eltype(lat_)
    lat = FT(lat)
    lon = FT(lon)
    
    # Find index (nearest neighbor, one could envision interpolation in space and time!):
    iLat = argmin(abs.(lat_ .- lat))
    iLon = argmin(abs.(lon_ .- lon))
    @show ds["T"][iLon,iLat, :, timeIndex]
    # Temperature profile
    T    = convert(Array{FT,1}, ds["T"][ iLon,iLat, :, timeIndex])
    # specific humidity profile
    q    = convert(Array{FT,1}, ds["QV"][iLon,iLat,  :, timeIndex])
    
    # Surafce pressure
    psurf = convert(FT, ds["PS"][iLon, iLat, timeIndex])
    
    # AK and BK global attributes (important to calculate pressure half-levels)
    ak = ds.attrib["HDF_GLOBAL.ak"][:]
    bk = ds.attrib["HDF_GLOBAL.bk"][:]

    p_half = (ak + bk * psurf)
    p_full = (p_half[2:end] + p_half[1:end - 1]) / 2
    close(ds)
    
    # Avogradro's number:
    Na = 6.0221415e23;
    # Dry and wet mass
    dryMass = 28.9647e-3  / Na  # in kg/molec, weighted average for N2 and O2
    wetMass = 18.01528e-3 / Na  # just H2O
    ratio = dryMass / wetMass 
    n_layers = length(T)
    # also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dryMass + vmr_h2o[i] * wetMass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o)
end;

"Read atmospheric profile (just works for our file, can be generalized"
function read_atmos_profile_MERRA2(file::String, lat::Real, lon::Real, timeIndex; g₀=9.8196)
    @assert 1 <= timeIndex <= 8

    ds = Dataset(file)

    # See how easy it is to actually extract data? Note the [:] in the end reads in ALL the data in one step
    lat_   = ds["lat"][:]
    lon_   = ds["lon"][:]
    
    FT = eltype(lat_)
    lat = FT(lat)
    lon = FT(lon)
    
    # Find index (nearest neighbor, one could envision interpolation in space and time!):
    iLat = argmin(abs.(lat_ .- lat))
    iLon = argmin(abs.(lon_ .- lon))
    @show ds["T"][iLon,iLat, :, timeIndex]
    # Temperature profile
    T    = convert(Array{FT,1}, ds["T"][ iLon,iLat, :, timeIndex])
    # specific humidity profile
    q    = convert(Array{FT,1}, ds["QV"][iLon,iLat,  :, timeIndex])
    
    # Surafce pressure
    psurf = convert(FT, ds["PS"][iLon, iLat, timeIndex])
    
    delp = ds["DELP"][iLon, iLat,:, timeIndex]
    # AK and BK global attributes (important to calculate pressure half-levels)
    #ak = ds.attrib["HDF_GLOBAL.ak"][:]
    #bk = ds.attrib["HDF_GLOBAL.bk"][:]
    a = cumsum(delp) .-0.1
    p_half = [0.0; a]
    p_full = (p_half[2:end] + p_half[1:end - 1]) / 2
    close(ds)
    
    # Avogradro's number:
    Na = 6.0221415e23;
    # Dry and wet mass
    dryMass = 28.9647e-3  / Na  # in kg/molec, weighted average for N2 and O2
    wetMass = 18.01528e-3 / Na  # just H2O
    ratio = dryMass / wetMass 
    n_layers = length(T)
    # also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dryMass + vmr_h2o[i] * wetMass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o)
end;

"Read atmospheric profile (just works for our file, can be generalized"
function generate_atmos_profile(T,p,q; g₀=9.807)
    FT = eltype(T)
    @assert length(p) ==  length(T)+1 == length(q)+1
    
    # Surafce pressure
    psurf = p[end]
    p_half = p 
    p_full = (p_half[2:end] + p_half[1:end - 1]) / 2
    
    # Avogradro's number:
    Na = FT(6.0221415e23);
    # Dry and wet mass
    dryMass = FT(28.9647e-3)  / Na  # in kg/molec, weighted average for N2 and O2
    wetMass = FT(18.01528e-3) / Na  # just H2O
    ratio = dryMass / wetMass 
    n_layers = length(T)
    # also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dryMass + vmr_h2o[i] * wetMass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end
    return AtmosphericProfile(FT(0.0), FT(0.0), psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o)
end;

"Computes cross section matrix for arbitrary number of absorbers"
function compute_profile_crossSections(profile::AtmosphericProfile, hitranModels, ν::AbstractRange{<:Real})
    nGases   = length(hitranModels)
    nProfile = length(profile.p)
    FT = eltype(profile.T)
    n_layers = length(profile.T)

    σ_matrix = zeros(FT, (length(ν), n_layers, nGases))
    for i = 1:n_layers
        p_ = profile.p[i] / 100 # in hPa
        T_ = profile.T[i]
        for j = 1:nGases
            σ_matrix[:,i,j] = absorption_cross_section(hitranModels[j], ν, p_, T_; wavelength_flag=true);
        end
    end
    return σ_matrix
end;

"Creates a Gaussian Kernel"
function gaussian_kernel(FWHM, res; fac=5)
    co = 2.355
    width = FWHM / res / co
    extent = ceil(fac * width)
    d = Normal(0, width)
    kernel = centered(pdf.(d, -extent:extent))
    return kernel ./ sum(kernel)
end;

"Creates a box Kernel"
function box_kernel(FWHM, res; fac=5)
    
    width = FWHM/res 
    extent = ceil(fac * width)
    d = Uniform(-width/2, width/2)
    kernel = centered(pdf.(d, -extent:extent))
    return kernel ./ sum(kernel)
end;

"Rescales x-axis to go from [-1,1]"
function rescale_x(a)
    a = a .- mean(a);
    a = a / (a[end] - a[1]) * 2;
end;

"Convolve and resample"
function conv_spectra(m::KernelInstrument, ν, spectrum)
    s = imfilter(spectrum, m.kernel, Algorithm.FIR())
    interp_cubic = CubicSplineInterpolation(ν, s)
    return interp_cubic(m.ν_out)
end;
    
function generate_conv_matrix(m::KernelInstrument, wl_in, Δwl)
    s = zeros(length(m.ν_out), length(wl_in))
    x_ = Δwl .* ((1:length(m.kernel.parent)) .+ m.kernel.offsets[1])
    interp_cubic = CubicSplineInterpolation(x_, m.kernel.parent,extrapolation_bc = 0)
    for i in axes(s,2), j in axes(s,1)
        s[j,i] = interp_cubic(m.ν_out[j] .- wl_in[i])
    end
    return s
end;

"Reduce profile dimensions"
function reduce_profile(n::Int, profile::AtmosphericProfile{FT}, σ_matrix, gasProfiles) where {FT}
    @assert n < length(profile.T)
    @unpack lat, lon, psurf = profile
    # New rough half levels (boundary points)
    #a = range(0, maximum(profile.p_levels), length=n + 1)
    a = reduce_pressure_levels(profile.p_levels, n+1)
    n = length(a)-1
    dims = size(σ_matrix)
    #FT = eltype(σ_matrix)
    σ_matrix_lr = zeros(FT, dims[1], n, dims[3])
    T = zeros(FT, n);
    q = zeros(FT, n);
    p_full = zeros(FT, n);
    p_levels = zeros(FT, n + 1);
    vmr_h2o  = zeros(FT, n);
    vcd_dry  = zeros(FT, n);
    vcd_h2o  = zeros(FT, n);
    indis = []
    gasProfile_low_res = [zeros(n) for i in 1:length(gasProfiles)]
    for i = 1:n
        ind = findall(a[i] .< profile.p .<= a[i + 1]);
        push!(indis,ind)
        #@show ind
        h = profile.vcd_dry[ind]
        h ./= sum(h)
        # This has to be weighted by vcd_dry!
        for iGas=1:dims[3]
            #@show iGas, ind
            #@show size(gasProfiles[iGas])
            # Weigh by gas VCD high res profile:
            gasP = gasProfiles[iGas][ind] .* profile.vcd_dry[ind]
            gas_vmr = sum(gasP)/sum(profile.vcd_dry[ind])
            gasProfile_low_res[iGas][i] = gas_vmr
            gasP ./= sum(gasP)
            #@show gas_vmr, mean(gasProfiles[iGas][ind])
            σ_matrix_lr[:,i,iGas] = σ_matrix[:,ind,iGas] *gasP
            #@show sum(gasP)
        end
        p_levels[i] = a[i]
        p_levels[i + 1] = a[i+1]
        p_full[i] = mean(profile.p[ind])
        T[i] = profile.T[ind]' * h
        q[i] = profile.q[ind]' * h
        vmr_h2o[i] = profile.vmr_h2o[ind]' * h
        vcd_dry[i] = sum(profile.vcd_dry[ind])
        vcd_h2o[i] = sum(profile.vcd_h2o[ind])
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_levels, vmr_h2o, vcd_dry, vcd_h2o), σ_matrix_lr, indis, gasProfile_low_res
end;

"Reduce profile dimensions"
function reduce_profile_vcd(n::Int, profile::AtmosphericProfile{FT}, σ_matrix, gasProfiles) where {FT}
    @assert n < length(profile.T)
    @unpack lat, lon, psurf = profile
    # New rough half levels (boundary points)
    #a = range(0, maximum(profile.p_levels), length=n + 1)
    a = reduce_pressure_levels(profile.p_levels, n+1)
    n = length(a)-1
    dims = size(σ_matrix)
    #FT = eltype(σ_matrix)
    σ_matrix_lr = zeros(FT, dims[1], n, dims[3])
    T = zeros(FT, n);
    q = zeros(FT, n);
    p_full = zeros(FT, n);
    p_levels = zeros(FT, n + 1);
    vmr_h2o  = zeros(FT, n);
    vcd_dry  = zeros(FT, n);
    vcd_h2o  = zeros(FT, n);
    indis = []
    gasProfile_low_res = [zeros(n) for i in 1:length(gasProfiles)]
    for i = 1:n
        ind = findall(a[i] .< profile.p .<= a[i + 1]);
        push!(indis,ind)
        #@show ind
        h = profile.vcd_dry[ind]
        h ./= sum(h)
        # This has to be weighted by vcd_dry!
        for iGas=1:dims[3]
            @show iGas, ind
            #@show size(gasProfiles[iGas])
            # Weigh by gas VCD high res profile:
            gasP = gasProfiles[iGas][ind] .* profile.vcd_dry[ind]
            gas_vmr = sum(gasP)/sum(profile.vcd_dry[ind])
            gasProfile_low_res[iGas][i] = sum(gasP)
            gasP ./= sum(gasP)
            @show gas_vmr, mean(gasProfiles[iGas][ind])

            # Weighted average of cross section:
            σ_matrix_lr[:,i,iGas] = σ_matrix[:,ind,iGas] *gasP
            #@show sum(gasP)
        end
        p_levels[i] = a[i]
        p_levels[i + 1] = a[i+1]
        p_full[i] = mean(profile.p[ind])
        T[i] = profile.T[ind]' * h
        q[i] = profile.q[ind]' * h
        vmr_h2o[i] = profile.vmr_h2o[ind]' * h
        vcd_dry[i] = sum(profile.vcd_dry[ind])
        vcd_h2o[i] = sum(profile.vcd_h2o[ind])
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_levels, vmr_h2o, vcd_dry, vcd_h2o), σ_matrix_lr, indis, gasProfile_low_res
end;

function reduce_pressure_levels(pressures, n::Int)
    # Ensure the input pressures are sorted
    sorted_pressures = sort(pressures)
    # Number of pressure levels in the input array
    m = length(sorted_pressures)

    # Check if the desired number of levels is greater than the available levels
    if n > m
        error("Desired number of levels (n) is greater than the available levels in the input array")
    end

    # Create a range of equally spaced values between the minimum and maximum pressures
    p_min = minimum(sorted_pressures)
    p_max = maximum(sorted_pressures)
    target_pressures = range(p_min, stop=p_max, length=n)

    # Find the closest pressure levels in the sorted array to the target pressures
    new_pressures = Array{Float64,1}(undef, n)
    for i in 1:n
        diff = sorted_pressures .- target_pressures[i]
        #@show diff, sorted_pressures
        indi = findall(diff .>= 0)
        closest_index = indi[1]
        #@show i, target_pressures[i]/100, closest_index
        new_pressures[i] = sorted_pressures[closest_index]
        #@show new_pressures[i]
    end
    #@show new_pressures
    # Remove duplicates in case multiple target pressures are closest to the same pressure level
    new_pressures = unique(new_pressures)

    #@show length(new_pressures)
    # If we have fewer levels than desired due to duplicates, refill the remaining levels
    #while length(new_pressures) < n
    #    remaining_pressures = setdiff(sorted_pressures, new_pressures)
    #    new_pressures = vcat(new_pressures, remaining_pressures[1])
    #    @show new_pressures
    #end

    return sort(new_pressures)
end


function jacobs_eq(gtype, wind_speed, enh, pixel_size; style="req")

    m = CarbonI.molar_mass()
    C = CarbonI.Constants()

    #             kg/(m s^2)   m/s        s/hr *  ppb    * s^2 / m * m
    kg_per_hour = C.p * wind_speed * 3600u"s/hr" * enh / 1e9 / C.g * pixel_size


    if gtype == "ch4"
        kg_per_hour *= m.ch4/m.air 
    elseif gtype == "co2"
        kg_per_hour *= m.co2/m.air 
    elseif gtype == "co"
        kg_per_hour *= m.co/m.air 
    elseif gtype == "h2o"
        kg_per_hour *= m.h2o/m.air
    elseif gtype == "n2o"
        kg_per_hour *= m.n2o/m.air
    elseif gtype == "hdo"
        kg_per_hour *= m.hdo/m.air
    elseif gtype == "c2h6"
        kg_per_hour *= m.c2h6/m.air
    elseif gtype == "co2_iso2"
        kg_per_hour *= m.co2_iso2/m.air
    else
        error("Gas type not supported:: $gtype")
    end

    if occursin("detect",style)
        kg_per_hour *= 2
    elseif occursin("req",style)
        kg_per_hour *= 1
    elseif occursin("quant",style)
        kg_per_hour *= 5
    else
        error("Style not supported:: $style")
    end
    
    
    return kg_per_hour



end