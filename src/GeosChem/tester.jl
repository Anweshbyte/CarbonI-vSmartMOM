using CarbonI, Plots

geos, aero = CarbonI.loadGeos("src/GeosChem/GeosChem.yaml")
aod_salc = geos.data["aod_salc"];
aod_sala = geos.data["aod_sala"];
aod_ocpi = geos.data["aod_ocpi"];
aod_bcpi = geos.data["aod_bcpi"];
aod_so4 = geos.data["aod_so4"];
aod_strat = geos.data["aod_strat"];
rad_so4 = geos.data["radius_so4"];

heatmap(sum(aod_so4[:,:,:,1], dims=3)[:,:,1]',color=:viridis); title!("AOD_SO4");savefig("/home/cfranken/AOD_SO4.pdf")
heatmap(sum(aod_sala[:,:,:,1], dims=3)[:,:,1]',color=:viridis); title!("AOD_SAL Accumlation mode");savefig("/home/cfranken/AOD_SALA.pdf")
heatmap(sum(aod_salc[:,:,:,1], dims=3)[:,:,1]',color=:viridis); title!("AOD_SAL Coarse mode");savefig("/home/cfranken/AOD_SALC.pdf")
heatmap(sum(aod_ocpi[:,:,:,1], dims=3)[:,:,1]',color=:viridis); title!("AOD_OCPI");savefig("/home/cfranken/AOD_OCPI.pdf")
heatmap(sum(aod_bcpi[:,:,:,1], dims=3)[:,:,1]',color=:viridis);title!("AOD_BCPI");savefig("/home/cfranken/AOD_BCPI.pdf")
heatmap(sum(aod_strat[:,:,:,1], dims=3)[:,:,1]',color=:viridis);title!("AOD Stratosphere (550nm)");savefig("/home/cfranken/AODstrat.pdf")

aod_2200nm = sum(aod_so4[:,:,:,1], dims=3)[:,:,1]' ./10.0 + sum(aod_ocpi[:,:,:,1], dims=3)[:,:,1]' ./10.0 + sum(aod_strat[:,:,:,1], dims=3)[:,:,1]' ./10 + sum(aod_sala[:,:,:,1], dims=3)[:,:,1]' ./5 + sum(aod_salc[:,:,:,1], dims=3)[:,:,1]'  + sum(aod_bcpi[:,:,:,1], dims=3)[:,:,1]' ./10

total_aod_all = sum(aod_so4[:,:,:,1], dims=3)[:,:,1]' + sum(aod_sala[:,:,:,1], dims=3)[:,:,1]' + sum(aod_salc[:,:,:,1], dims=3)[:,:,1]' + sum(aod_ocpi[:,:,:,1], dims=3)[:,:,1]' + sum(aod_bcpi[:,:,:,1], dims=3)[:,:,1]' + sum(aod_strat[:,:,:,1], dims=3)[:,:,1]'

dp = CarbonI.computeColumnAveragingOperator(geos);
tropo_lev = geos.data["tropopause_level"];
xn2o = CarbonI.getColumnAverage(geos.data["N2O"], dp);
xch4 = CarbonI.getColumnAverage(geos.data["CH4"], dp);
xch4_trop = CarbonI.getTroposphericColumnAverage(geos.data["CH4"], dp, tropo_lev);
xn2o_trop = CarbonI.getTroposphericColumnAverage(geos.data["N2O"], dp, tropo_lev);

c2h6 = geos.data["C2H6"];
xc2h6 = CarbonI.getColumnAverage(c2h6, dp)

#using CarbonI, Plots

#lat = -3.5 # 27.2
#lon = -63.0 # 82.5
lat = 27.2
lon = 82.5

iLat = argmin(abs.(geos.data["lat"] .- lat))
iLon = argmin(abs.(geos.data["lon"] .- lon))
iTime = 1

FT = Float32
# Test profile extractionH2O,CO2,CH4, N2O, CO
q, T, p, p_full, new_vmr = CarbonI.get_profiles_from_geoschem(geos, iLon, iLat, iTime, ["H2O", "CO2", "CH4", "N2O","CO", "C2H6"])
aod = CarbonI.get_aerosols_from_geoschem(geos, iLon, iLat, iTime, ["aod_ocpi", "aod_so4", "aod_strat", "aod_dust", "aod_sala", "aod_salc"])
aod_radii = CarbonI.get_aerosols_from_geoschem(geos, iLon, iLat, iTime, ["radius_ocpi", "radius_so4", "radius_sala","radius_salc"])
p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr2, Δz = CarbonI.compute_atmos_profile_fields(T,convert.(FT,p),q,new_vmr)

plot(aod["aod_so4"], p_full, yflip=true, label="Sulfate Aerosols")
plot!(aod["aod_ocpi"], p_full, yflip=true, label="Organic Carbon Aerosols")
plot!(aod["aod_strat"], p_full, yflip=true, label="Stratospheric Aerosols")
plot!(aod["aod_sala"], p_full, yflip=true, label="SALA Aerosols")
plot!(aod["aod_salc"], p_full, yflip=true, label="SALC Aerosols")
plot!(aod["aod_dust"], p_full, yflip=true, label="Dust Aerosols",legend=:topright)

aeros = ("aod_so4","aod_ocpi","aod_strat","aod_sala","aod_salc","aod_dust")
for aeroi in aeros
    @show aeroi, sum(aod[aeroi])
end
Total_AOD = sum(aod["aod_so4"]) + sum(aod["aod_ocpi"]) + sum(aod["aod_strat"]) + sum(aod["aod_dust"])

heatmap(geos.data["C2H6"][:,:,1,1]'*1e9)

function compute_aod_half_level_and_weighted_radius(aod::Array{FT,3}, radius::Array{FT,3}) where {FT}
    # Get the dimensions of the input arrays
    n_lon, n_lat, n_height = size(aod)
    
    # Initialize arrays to store the results
    half_aod_levels = Array{Int, 2}(undef, n_lon, n_lat)
    weighted_avg_radius = Array{FT, 2}(undef, n_lon, n_lat)

    # Loop through each (longitude, latitude) point
    for i in 1:n_lon
        for j in 1:n_lat
            # Get the AOD and radius values for this grid point across all height levels
            aod_vals = aod[i, j, :]
            radius_vals = radius[i, j, :]
            
            # ---- Part 1: Compute the height level where AOD reaches half of total ----
            # Total AOD for this grid point
            total_aod = sum(aod_vals)
            
            # Cumulative AOD
            cumulative_aod = cumsum(aod_vals)
            
            # Find the height level where cumulative AOD reaches half the total AOD
            half_aod = total_aod / 2
            level = findfirst(cumulative_aod .>= half_aod)
            half_aod_levels[i, j] = level
            
            # ---- Part 2: Compute the AOD-weighted average radius ----
            if total_aod > 0
                # Compute weighted sum of radius using AOD as weights
                weighted_sum_radius = sum(radius_vals .* aod_vals)
                # Compute the weighted average radius
                weighted_avg_radius[i, j] = weighted_sum_radius / total_aod
            else
                # Handle cases where total AOD is zero
                weighted_avg_radius[i, j] = NaN  # or 0 if more appropriate
            end
        end
    end

    return half_aod_levels, weighted_avg_radius
end
