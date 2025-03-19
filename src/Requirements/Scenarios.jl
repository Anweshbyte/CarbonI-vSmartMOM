
using Unitful
using CarbonI

Base.@kwdef mutable struct Scenario
    ### Characteristics of the scenario 
    lat::Float64
    lon::Float64
    sza::Float64
    profile_hr
    wind_speed
    broadband_albedo::Float64

end


function reference_scenario()
    lat = 40.0
    lon = -100.0
    sza = 50.0
    wind_speed = 2.0u"m/s"
    broadband_albedo = 0.15
    # make this flat (or quartz like)
    profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, lat, lon, 7);
    sce = Scenario(lat=lat, lon=lon, sza=sza, profile_hr=profile_hr, wind_speed=wind_speed,broadband_albedo=broadband_albedo)
    return sce
end

function stressing_scenario()
    lat = 0.0
    lon = -62.0
    sza = 30.0
    wind_speed = 2.5u"m/s"
    broadband_albedo = 0.06
    profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, lat, lon, 7);
    sce = Scenario(lat=lat, lon=lon, sza=sza, profile_hr=profile_hr, wind_speed=wind_speed, broadband_albedo=broadband_albedo)
    return sce
end

function generic_scenario(lat::Float64, lon::Float64, sza::Float64, wind_speed, broadband_albedo)
    profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, lat, lon, 7);
    sce = Scenario(lat=lat, lon=lon, sza=sza, profile_hr=profile_hr, wind_speed=wind_speed, broadband_albedo=broadband_albedo)
    return sce
end


Base.@kwdef mutable struct Constants
    g = 9.81u"m/s^2"
    p = 100000u"kg/m/s^2" # aka Pa
end

Base.@kwdef mutable struct molar_mass
    co2 = 0.04401u"kg/mol"
    ch4 = 0.01604u"kg/mol"
    h2o = 0.018015u"kg/mol"
    hdo = 0.01902u"kg/mol"
    n2o = 0.044013u"kg/mol"
    co = 0.02801u"kg/mol"
    co2_iso2 = 0.04599u"g/mol" # need
    c2h6 = 0.03007u"kg/mol"
    air = 0.02897u"kg/mol"
end

