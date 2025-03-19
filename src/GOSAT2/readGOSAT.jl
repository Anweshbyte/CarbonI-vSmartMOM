using NCDatasets

function read_gosat2(file)
    g = Dataset(file)
    rad = g.group["SoundingData"].group["Radiance"]
    info = g.group["SoundingData"].group["WavenumberInfo"]

    sza = g.group["SoundingGeometry"]["solarZenith"][:]
    lat = g.group["SoundingGeometry"]["latitude"][:]
    lon = g.group["SoundingGeometry"]["longitude"][:]
    # Load more later
    #vza = g.group["SoundingGeometry"].group["sensorZenith"][:]
    wn = range(info["beginWN"][5], step=info["deltaWN"][5], length=info["numWN"][5])
    ind = findall(2390 .> 1e7./wn .> 1950)
    rad_band3 = 0.5*rad["band3P"][1,:,ind] + 0.5*rad["band3S"][1,:,ind]
    return rad_band3, sza, wn[ind], lat, lon
end

