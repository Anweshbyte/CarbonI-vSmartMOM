using NCDatasets
using LinearAlgebra
using Plots
using Statistics
using Missings
using Glob

# Define the target latitude and longitude
target_lat_lon = [34.0, -118.25]

# Directory containing NetCDF files
directory_path = "rads"

# Get a list of all NetCDF files in the directory
file_paths = glob("*.nc", directory_path)

# Function to convert latitude and longitude to Cartesian coordinates
function latlon_to_cartesian(lat, lon, radius=6371.0)
    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(lon)
    x = radius * cos(lat_rad) * cos(lon_rad)
    y = radius * cos(lat_rad) * sin(lon_rad)
    z = radius * sin(lat_rad)
    return x, y, z
end

# Convert the target location to Cartesian coordinates
target_x, target_y, target_z = latlon_to_cartesian(target_lat_lon[1], target_lat_lon[2])

# Initialize a 2D array to store distances
distance_matrix = []

# Initialize variables to store latitude and longitude data
lat = nothing
lon = nothing

# Process the first file to get the latitude and longitude data
first_file_path = file_paths[1]
ds = Dataset(first_file_path, "r")
lat = ds.group["location"]["lat"][:]
lon = ds.group["location"]["lon"][:]
close(ds)

# Get the dimensions of the latitude and longitude arrays
n_rows, n_cols = size(lat)

# Initialize the distance matrix with the appropriate size
distance_matrix = zeros(Float64, n_rows, n_cols)

# Convert all locations to Cartesian coordinates and calculate distances
for i in 1:n_rows
    for j in 1:n_cols
        x, y, z = latlon_to_cartesian(lat[i, j], lon[i, j])
        distance_matrix[i, j] = sqrt((x - target_x)^2 + (y - target_y)^2 + (z - target_z)^2)
    end
end

# Find the index of the closest location
closest_index = argmin(distance_matrix)

pixel_x = closest_index[1]
pixel_y = closest_index[2]
start_band = 49
end_band = 55

# Function to compute normalized radiance for a given file and pixel
function radiance_dep(file_path::AbstractString, pixel_x::Int, pixel_y::Int, start_band::Int, end_band::Int)
    ds = Dataset(file_path, "r")
    radiance = ds["radiance"][:]
    point_radiance = radiance[:, pixel_x, pixel_y]
    close(ds)
    
    num_bands = length(point_radiance)
    group_size = 5
    
    normalized_radiance = similar(point_radiance)
    
    for i in 1:num_bands
        group_start = max(1, i - div(group_size, 2))
        group_end = min(num_bands, i + div(group_size, 2))
        group_values = point_radiance[group_start:group_end]
        max_indices = sortperm(group_values, rev=true)[1:2]
        mean_max_values = mean(group_values[max_indices])
        normalized_radiance[i] = point_radiance[i] / mean_max_values
    end
    return normalized_radiance[start_band:end_band]
end

# Extract dates from file names
file_dates = [match(r"(\d{8})", file_path).match for file_path in file_paths]

# Plot initialization
plot()

# Iterate over each file path
normalized_radiances = []
for (i, file_path) in enumerate(file_paths)
    # Compute normalized radiance
    norm_radiance = radiance_dep(file_path, pixel_x, pixel_y, start_band, end_band)
    push!(normalized_radiances, norm_radiance)
    
    # Calculate wavelengths for current file
    band_numbers = start_band:end_band
    wavelengths = map(x -> x, band_numbers)  # Replace with actual function if needed
    
    # Extract date from file name
    file_date = file_dates[i]
    
    # Plot normalized radiance with date as label
    plot!(wavelengths, norm_radiance, label="$file_date")
end

# Customize plot labels and save the plot
xlabel!("Wavelength (nm)")
ylabel!("Normalized Radiance")
title!("Normalized Radiance for Multiple Files")
plot!(size=(600, 1000), yticks=:auto)
savefig("/Users/arpitasen/Documents/GitHub/EMIT-Caltech/data/timeseriesplt")