using Plots

function band_to_wavelength(band_number::Int)
    min_wavelength = 381
    max_wavelength = 2493
    total_bands = 285
    wavelength_range = max_wavelength - min_wavelength
    wavelength_per_band = wavelength_range / total_bands
    wavelength = min_wavelength + (band_number - 1) * wavelength_per_band
    return wavelength
end

#User Params =>

file_paths = [
    "/Users/arpitasen/Desktop/Sanghavi/Data/EMIT_L1B_RAD_001_20230430T164419_2312011_002.nc",
    # "/Users/arpitasen/Desktop/Sanghavi/Data/EMIT_L1B_RAD_001_20230614T224902_2316515_008.nc",
    # "/Users/arpitasen/Desktop/Sanghavi/Data/EMIT_L1B_RAD_001_20231022T192739_2329513_007.nc",
    # "/Users/arpitasen/Desktop/Sanghavi/Data/EMIT_L1B_RAD_001_20231026T175222_2329912_010.nc"
]
pixel_x = 1050
pixel_y = 500
start_band = 49
end_band = 55

normalized_radiances = []

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
for (i, file_path) in enumerate(file_paths)
    # Compute normalized radiance (replace with actual function call)
    norm_radiance = radiance_dep(file_path, pixel_x, pixel_y, start_band, end_band)
    push!(normalized_radiances, norm_radiance)
    
    # Calculate wavelengths for current file
    band_numbers = start_band:end_band
    wavelengths = map(band_to_wavelength, band_numbers)
    
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
# # savefig("/Users/arpitasen/Desktop/Sanghavi/Data/normalized_radiances_plot_new.png")
