using Pkg
using CSV
using DataFrames
using Dates
using HTTP
using JSON
using Base.Iterators
using IterTools
using GeometryBasics
using Shapefile
using GeoJSON
using GeoInterface
using DataFrames
using SplitApplyCombine

page_num = 0
page_size = 2000 # CMR page size limit

function poltocart(lat, lon, radius=6371.0)
    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(lon)
    x = radius * cos(lat_rad) * cos(lon_rad)
    y = radius * cos(lat_rad) * sin(lon_rad)
    z = radius * sin(lat_rad)
    return x, y, z
end

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

function band_to_wavelength(band_number::Int)
    min_wavelength = 381
    max_wavelength = 2493
    total_bands = 285
    wavelength_range = max_wavelength - min_wavelength
    wavelength_per_band = wavelength_range / total_bands
    wavelength = min_wavelength + (band_number - 1) * wavelength_per_band
    return wavelength
end

function extract_emit_data(lat::Float64, lon::Float64, start_date::DateTime, end_date::DateTime)

    doi = "10.5067/EMIT/EMITL1BRAD.001" #Required, you may find the DOI for the dataset here: https://lpdaac.usgs.gov/product_search/?query=emit&view=cards&sort=title"
    file_path = "/Users/arpitasen/Desktop/EMIT-Caltech/data/EMIT_extract.csv" # File path where we store the extracted data

    cmrurl = "https://cmr.earthdata.nasa.gov/search/"
    doisearch = cmrurl * "collections.json?doi=" * doi
    response1 = HTTP.get(doisearch)
    json_response = JSON.parse(String(response1.body))
    concept_id = json_response["feed"]["entry"][1]["id"]

    dt_format = dateformat"yyyy-mm-ddTHH:MM:SSZ"
    start_str = Dates.format(start_date, dt_format)
    end_str = Dates.format(end_date, dt_format)
    temporal_str = start_str * "," * end_str

    point_str = string(lon) * "," * string(lat)

    granule_arr = []
    output = DataFrame(Asset = String[], URL = [], CC = String[], Poly = [])
    while true
        global page_num += 1
        cmr_param = Dict(
            "collection_concept_id" => concept_id,
            "page_size" => string(page_size),
            "page_num" => string(page_num),
            "temporal" => temporal_str,
            "point" => point_str
        )
        granulesearch = cmrurl * "granules.json"
        response = HTTP.post(granulesearch, [], HTTP.Form(cmr_param))
        granules = JSON.parse(String(response.body))["feed"]["entry"] # Returns a dictionary with all the measurements and parameters
        if isempty(granules)  # Exit loop if no granules are returned
            break
        end
        for g in granules
            granule_urls = []
            global cloud_cover = get(g, "cloud_cover", "")
            if haskey(g, "polygons") # Extract and process the polygons
                polygons = get(g, "polygons", "")
                for poly in polygons
                    terms = split(poly[1])
                    ltln = []
                    for i in 1:2:length(terms)-1
                        push!(ltln, terms[i] * " " * terms[i+1])
                    end
                    global granule_poly = [join(reverse(split(s, ' ')), ' ') for s in ltln]
                end
            end
            global granule_urls = [x["href"] for x in g["links"] if occursin("https", x["href"]) && occursin(".nc", x["href"]) && !occursin(".dmrpp", x["href"])]
            df = DataFrame(URL = [granule_urls], CC = cloud_cover, Poly = [granule_poly])
            df = DataFrames.flatten(df,[:URL])
            df[!,"Asset"] = last.(split.(df[:, "URL"], '/'))
            df = df[:, ["Asset", setdiff(names(df), ["Asset"])...]]
            append!(output, df)
        end
    end
    CSV.write(file_path, output)
end
