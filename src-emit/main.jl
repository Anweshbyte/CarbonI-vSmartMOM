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
using Downloads
using Libdl

include("utils.jl")

# Start by extracting the data into a CSV file
extract_emit_data(34.0, -118.25, DateTime(2022, 9, 3), DateTime(2024, 3, 9, 23, 23, 59))

csv_file = "data/EMIT_extract.csv"
df = CSV.read(csv_file, DataFrame)

for row in eachrow(df)
    url = row.URL
    if occursin("RAD", url)
        run(`open $url`)
    end
end





