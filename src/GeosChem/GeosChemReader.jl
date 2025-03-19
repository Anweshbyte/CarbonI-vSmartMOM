using NCDatasets
using YAML

"""
    struct GeosData
    A struct to hold the data read from netCDF files.

    Fields:
    - data::Dict{String, Any}: A dictionary where the keys are internal variable names and the values are their corresponding data.
"""
struct GeosData
    data::Dict{String, Any}
end

"""
    loadGeos(config_path::String) -> GeosData

Load data from netCDF files based on the configuration specified in a YAML file.

Arguments:
- `config_path::String`: Path to the YAML configuration file.

Returns:
- `GeosData`: A struct containing the loaded data.

The YAML configuration file should specify the file paths and the variables to read, along with their internal variable names. It also ensures that the filenames are identical apart from the second string when split using '.'.
"""
function loadGeos(config_path::String)
    config = YAML.load_file(config_path)
    all_data = Dict{String, Any}()
    file_paths = [file_info["path"] for file_info in config["files"]]

    # Test if all filenames are identical apart from the second string if split using '.'
    base_names = [split(filepath, ".") for filepath in file_paths]
    for i in 1:length(base_names) - 1
        @assert length(base_names[i]) == length(base_names[i+1]) "Filenames have different lengths when split."
        for j in 1:length(base_names[i])
            if j != 2
                @assert base_names[i][j] == base_names[i+1][j] "Filenames differ at component $j: $(base_names[i][j]) vs $(base_names[i+1][j])"
            end
        end
    end

    for file_info in config["files"]
        file_path = file_info["path"]
        variables = file_info["variables"]
        ds = NCDataset(file_path)

        for (file_var, internal_var) in variables
            all_data[internal_var] = ds[file_var][:]
        end

        close(ds)
    end

    return GeosData(all_data), config["aerosols"]
end

"""
    getColumnAverage(gas, dp) -> Array{Float64, 2}

Calculate the column average of a gas species.

Arguments:
- `gas`: The gas concentration array with dimensions (latitude, longitude, levels, time).
- `dp`: The pressure difference array.

Returns:
- `Array{Float64, 2}`: The column average of the gas species with dimensions (latitude, longitude).
"""
function getColumnAverage(gas, dp)
    xgas = zeros(size(gas, 1), size(gas, 2))
    for i in 1:size(gas, 1)
        for j in 1:size(gas, 2)
            xgas[i, j] = gas[i, j, :, 1]' * dp
        end
    end
    return xgas
end


"""
    getTroposphericColumnAverage(gas, dp, tropoLev) -> Array{Float64, 2}

Calculate the tropospheric column average of a gas species up to the tropopause level, accounting for non-integer tropopause levels.

Arguments:
- `gas`: The gas concentration array with dimensions (latitude, longitude, levels, time).
- `dp`: The pressure difference array.
- `tropoLev`: The tropopause level array with dimensions (latitude, longitude, time).

Returns:
- `Array{Float64, 2}`: The tropospheric column average of the gas species with dimensions (latitude, longitude).

The function handles non-integer tropopause levels by including the last level, weighted according to the fractional part of the tropopause level index.
"""
function getTroposphericColumnAverage(gas, dp, tropoLev)
    xgas = zeros(size(gas, 1), size(gas, 2))
    for i in 1:size(gas, 1)
        for j in 1:size(gas, 2)
            tropo_level = tropoLev[i, j, 1]
            int_level = floor(Int, tropo_level)
            frac_level = tropo_level - int_level

            if int_level > 0
                # Include all levels up to the integer part of tropo_level
                xgas[i, j] = gas[i, j, 1:int_level, 1]' * dp[1:int_level]
                # Add the fractional part of the next level
                if int_level < size(gas, 3)
                    xgas[i, j] += gas[i, j, int_level + 1, 1] * dp[int_level + 1] * frac_level
                end
                # Normalize by the sum of the included dp values
                dp_sum = sum(dp[1:int_level]) + dp[int_level + 1] * frac_level
                xgas[i, j] /= dp_sum
            else
                # If the tropo_level is less than or equal to the first level, use only the fractional part
                xgas[i, j] = gas[i, j, 1, 1] * dp[1] * frac_level
                xgas[i, j] /= dp[1] * frac_level
            end
        end
    end
    return xgas
end


"""
    computeColumnAveragingOperator(geos) -> Array{Float64}

Compute the column averaging operator `dp` from the Geos files.

Arguments:
- `geos::GeosData`: All Geos Data.

Returns:
- `Array{Float64}`: The column averaging operator `dp`.

This function computes the following:
1. `levi = ai / 1000 + bi`
2. `dp = abs.(diff(levi))` (ignore the need to convert to dry pressure so far))
3. Normalize `dp` by dividing by the sum of its elements (note, we currently use the "wet" dp)
"""
function computeColumnAveragingOperator(geos::GeosData) 
    ai = geos.data["ai"][:];
    bi = geos.data["bi"][:];
    levi = ai / 1000 .+ bi
    dp = abs.(diff(levi))
    dp /= sum(dp)
    return dp
end