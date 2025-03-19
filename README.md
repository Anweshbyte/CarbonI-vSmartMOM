# Carbon-I

Carbon-I Julia and Python support packages.  This repo is under active development and subject to lots of change during the development phase.


> ### Installation Instructions
> 
> A future version of this repository will include all components as a fully installable package.  For now, we recommend:
> - [Install Julia 1.11.2 ](https://julialang.org/downloads/)
> - Optionally set your JULIA_DEPOT_PATH to a specific location for installation (defaults to ~/.julia)
> - Clone this repository to a location (referred to as CARBONI_HOME below)
> - Then run:
```
export JULIA_PROJECT=${CARBONI_HOME}
julia -e 'using Pkg; Pkg.precompile()'
```

# EMIT Data extraction using the CMR APIs

## Introduction
The Earth Surface Mineral Dust Source Investigation (EMIT) data can be accessed programmatically using NASA's Common Metadata Repository (CMR) API. The CMR is a metadata system that catalogs Earth Science data and associated metadata records. The CMR API allows users to search through its vast metadata holdings using various parameters and keywords. The following code builds the concept ID:

        cmrurl = "https://cmr.earthdata.nasa.gov/search/"
        doisearch = cmrurl * "collections.json?doi=" * doi
        response1 = HTTP.get(doisearch)
        json_response = JSON.parse(String(response1.body))
        concept_id = json_response["feed"]["entry"][1]["id"]

## Concept ID and DOI
To search for the EMIT dataset, you need NASA EarthData's unique Concept ID for the dataset. This Concept ID can be obtained using the Digital Object Identifier (DOI) for the dataset. DOIs can be found by clicking the `Citation` link on the [LP DAAC's EMIT Product Pages](https://lpdaac.usgs.gov/product_search/?query=emit&view=cards&sort=title).

### Concept ID for EMIT L2A Reflectance Dataset
The unique NASA-given Concept ID for the EMIT L2A Reflectance dataset is essential for retrieving relevant files (or granules).

## Handling Pagination
The CMR API has a limit of 2000 results per page and a maximum of 1 million granules matched. To handle pagination, use the `page_num` parameter to loop through the search result pages.

## Creating a DataFrame with the Resulting Links
A 'DataFrame` can be used to store the download URLs and geometries of each file. The EMIT L2A Reflectance and Uncertainty and Mask collection contains three assets per granule: reflectance, reflectance uncertainty, and masks. By printing the list, you can see the three assets corresponding to a single polygon.

## EMIT product naming conventions
Level 1 and 2 data product filenames will adhere to the following naming convention:

    EMIT_L1B_RAD_001_20220101T083015_2200105_007
    EMIT – Sensor
    L1B – Product Level
    RAD – Product Type
    001 – Product Version Number
    20220101 – Date of Acquisition (YYYYMMDD)
    T083015 – Time of Acquisition (HHMMSS) (in UTC)
    2200105 – Orbit identification number (YYDOYNN), where NN is a two-digit incrementing orbit number for each day
    007 – Scene identification number

An orbit is referred to as a scene and can cover thousands of kilometers down-track. For the sake of convenience, each scene is broken into granules of 1280 down-track lines which can be reassembled seamlessly into scenes. The last granule can be extended up to 2559 lines to prevent having a granule with very few lines.
