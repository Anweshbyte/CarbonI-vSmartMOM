using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA, Plots, Distributions

Δwl = 0.005
wl = 2035:Δwl:2385

# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix(wl)

R_conv_carbonI_dict = Dict{Tuple{Float64, Float64,Float64, Float64}, Vector{Float64}}()
#R_conv_emit_dict = Dict{Tuple{Float64, Float64}, Vector{Float64}}()

parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i_nowater.yaml")

# Set of runs:
aods = exp.(-8:0.2:3)
albedos = exp.(-6:0.2:0)
szas = 20.0:20:60.0
p_aeros = 750.0:100:950.0
all = length(aods) * length(albedos) *  length(szas) * length(p_aeros)
@show "Runtime in hours" all*2/60/60
for sza in szas
    for p_aero in p_aeros
        parameters.scattering_params.rt_aerosols[1].profile = Normal(p_aero, 50.0)
        parameters.sza = sza
        model = model_from_parameters(parameters);
        aod_profile_normalized = model.τ_aer[1][1,:] / sum(model.τ_aer[1][1,:])
        model.τ_rayl[1] .= 0.0
        for aod in aods
            for albedo in albedos
                println("Iter for SZA, p_aero, albedo, aod: ")
                @show sza, p_aero,albedo, aod
                model.τ_aer[1][1,:] = aod * aod_profile_normalized
                model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(albedo)
                a = rt_run(model)
                aod_band = sum(model.τ_aer[1]); 

                ν = parameters.spec_bands[1];
                Δν = mean(diff(ν));
                gridi = ν[1]:Δν:ν[end]+10eps()
                
                R = a[1][1,1,:];
                
                rad_inter = CubicSplineInterpolation(gridi, R);
                R_conv_carbonI = cM*rad_inter(1e7./wl);
                

                key = (sza, p_aero, aod, albedo) 
                R_conv_carbonI_dict[key] = R_conv_carbonI
            end
        end
    end
end

@save "simulated_rads_all_noH2OnoRayleigh.jld2" R_conv_carbonI_dict 
#  @load "simulated_rads_all.jld2" R_conv_carbonI_dict

