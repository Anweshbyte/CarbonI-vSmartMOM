using CairoMakie
set_theme!(theme_ggplot2())

using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA, Plots, Distributions
include("./src/ProxyTests/setupTestFits_generic.jl")
include("./src/ProxyTests/createMicroWindows.jl")
include("./src/Plots/CI_colors.jl")

Δwl = 0.005
wl = 2035:Δwl:2385

# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)


parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")
aod = 0.1
p_aero = 900.0
albedo = 0.1
sza = 30.0

parameters.scattering_params.rt_aerosols[1].profile = Normal(p_aero, 50.0)
parameters.sza = sza
model = model_from_parameters(parameters);
aod_profile_normalized = model.τ_aer[1][1,:] / sum(model.τ_aer[1][1,:])
#model.τ_rayl[1] .= 0.0
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

#indLR1 = 8:260
#indLR2 = 261:495
x, yy, S, A, K = ii_mw1(R_conv_carbonI[indLR1]);
x2, yy2, S2, A2, K2 = ii_mw2(R_conv_carbonI[indLR2]);

# Plot starts here:
f = Figure(resolution=(600,700),  fontsize=13)#, backgroundcolor = :transparent)
ax1 = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Reflected Radiance (I/I₀)",  title="RT Simulations (Multiple-Scattering)")
lines!(ax1,wl, rad_inter(1e7./wl), alpha=0.2)
lines!(ax1,wl_ci, R_conv_carbonI, color = :red)
CairoMakie.xlims!(ax1, 2040, 2378)
CairoMakie.ylims!(ax1, 0, 0.1)
hidexdecorations!(ax1, grid=false)
ax = Axis(f[2,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Derivative wrt trace gas changes",  title="Jacobians")
lines!(ax,wl_ci[indLR2], K2[:,18]*100 , alpha=0.72,color = CarbonI_colors[1], linewidth=2, label="CO₂ * 100")
lines!(ax,wl_ci[indLR2], K2[:,9] , alpha=0.72,color = CarbonI_colors[11], linewidth=2,label="CH₄")
lines!(ax,wl_ci[indLR2], K2[:,27] , alpha=0.72,color = CarbonI_colors[7], linewidth=2,label="N₂O")
lines!(ax,wl_ci[indLR2], K2[:,36]*1000 , alpha=0.72,color = CarbonI_colors[10], linewidth=2,label="H₂O * 1000")
lines!(ax,wl_ci[indLR2], K2[:,45] , alpha=0.72,color = CarbonI_colors[13], linewidth=2,label="CO")

lines!(ax,wl_ci[indLR1], K[:,9]*100 , alpha=0.72,linewidth=2,color = CarbonI_colors[1])
lines!(ax,wl_ci[indLR1], K[:,18] , alpha=0.72,linewidth=2,color = CarbonI_colors[11])
lines!(ax,wl_ci[indLR1], K[:,27] , alpha=0.72,linewidth=2,color = CarbonI_colors[7])
lines!(ax,wl_ci[indLR1], K[:,36]*1000 , alpha=0.72,linewidth=2,color = CarbonI_colors[10])
axislegend(ax, position=:cb, labelsize =12)
#lines!(ax,wl_ci[indLR1], K2[:,50] , alpha=0.72)
#lines!(ax,wl_ci, R_conv_carbonI, color = :red)
CairoMakie.xlims!(ax, 2040, 2378)
#CairoMakie.ylims!(ax1, 0, 0.1)
hidexdecorations!(ax, grid=false)
ax2 = Axis(f[3,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Reflected Radiance (I/I₀)",  title="Non-Scattering retrievals")

lines!(ax2,wl_ci, R_conv_carbonI, color = :red, label="Simulated Measurement")
lines!(ax2, wl_ci[indLR1], yy,color = CarbonI_colors[1], linestyle = :dash, label="Fit")
lines!(ax2, wl_ci[indLR2], yy2,color = CarbonI_colors[1], linestyle = :dash)
lines!(ax2, wl_ci[indLR1], 100(yy .- R_conv_carbonI[indLR1]),color = CarbonI_colors[4], label="Fit - Measurement (*100)")
lines!(ax2, wl_ci[indLR2], 100(yy2 .- R_conv_carbonI[indLR2]),color = CarbonI_colors[4])
axislegend(ax2, position=:cc, labelsize =12)

CairoMakie.xlims!(ax2, 2040, 2378)

f
save("plots/ProxyFits.pdf", f)