using CairoMakie
using JLD2, Statistics, vSmartMOM, InstrumentOperator, CarbonI 
using ImageFiltering, Interpolations, Statistics, CUDA, Plots, Distributions

# Define some colors and themes
set_theme!(CairoMakie.theme_black())
# Define the color using the hex code
my_color = RGB(0xDC / 255, 0xCA / 255, 0xB0 / 255)  # Normalize the hex values to [0, 1]


# Define our standard fitting windows:
include("./src/ProxyTests/setupTestFits_generic.jl")
include("./src/ProxyTests/createMicroWindows.jl")
include("./src/Plots/CI_colors.jl")

# Wavelength grid for high resolution simulations:
Δwl = 0.005
wl = 2035:Δwl:2385

# Define an instrument operator (for Carbon-I and EMIT):
cM, wl_ci = CarbonI.create_carbonI_conv_matrix_cbe(wl)
cM_emit, wl_emit = CarbonI.create_emit_conv_matrix(wl)

# Load default scenario and run the RT model (the fit has the same input, T, H2O and other profiles could be shuffled here to create inconsistencies)
parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")
aod = 0.1
p_aero = 900.0
albedo = 0.1
sza = 30.0

# Update parameters with the variables set above:
parameters.scattering_params.rt_aerosols[1].profile = Normal(p_aero, 50.0)
parameters.sza = sza

# Generate the model from the parameters:
model = model_from_parameters(parameters);
aod_profile_normalized = model.τ_aer[1][1,:] / sum(model.τ_aer[1][1,:])
model.τ_aer[1][1,:] = aod * aod_profile_normalized
# Lambertian albedo, this could be changed to a polynomial or a more complex model as well
model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(albedo)

# Run the model:
a = rt_run(model)

# Generate the convolved radiances for EMIT and Carbon-I 
ν = parameters.spec_bands[1];
Δν = mean(diff(ν));
gridi = ν[1]:Δν:ν[end]+10eps()
R = a[1][1,1,:];
rad_inter = CubicSplineInterpolation(gridi, R);
R_conv_carbonI = cM*rad_inter(1e7./wl);
R_conv_EMIT = cM_emit*rad_inter(1e7./wl);

# Perform fits in Microwindows:
x, yy, S, A, K = ii_mw1(R_conv_carbonI[indLR1]);
x2, yy2, S2, A2, K2 = ii_mw2(R_conv_carbonI[indLR2]);

# Plot starts here:
f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax1 = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Reflected Radiance (I/I₀)",  title="RT Simulations (Multiple-Scattering)")
lines!(ax1,wl, rad_inter(1e7./wl), color=my_color, linewidth=0.2)
CairoMakie.xlims!(ax1, 2040, 2378)
CairoMakie.ylims!(ax1, 0, 0.1)
save("plots/Fits_1.pdf", f)
#lines!(ax1,wl_ci, R_conv_carbonI, color = :orange,linewidth=3)
#save("plots/Fits_2.pdf", f)


f

f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax1 = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Reflected Radiance (I/I₀)",  title="RT Simulations (Multiple-Scattering)")
lines!(ax1,wl, rad_inter(1e7./wl), color=my_color, linewidth=0.2, alpha=0.4)
CairoMakie.xlims!(ax1, 2040, 2378)
CairoMakie.ylims!(ax1, 0, 0.1)

lines!(ax1,wl_ci, R_conv_carbonI, color = :orange,linewidth=3)
save("plots/Fits_2.pdf", f)


f

f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax1 = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Reflected Radiance (I/I₀)",  title="RT Simulations (Multiple-Scattering)")
lines!(ax1,wl, rad_inter(1e7./wl), color=my_color, linewidth=0.2, alpha=0.4)
CairoMakie.xlims!(ax1, 2040, 2378)
CairoMakie.ylims!(ax1, 0, 0.1)

lines!(ax1,wl_emit, R_conv_EMIT, color = :red,linewidth=3, alpha=0.5)
save("plots/Fits_3.pdf", f)
lines!(ax1,wl_ci, R_conv_carbonI, color = :orange,linewidth=2.5, alpha=0.75)
save("plots/Fits_4.pdf", f)


f

# Jacobians:
f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Derivative wrt trace gas changes",  title="Jacobians")
#lines!(ax,wl_ci[indLR2], K2[:,18]*100 , alpha=0.72,color = CarbonI_colors[1], linewidth=2, label="CO₂ * 100")
lines!(ax,wl_ci[indLR2], K2[:,9] , alpha=0.72,color = CarbonI_colors[11], linewidth=2,label="CH₄")
#lines!(ax,wl_ci[indLR2], K2[:,27] , alpha=0.72,color = CarbonI_colors[7], linewidth=2,label="N₂O")
#lines!(ax,wl_ci[indLR2], K2[:,36]*1000 , alpha=0.72,color = CarbonI_colors[10], linewidth=2,label="H₂O * 1000")
#lines!(ax,wl_ci[indLR2], K2[:,45] , alpha=0.72,color = CarbonI_colors[13], linewidth=2,label="CO")

#lines!(ax,wl_ci[indLR1], K[:,9]*100 , alpha=0.72,linewidth=2,color = CarbonI_colors[1])
lines!(ax,wl_ci[indLR1], K[:,18] , alpha=0.72,linewidth=2,color = CarbonI_colors[11])
#lines!(ax,wl_ci[indLR1], K[:,27] , alpha=0.72,linewidth=2,color = CarbonI_colors[7])
#lines!(ax,wl_ci[indLR1], K[:,36]*1000 , alpha=0.72,linewidth=2,color = CarbonI_colors[10])
CairoMakie.xlims!(ax, wl_ci[1], wl_ci[end])
axislegend(ax, position=:cb, labelsize =12)
save("plots/Fits_Jac_CH4.pdf", f)

f

# Jacobians:
f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Derivative wrt trace gas changes",  title="Jacobians")
#lines!(ax,wl_ci[indLR2], K2[:,18]*100 , alpha=0.72,color = CarbonI_colors[1], linewidth=2, label="CO₂ * 100")
lines!(ax,wl_ci[indLR2], K2[:,9] , alpha=0.72,color = CarbonI_colors[11], linewidth=2,label="CH₄")
#lines!(ax,wl_ci[indLR2], K2[:,27] , alpha=0.72,color = CarbonI_colors[7], linewidth=2,label="N₂O")
#lines!(ax,wl_ci[indLR2], K2[:,36]*1000 , alpha=0.72,color = CarbonI_colors[10], linewidth=2,label="H₂O * 1000")
#lines!(ax,wl_ci[indLR2], K2[:,45] , alpha=0.72,color = CarbonI_colors[13], linewidth=2,label="CO")

#lines!(ax,wl_ci[indLR1], K[:,9]*100 , alpha=0.72,linewidth=2,color = CarbonI_colors[1])
lines!(ax,wl_ci[indLR1], K[:,18] , alpha=0.72,linewidth=2,color = CarbonI_colors[11])
#lines!(ax,wl_ci[indLR1], K[:,27] , alpha=0.72,linewidth=2,color = CarbonI_colors[7])
#lines!(ax,wl_ci[indLR1], K[:,36]*1000 , alpha=0.72,linewidth=2,color = CarbonI_colors[10])
CairoMakie.xlims!(ax, wl_ci[1], wl_ci[end])
#axislegend(ax, position=:cb, labelsize =12)
#save("plots/Fits_Jac_CH4.pdf", f)
lines!(ax,wl_ci[indLR2], K2[:,18]*100 , alpha=0.72,color = CarbonI_colors[1], linewidth=2, label="CO₂ * 100")
lines!(ax,wl_ci[indLR1], K[:,9]*100 , alpha=0.72,linewidth=2,color = CarbonI_colors[1])
axislegend(ax, position=:cb, labelsize =12)
save("plots/Fits_Jac_CH4_CO2.pdf", f)
f

# Jacobians:
f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Derivative wrt trace gas changes",  title="Jacobians")
lines!(ax,wl_ci[indLR2], K2[:,18]*100 , alpha=0.72,color = CarbonI_colors[1], linewidth=2, label="CO₂ * 100")
#lines!(ax,wl_ci[indLR2], K2[:,9] , alpha=0.72,color = CarbonI_colors[11], linewidth=2,label="CH₄")
#lines!(ax,wl_ci[indLR2], K2[:,27] , alpha=0.72,color = CarbonI_colors[7], linewidth=2,label="N₂O")
#lines!(ax,wl_ci[indLR2], K2[:,36]*1000 , alpha=0.72,color = CarbonI_colors[10], linewidth=2,label="H₂O * 1000")
lines!(ax,wl_ci[indLR2], K2[:,45] , alpha=0.72,color = CarbonI_colors[13], linewidth=3,label="CO")

lines!(ax,wl_ci[indLR1], K[:,9]*100 , alpha=0.72,linewidth=2,color = CarbonI_colors[1])
#lines!(ax,wl_ci[indLR1], K[:,18] , alpha=0.72,linewidth=2,color = CarbonI_colors[11])
#lines!(ax,wl_ci[indLR1], K[:,27] , alpha=0.72,linewidth=2,color = CarbonI_colors[7])
#lines!(ax,wl_ci[indLR1], K[:,36]*1000 , alpha=0.72,linewidth=2,color = CarbonI_colors[10])
CairoMakie.xlims!(ax, wl_ci[1], wl_ci[end])
#axislegend(ax, position=:cb, labelsize =12)
#save("plots/Fits_Jac_CH4.pdf", f)

axislegend(ax, position=:cb, labelsize =12)
save("plots/Fits_Jac_CO2_CO.pdf", f)
f

# Jacobians:
f = Figure(resolution=(600,400),  fontsize=13, backgroundcolor = :black, fontcolor=:white)
ax = Axis(f[1,1],yminorgridvisible = true,  yminorticks = IntervalsBetween(5), xlabel="Wavelength (nm)",ylabel="Derivative wrt trace gas changes",  title="Jacobians")
lines!(ax,wl_ci[indLR2], K2[:,18]*100 , alpha=0.72,color = CarbonI_colors[1], linewidth=2, label="CO₂ * 100")
#lines!(ax,wl_ci[indLR2], K2[:,9] , alpha=0.72,color = CarbonI_colors[11], linewidth=2,label="CH₄")
lines!(ax,wl_ci[indLR2], K2[:,27] , alpha=0.72,color = CarbonI_colors[7], linewidth=2,label="N₂O")
#lines!(ax,wl_ci[indLR2], K2[:,36]*1000 , alpha=0.72,color = CarbonI_colors[10], linewidth=2,label="H₂O * 1000")
lines!(ax,wl_ci[indLR2], K2[:,45] , alpha=0.72,color = CarbonI_colors[13], linewidth=3,label="CO")

lines!(ax,wl_ci[indLR1], K[:,9]*100 , alpha=0.72,linewidth=2,color = CarbonI_colors[1])
#lines!(ax,wl_ci[indLR1], K[:,18] , alpha=0.72,linewidth=2,color = CarbonI_colors[11])
lines!(ax,wl_ci[indLR1], K[:,27] , alpha=0.72,linewidth=2,color = CarbonI_colors[7])
#lines!(ax,wl_ci[indLR1], K[:,36]*1000 , alpha=0.72,linewidth=2,color = CarbonI_colors[10])
CairoMakie.xlims!(ax, wl_ci[1], wl_ci[end])
#axislegend(ax, position=:cb, labelsize =12)
#save("plots/Fits_Jac_CH4.pdf", f)

axislegend(ax, position=:cb, labelsize =12)
save("plots/Fits_Jac_CO2_CO_N2O.pdf", f)
f

