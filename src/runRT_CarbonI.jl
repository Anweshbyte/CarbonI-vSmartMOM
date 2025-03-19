using vSmartMOM, InstrumentOperator, CarbonI, ImageFiltering, Interpolations, Statistics, CUDA

# Define wavelength grid for forward model:
Δwl = 0.005
wl = 2035:Δwl:2385

# Define an instrument:
cM_CarbonI, wl_carbonI = CarbonI.create_carbonI_conv_matrix(wl)
cM_EMIT, wl_emit    = CarbonI.create_emit_conv_matrix(wl)

# Load parametes and create model
parameters = parameters_from_yaml("src/yaml/carbon-i.yaml")
model = model_from_parameters(parameters);

# # Edit parameters here:
# model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(0.05)
# ref = sum(model.τ_aer[1][1,:]);
# model.τ_aer[1][1,:] .*= 0.1/ref

# # run model
# a = rt_run(model)

# ν = parameters.spec_bands[1];
# gridi = LinRange(ν[1], ν[end], length(ν))

# # Get reflected radiance:
# R = a[1][1,1,:];

# # Interpolate to the grid:
# rad_inter = CubicSplineInterpolation(gridi, R);
# R_wl = rad_inter(1e7./wl);
# # Convolve with instrument:
# R_conv_carbonI = cM_CarbonI * R_wl
# R_conv_emit =    cM_EMIT * R_wl

# # Plot simulated radiances:
# plot(wl, R_wl, label="HighRes", alpha=0.2)
# plot!(wl_carbonI, R_conv_carbonI, label="Carbon-I")
# plot!(wl_emit, R_conv_emit, label="EMIT")

# # Plots ILS:
# plot(wl .- wl_carbonI[100], cM_CarbonI[100,:], label="Carbon-I", alpha=0.5, xlims=(-10,10))
# plot!(wl .- wl_emit[20], cM_EMIT[20,:], label="EMIT", alpha=0.5, xlims=(-10,10))

# plot(1e7./gridi, compute_absorption_cross_section(a[1], gridi, 1013, 298.0),alpha=0.5, label=parameters.absorption_params.molecules[1][1], yscale=:log10)
# plot!(1e7./gridi, compute_absorption_cross_section(a[2], gridi, 1013, 298.0),alpha=0.5, label=parameters.absorption_params.molecules[1][2], yscale=:log10)
# plot!(1e7./gridi, compute_absorption_cross_section(a[3], gridi, 1013, 298.0),alpha=0.5, label=parameters.absorption_params.molecules[1][3], yscale=:log10)
# plot!(1e7./gridi, compute_absorption_cross_section(a[4], gridi, 1013, 298.0),alpha=0.5, label=parameters.absorption_params.molecules[1][4], yscale=:log10)

# plot(1e7./gridi, compute_absorption_cross_section(a[5], gridi, 1013, 298.0),alpha=0.5, label=parameters.absorption_params.molecules[1][5], yscale=:log10)
# plot!(1e7./gridi, compute_absorption_cross_section(a[6], gridi, 1013, 298.0),alpha=0.5, label=parameters.absorption_params.molecules[1][6], yscale=:log10)
# plot!(1e7./gridi, compute_absorption_cross_section(a[7], gridi, 1013, 298.0),alpha=0.25, label=parameters.absorption_params.molecules[1][7], yscale=:log10)
# plot!(1e7./gridi, compute_absorption_cross_section(a[8], gridi, 1013, 298.0),alpha=0.5, label=parameters.absorption_params.molecules[1][8], yscale=:log10)