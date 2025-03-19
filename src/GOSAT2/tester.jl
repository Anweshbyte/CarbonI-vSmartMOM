using ImageFiltering, Interpolations
# Test:
rad, sza, wn, lat, lon = read_gosat2("/home/cfranken/data/GOSAT2/07/06/GOSAT2TFTS220240706151105003_1BSDN00OB1D220221.h5")
Δwl = 0.005
wl = 2035:Δwl:2385
# Define an instrument:
cM_CarbonI, wl_carbonI = CarbonI.create_carbonI_conv_matrix(wl)
cM_EMIT, wl_emit       = CarbonI.create_emit_conv_matrix(wl)

wnn =  wn[1]:mean(diff(wn)):wn[end];
index = 120;
rad_inter = CubicSplineInterpolation(wnn, rad[index,:]);
R_conv_carbonI = cM_CarbonI * rad_inter(1e7./wl);
R_conv_EMIT = cM_EMIT * rad_inter(1e7./wl);

plot(1e7./wnn,rad[index,:]*1e7, alpha=0.7,linewidth=1, label="GOSAT-2")
xlabel!("Wavelength [nm]")
ylabel!("Reflected Radiance")
ylims!(0,3.95)
xlims!(2040,2385)
savefig("GOSAT2_Example1.pdf")

plot(1e7./wnn,rad[index,:]*1e7, alpha=0.3,linewidth=1, label="GOSAT-2")
plot!(wl_emit, R_conv_EMIT*1e7, alpha=0.7, linewidth=2, label="EMIT")
xlabel!("Wavelength [nm]")
ylabel!("Reflected Radiance")
ylims!(0,3.95)
xlims!(2040,2385)
savefig("GOSAT2_Example2.pdf")

plot(1e7./wnn,rad[index,:]*1e7, alpha=0.3,linewidth=1, label="GOSAT-2")
plot!(wl_emit, R_conv_EMIT*1e7, alpha=0.5, linewidth=1.5, label="EMIT")
plot!(wl_carbonI, R_conv_carbonI*1e7, linewidth=2,color=:black, alpha=0.7, label="Carbon-I")
xlabel!("Wavelength [nm]")
ylabel!("Reflected Radiance")
ylims!(0,3.95)
xlims!(2040,2385)
savefig("GOSAT2_Example3.pdf")


plot!(wl_emit, R_conv_EMIT*1e7, label="EMIT")