
using CarbonI
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials, DelimitedFiles
using CairoMakie
# Load spectroscopies:
co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();

#include(joinpath(@__DIR__, "readSun_DC.jl"))
include(joinpath(@__DIR__, "readSun.jl"))
include(joinpath(@__DIR__, "forwardModel.jl"))

# Load some profile:
MD = "/net/fluo/data1/ftp/XYZT_ESE156/Data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"
MD = CarbonI.default_merra_file
#MD = "./MERRA2_400.tavg3_3d_asm_Np.20200610.nc4"
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);


# What latitude do we want? Take Caltech as example
myLat = 36.604
myLon = -97.486

#myLat = 34.1478
#myLon = -118.1445

myLat = 0.0
myLon = -62
profile_hr = CarbonI.read_atmos_profile_MERRA2(MD, myLat, myLon, 7);

# Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
n_layers = 10

Δwl = 0.01
wl = 2000:Δwl:2400
σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);
nL = length(profile_hr.T)
	
vmr_co2 = zeros(nL) .+ 407e-6
vmr_ch4 = zeros(nL) .+ 1.8e-6
vmr_ch4[1:3] .= 1.4e-6
vmr_h2o = profile_hr.vcd_h2o ./ profile_hr.vcd_dry
vmr_co  = zeros(nL) .+ 100e-9
vmr_n2o = zeros(nL) .+ 337e-9
vmr_n2o[1:3] .= 100e-9
vmr_hdo = vmr_h2o * 0.9
vmr_c2h6 = zeros(nL) .+ 1.0e-9
vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, vmr_hdo, vmr_co2, vmr_c2h6];

sol  = CubicSplineInterpolation(range(wlSol[1],wlSol[end], length=length(wlSol)),solar_irr, extrapolation_bc=Interpolations.Flat());
# Reduce to fewer dimensions:
profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers, profile_hr, σ_matrix_hr,vmrs);



# Define a polynomial scaling
p = Polynomial([0.2,0.0001,0.000001]);

# Define an instrument:
FWHM  = 1.7  # 
SSI  = 0.7
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
lociBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));

# Define state vector:
#x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; vmr_hdo; vmr_co2 ; vmr_c2h6 ;zeros(10) ];
nLeg = 10
xPoly = zeros(nLeg).+eps()
xPoly[1] = 1.0
x = [reduce(vcat,gasProfiles) ; xPoly ];

@show size(x)
sza = 45.0

result = DiffResults.JacobianResult(zeros(length(lociBox.ν_out)),x);


# Define the instrument specs:
ET  = 57.0u"ms"         # Exposure time
SSI = (2*0.7)u"nm"      # Spectral resolution
Pitch = 18.0u"μm"       # Pixel pitch
FPA_QE = 0.85           # FPA quantum efficiency
Bench_efficiency = 0.65 # Bench efficiency
Fnumber = 2.0           # F-number
readout_noise = 70.0    # Readout noise
dark_current = 10e3u"1/s" # Dark current

ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);
clima_alb = readdlm(CarbonI.albedo_file,',', skipstart=1)
#soil = CubicSplineInterpolation(450:2500,r[:,140], extrapolation_bc=Interpolations.Flat());
soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.16, extrapolation_bc=Interpolations.Flat());
solarIrr = sol(wl);
refl   = soil(wl);

L_conv = CarbonI.forward_model_x_(x; sun=sol(wl),reflectance=soil(wl), sza=0.0, instrument=lociBox, profile=profile,σ_matrix=σ_matrix, wl=wl )

nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (L_conv)u"mW/m^2/nm/sr");
nesr_unitless = nesr./1u"mW/m^2/nm/sr";
plot(lociBox.ν_out, L_conv ./ nesr_unitless)
e = InstrumentOperator.photons_at_fpa(ins, (lociBox.ν_out)u"nm", (L_conv)u"mW/m^2/nm/sr");

# Get prior covariance matrix:
n_state = length(x);
Sₐ = zeros(n_state,n_state);
rel_error = 0.0001;
# vcd_ratio = profile_caltech.vcd_dry ./ mean(profile_caltech.vcd_dry)
	
# Fill the diagonal for the trace gases:
for i=1:80
	Sₐ[i,i] = (rel_error*x[i])^2   
end
# CO2 at surface, 100% error
Sₐ[10,10] = (20x[10])^2
Sₐ[20,20] = (20x[20])^2
Sₐ[30,30] = (20x[30])^2
Sₐ[40,40] = (22x[40])^2
Sₐ[50,50] = (22x[50])^2
Sₐ[60,60] = (22x[60])^2
Sₐ[70,70] = (122x[70])^2
Sₐ[80,80] = (1022x[80])^2
# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
for i=81:n_state
	Sₐ[i,i] = 1e2;
end
ratio = profile.vcd_dry/sum(profile.vcd_dry);

h_co2 = zeros(length(x));
h_co213 = zeros(length(x));
h_ch4 = zeros(length(x));
h_h2o = zeros(length(x));
h_co  = zeros(length(x));
h_hdo = zeros(length(x));
h_n2o = zeros(length(x));
h_c2h6 = zeros(length(x));

h_co2[1:10] .= ratio;
h_h2o[11:20] .= ratio;
h_ch4[21:30] .= ratio;
h_co[31:40] .= ratio;
h_n2o[41:50] .= ratio;
h_hdo[51:60] .= ratio;
h_co213[61:70] .= ratio;
h_c2h6[71:80] .= ratio;


# Run the SSA test:

sza = 30
alb = 0.1
refl = alb.+0.0*soil(wl)
ForwardDiff.jacobian!(result, forward_model_x_, x);
K = DiffResults.jacobian(result);
F = DiffResults.value(result);
# Adapt K for the legendre polynomials:
##ranger = range(-1,1,length(F))
#pp = zeros(nLeg)
#for ii=1:nLeg
#    pp = zeros(nLeg)
#    pp[ii] = 1.0
#    p = Legendre(pp)
#    K[:,80+ii] .= p.(ranger)
#end
    
# Define the instrument specs at 400m:
ET  = 40.0u"ms"         # Exposure time
SSI = (2*0.7)u"nm"      # Spectral resolution
Pitch = 18.0u"μm"       # Pixel pitch
FPA_QE = 0.85           # FPA quantum efficiency
Bench_efficiency = 0.65 # Bench efficiency
Fnumber = 2.2           # F-number
readout_noise = 100.0    # Readout noise
dark_current = 10e3u"1/s" # Dark current

ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);    
nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
nesr_ = nesr./1u"mW/m^2/nm/sr"
e = InstrumentOperator.photons_at_fpa(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
photon_flux =  F/1000 .* lociBox.ν_out * 1e-9/ (6.626e-34 * 2.998e8) 

using NCDatasets
f = Dataset("data/Amazon_global.nc")
la = f["data"][:]

f2 = Dataset("data/Amazon_global_oco.nc")
la2 = f2["data"][:]
la2_ = la2[:,end:-1:1]'


laa = la[:,end:-1:1]'

ll = 740
fpa = zeros(3072,512);
fpa[33:end-34,13:end-14] .= laa[ll,:].*e';
laa2 = laa[1400:1400+1680,2750:2750+333]
la_oco2 = la2_[convert(Int,ceil(1400/6.7)):convert(Int,ceil((1400+1680)/6.7)),convert(Int,ceil(2750/6.7)):convert(Int,ceil((2750+333)/6.7))]'
fpa1 = laa2[1,:].*e'
fpa2 = laa2[:,end].*e'
# Create the heatmap

fig = Figure(size=(700, 700))
# Use heatmap, and set axis aspect to `DataAspect()`
ax = Axis(fig[1, 1], aspect = DataAspect())
ax2 = Axis(fig[2, 1], aspect = DataAspect())
#ax3 = Axis(fig[3, 1], aspect = DataAspect())
#ax3 = Axis(fig[2, 2])

# Plot the heatmap with axis limits set to size of the data
CairoMakie.heatmap!(ax2, fpa1, colorrange=(0.0,16100),colormap=:viridis)
#CairoMakie.heatmap!(ax3, fpa2, colorrange=(0.0,18000),colormap=:tokyo)
CairoMakie.heatmap!(ax, laa', colorrange=(0.0,0.6), colormap=:grays)
#CairoMakie.hlines!(ax, [ll], color = :red, linewidth = 2, alpha=0.4)
#lines!(ax3, [zeros(13); photon_flux; zeros(13)], 1:512, color = :black, linewidth = 1)
#hide_axis!(ax)  
#hide_axis!(ax2)
rowsize!(fig.layout,1,Relative(0.828))
hidexdecorations!(ax2, grid=false)
hidexdecorations!(ax, grid=false)
hideydecorations!(ax, grid=false)
hideydecorations!(ax2, grid=false)
#hideydecorations!(ax3, grid=false)
#hidexdecorations!(ax3, grid=false)
rowgap!(fig.layout,5)
#linkyaxes!(ax, ax3)
fig
save("plots/LA_target1.png", fig)

fig = Figure(size=(400,2000))
# Use heatmap, and set axis aspect to `DataAspect()`
ax = Axis(fig[1, 1], aspect = DataAspect())
#ax2 = Axis(fig[2, 1], aspect = DataAspect())
#ax3 = Axis(fig[3, 1], aspect = DataAspect())
#ax3 = Axis(fig[2, 2])

# Plot the heatmap with axis limits set to size of the data
#CairoMakie.heatmap!(ax2, fpa1, colorrange=(0.0,18000),colormap=:tokyo)
#CairoMakie.heatmap!(ax3, fpa2, colorrange=(0.0,18000),colormap=:tokyo)
CairoMakie.heatmap!(ax, laa[1400:1400+1680,2750:2750+333]', colorrange=(0.0,0.6), colormap=:grays)
#CairoMakie.hlines!(ax, [ll], color = :red, linewidth = 2, alpha=0.4)
#lines!(ax3, [zeros(13); photon_flux; zeros(13)], 1:512, color = :black, linewidth = 1)
#hide_axis!(ax)  
#hide_axis!(ax2)
#rowsize!(fig.layout,1,Relative(0.828))
#hidexdecorations!(ax2, grid=false)
hidexdecorations!(ax, grid=false)
hideydecorations!(ax, grid=false)
#hideydecorations!(ax2, grid=false)
#hideydecorations!(ax3, grid=false)
#hidexdecorations!(ax3, grid=false)
#rowgap!(fig.layout,5)
#linkyaxes!(ax, ax3)
fig
save("plots/Amazon_global_map.png", fig)

fig = Figure(size=(400,2000))
# Use heatmap, and set axis aspect to `DataAspect()`
ax = Axis(fig[1, 1], aspect = DataAspect())
#ax2 = Axis(fig[2, 1], aspect = DataAspect())
#ax3 = Axis(fig[3, 1], aspect = DataAspect())
#ax3 = Axis(fig[2, 2])

# Plot the heatmap with axis limits set to size of the data
#CairoMakie.heatmap!(ax2, fpa1, colorrange=(0.0,18000),colormap=:tokyo)
#CairoMakie.heatmap!(ax3, fpa2, colorrange=(0.0,18000),colormap=:tokyo)
CairoMakie.heatmap!(ax, la_oco2[20:27,:], colorrange=(0.0,0.6), colormap=:grays)
#CairoMakie.hlines!(ax, [ll], color = :red, linewidth = 2, alpha=0.4)
#lines!(ax3, [zeros(13); photon_flux; zeros(13)], 1:512, color = :black, linewidth = 1)
#hide_axis!(ax)  
#hide_axis!(ax2)
#rowsize!(fig.layout,1,Relative(0.828))
#hidexdecorations!(ax2, grid=false)
hidexdecorations!(ax, grid=false)
hideydecorations!(ax, grid=false)
#hideydecorations!(ax2, grid=false)
#hideydecorations!(ax3, grid=false)
#hidexdecorations!(ax3, grid=false)
#rowgap!(fig.layout,5)
#linkyaxes!(ax, ax3)
fig
save("plots/Amazon_global_map_OCO2.png", fig)

fig = Figure()
# Use heatmap, and set axis aspect to `DataAspect()`
ax = Axis(fig[1, 1], aspect = DataAspect())
#ax2 = Axis(fig[2, 1], aspect = DataAspect())
#ax3 = Axis(fig[3, 1], aspect = DataAspect())
#ax3 = Axis(fig[2, 2])

# Plot the heatmap with axis limits set to size of the data
#CairoMakie.heatmap!(ax2, fpa1, colorrange=(0.0,18000),colormap=:tokyo)
#CairoMakie.heatmap!(ax3, fpa2, colorrange=(0.0,18000),colormap=:tokyo)
CairoMakie.heatmap!(ax, fpa1, colorrange=(0.0,16100),colormap=:grays)
#CairoMakie.hlines!(ax, [ll], color = :red, linewidth = 2, alpha=0.4)
#lines!(ax3, [zeros(13); photon_flux; zeros(13)], 1:512, color = :black, linewidth = 1)
#hide_axis!(ax)  
#hide_axis!(ax2)
#rowsize!(fig.layout,1,Relative(0.828))
#hidexdecorations!(ax2, grid=false)
hidexdecorations!(ax, grid=false)
hideydecorations!(ax, grid=false)
#hideydecorations!(ax2, grid=false)
#hideydecorations!(ax3, grid=false)
#hidexdecorations!(ax3, grid=false)
#rowgap!(fig.layout,5)
#linkyaxes!(ax, ax3)
fig
save("plots/Amazon_gobal_side1.png", fig)

fig = Figure()
# Use heatmap, and set axis aspect to `DataAspect()`
ax = Axis(fig[1, 1], aspect = DataAspect())
#ax2 = Axis(fig[2, 1], aspect = DataAspect())
#ax3 = Axis(fig[3, 1], aspect = DataAspect())
#ax3 = Axis(fig[2, 2])

# Plot the heatmap with axis limits set to size of the data
#CairoMakie.heatmap!(ax2, fpa1, colorrange=(0.0,18000),colormap=:tokyo)
#CairoMakie.heatmap!(ax3, fpa2, colorrange=(0.0,18000),colormap=:tokyo)
CairoMakie.heatmap!(ax, fpa2, colorrange=(0.0,16100),colormap=:grays)
#CairoMakie.hlines!(ax, [ll], color = :red, linewidth = 2, alpha=0.4)
#lines!(ax3, [zeros(13); photon_flux; zeros(13)], 1:512, color = :black, linewidth = 1)
#hide_axis!(ax)  
#hide_axis!(ax2)
#rowsize!(fig.layout,1,Relative(0.828))
#hidexdecorations!(ax2, grid=false)
hidexdecorations!(ax, grid=false)
hideydecorations!(ax, grid=false)
#hideydecorations!(ax2, grid=false)
#hideydecorations!(ax3, grid=false)
#hidexdecorations!(ax3, grid=false)
#rowgap!(fig.layout,5)
#linkyaxes!(ax, ax3)
fig
save("plots/Amazon_gobal_side2.png", fig)

fig = Figure(size=(300, 500))
ax = Axis(fig[1, 1])
lines!(ax,[zeros(13); e; zeros(13)], 1:512, color = :black, linewidth = 1)
hidedecorations!(ax)
ylims!(1,512)
#hideydecorations!(ax, grid=false)
fig
save("plots/LA_target_spec.pdf", fig)

fig = Figure(size=(500, 80))
ax = Axis(fig[1, 1])
lines!(ax,fpa[:,300], color = :black, linewidth = 1)
hidedecorations!(ax)
xlims!(1,3072)
#hideydecorations!(ax, grid=false)
fig
save("plots/LA_target_space.pdf", fig)

