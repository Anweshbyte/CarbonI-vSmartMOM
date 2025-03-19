using CarbonI
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials, DelimitedFiles
using CairoMakie
using Artifacts, LazyArtifacts
# Load spectroscopies:

co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();

#include(joinpath(@__DIR__, "readSun_DC.jl"))
include(joinpath(@__DIR__, "src/readSun.jl"))
include(joinpath(@__DIR__, "src/forwardModel.jl"))

# Load some profile:
MD = CarbonI.default_merra_file
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);

# What latitude do we want? Take Amazon with lots of H2O for that requirement
myLat = 0.0
myLon = -62
# Read in high resolution profile:
profile_hr = CarbonI.read_atmos_profile_MERRA2(MD, myLat, myLon, 7);

# Reduce dimensions, group layers together to get  layers of equal pressure difference:
n_layers = 10

# Define wavelength grid:
Δwl = 0.01
wl = 2000:Δwl:2400
# Precompute the cross sections:
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

# Reduce profile and cross sections to fewer dimensions:
profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers, profile_hr, σ_matrix_hr,vmrs);


# Define an instrument:
# LSF(λ)
FWHM  = 2.2  # FWHM  = 2.2 creates an effective FWHM of 2.5nm, i.e. required. (1.5 creates an effective kernel of 1.88)
# Spectral Sampling Interval
SSI  = 0.7
# Slit blur (2*SSI box)
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
# LSF (Gaussian)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
# Pixel response (1*SSI)
kern3 = CarbonI.box_kernel(SSI, Δwl)

# Combine the kernels:
kernf = imfilter(imfilter(kern1, kern2), kern3)
# Create the instrument response:
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
ET  = 44.0u"ms"          # Exposure time
SSI = (2*0.7)u"nm"       # Spectral resolution
Pitch = 18.0u"μm"        # Pixel pitch
FPA_QE = 0.8             # FPA quantum efficiency
Bench_efficiency = 0.6*1.3   # Bench efficiency (0.6 Required, 0.6*1.5 goal, include QE bump)
Fnumber = 2.2            # F-number
readout_noise = 100.0    # Readout noise (100 CBE, 120 required
dark_current = 5e3u"1/s" # Dark current

ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);
clima_alb = readdlm(CarbonI.albedo_file,',', skipstart=1)
#soil = CubicSplineInterpolation(450:2500,r[:,140], extrapolation_bc=Interpolations.Flat());
soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.16, extrapolation_bc=Interpolations.Flat());
solarIrr = sol(wl);
refl   = soil(wl);

L_conv = CarbonI.forward_model_x_(x; sun=sol(wl),reflectance=refl, sza=30.0, instrument=lociBox, profile=profile,σ_matrix=σ_matrix, wl=wl )

nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (L_conv)u"mW/m^2/nm/sr");
nesr_unitless = nesr./1u"mW/m^2/nm/sr";
#plot(lociBox.ν_out, L_conv ./ nesr_unitless)
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
Sₐ[10,10] = (200x[10])^2
Sₐ[20,20] = (200x[20])^2
Sₐ[30,30] = (200x[30])^2
Sₐ[40,40] = (200x[40])^2
Sₐ[50,50] = (200x[50])^2
Sₐ[60,60] = (200x[60])^2
Sₐ[70,70] = (200x[70])^2
Sₐ[80,80] = (20000x[80])^2
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
refl = soil(wl)#*0.6; #alb.+0.0*soil(wl)
ForwardDiff.jacobian!(result, forward_model_x_, x);
K = DiffResults.jacobian(result);
F = DiffResults.value(result);


ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);    
nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
nesr_ = nesr./1u"mW/m^2/nm/sr"
e = InstrumentOperator.photons_at_fpa(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
photon_flux =  F/1000 .* lociBox.ν_out * 1e-9/ (6.626e-34 * 2.998e8) 
Se = Diagonal(nesr_.^2);
G = inv(K'inv(Se)K + inv(Sₐ))K'inv(Se);

# Posterior covariance matrix (Profile retrievals, 10 layers):
Ŝ = inv(K'inv(Se)K + inv(Sₐ));

# Apply column operator to get to total column uncertainties:
ch4_error    = sqrt(h_ch4' * Ŝ * h_ch4)*1e9
co2_error    = sqrt(h_co2' * Ŝ * h_co2)*1e6
h2o_error    = sqrt(h_h2o' * Ŝ * h_h2o)*1e6
hdo_error    = sqrt(h_hdo' * Ŝ * h_hdo)*1e6
n2o_error    = sqrt(h_n2o' * Ŝ * h_n2o)*1e9
co_error     = sqrt(h_co'  * Ŝ * h_co)*1e9
co213_error  = sqrt(h_co213'  * Ŝ * h_co213)*1e6
c2h6_error   = sqrt(h_c2h6' * Ŝ * h_c2h6)*1e9

# For co-adding:
@show ch4_error/sqrt(10)
@show co2_error/sqrt(10)
@show n2o_error/sqrt(10)

@show n2o_error/sqrt(12)


FWHMs = 1.5:0.1:2.8
throughputs = 0.4:0.02:0.8
readout_noises = 80:5:125
n2o_errors_relative = zeros(length(FWHMs),length(throughputs), length(readout_noises))
SSI = (2*0.7)u"nm" 

for (i,FWHM) in enumerate(FWHMs)
	SSI_  = 0.7
	# Slit blur (2*SSI box)
	#kern1 = CarbonI.box_kernel(2*SSI_, Δwl)
	# LSF (Gaussian)
	kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
	# Pixel response (1*SSI)
	#kern3 = CarbonI.box_kernel(SSI_, Δwl)

	# Combine the kernels:
	kernf = imfilter(imfilter(kern1, kern2), kern3)
	# Create the instrument response:
	lociBox = CarbonI.KernelInstrument(kern2, collect(2040:SSI_:2380));
	ForwardDiff.jacobian!(result, forward_model_x_, x);
	K = DiffResults.jacobian(result);
	F = DiffResults.value(result);
	for (j,throughput) in enumerate(throughputs)
		for (k,readout_noise) in enumerate(readout_noises)
			FPA_QE = 1.0             # FPA quantum efficiency
			Bench_efficiency = throughput   # Bench efficiency (0.6 Required, 0.6*1.5 goal, include QE bump)
			
			#readout_noise = 120.0
			ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);
			nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
			nesr_ = nesr./1u"mW/m^2/nm/sr"
			e = InstrumentOperator.photons_at_fpa(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
			photon_flux =  F/1000 .* lociBox.ν_out * 1e-9/ (6.626e-34 * 2.998e8) 
			Se = Diagonal(nesr_.^2);
			G = inv(K'inv(Se)K + inv(Sₐ))K'inv(Se);

			# Posterior covariance matrix (Profile retrievals, 10 layers):
			Ŝ = inv(K'inv(Se)K + inv(Sₐ));

			# Apply column operator to get to total column uncertainties:
			
			n2o_error    = sqrt(h_n2o' * Ŝ * h_n2o)*1e9
			n2o_errors_relative[i,j,k] = n2o_error/sqrt(12)/337*100
			@show n2o_error/sqrt(12)/337*100
		end
	end
end

f = Figure(resolution=(500,400), title="Precision Sensitivity Analysis", fontsize=16)
ax1 = Axis(f[1,1], ylabel="Total System Throughput",xlabel="FWHM (nm)",  title="N₂O Precision Error (%)")
contourf!(ax1,FWHMs, throughputs, n2o_errors_relative[:,:,5], labels=true, levels=2.0:0.125:4, colorrange=(2.4,4), colormap = (:viridis, 0.9));
contour!(ax1, FWHMs, throughputs, n2o_errors_relative[:,:,5], levels=2.5:0.25:3.75, labels=true,  labelsize = 14,labelfont = :bold, labelcolor = :white);
contour!(ax1, FWHMs, throughputs, n2o_errors_relative[:,:,5], levels=[3.5], labels=true,  labelsize = 14,labelfont = :bold, labelcolor = :red);
CairoMakie.xlims!(FWHMs[1],FWHMs[end])
CairoMakie.ylims!(throughputs[1],throughputs[end])
x1 = 2.5
y1 = 0.48
scatter!(ax1,x1, y1)
text!(x1, y1, text = "Requirement", align = (:left,:bottom))

x2 = 1.7
y2 = 0.61
scatter!(ax1,x2, y2)
text!(x2, y2, text = "CBE", align = (:left,:bottom))

f
save("plots/PrecisionSensitivity_study_v1.pdf", f)

f = Figure(resolution=(500,400), title="Precision Sensitivity Analysis", fontsize=16)
ax1 = Axis(f[1,1], xlabel="Total System Throughput",ylabel="Readout noise (e⁻¹)",  title="N₂O Precision Error (%)")
contourf!(ax1,throughputs, readout_noises,n2o_errors_relative[4,:,:], labels=true, levels=2.0:0.125:4, colorrange=(2.4,4), colormap = (:viridis, 0.9));
contour!(ax1, throughputs, readout_noises,n2o_errors_relative[4,:,:], levels=2.5:0.25:3.75, labels=true,  labelsize = 14,labelfont = :bold, labelcolor = :white);
contour!(ax1, throughputs, readout_noises,n2o_errors_relative[4,:,:], levels=[3.5], labels=true,  labelsize = 14,labelfont = :bold, labelcolor = :red);
CairoMakie.ylims!(readout_noises[1],readout_noises[end])
CairoMakie.xlims!(throughputs[1],throughputs[end])
x1 = 0.48
y1 = 120
scatter!(ax1,x1, y1)
text!(x1, y1, text = "Requirement", align = (:left,:bottom))

x2 = 0.61
y2 = 100
scatter!(ax1,x2, y2)
text!(x2, y2, text = "CBE", align = (:right,:top))

f
save("plots/PrecisionSensitivity_study_v2.pdf", f)
