
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
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);


# What latitude do we want? Take Caltech as example
myLat = 36.604
myLon = -97.486

myLat = 34.1478
myLon = -118.1445

# Uncomment for high water scenario:
#myLat = 0.0
#myLon = -62
profile_hr = CarbonI.read_atmos_profile_MERRA2(CarbonI.default_merra_file, myLat, myLon, 7);

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
ET  = 44.0u"ms"         # Exposure time
SSI = (2*0.7)u"nm"      # Spectral resolution
Pitch = 18.0u"μm"       # Pixel pitch
FPA_QE = 0.85           # FPA quantum efficiency
Bench_efficiency = 0.6 # Bench efficiency
Fnumber = 2.2           # F-number
readout_noise = 100.0    # Readout noise
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
alb = 0.05
refl = soil(wl)#*0.6; #alb.+0.0*soil(wl)
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
ET  = 44.0u"ms"         # Exposure time
SSI = (2*0.7)u"nm"      # Spectral resolution
Pitch = 18.0u"μm"       # Pixel pitch
FPA_QE = 0.85           # FPA quantum efficiency
Bench_efficiency = 0.63 # Bench efficiency
Fnumber = 2.2           # F-number
readout_noise = 100    # Readout noise
dark_current = 10e3u"1/s" # Dark current

ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);    
nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
nesr_ = nesr./1u"mW/m^2/nm/sr"
e = InstrumentOperator.photons_at_fpa(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
photon_flux =  F/1000 .* lociBox.ν_out * 1e-9/ (6.626e-34 * 2.998e8) 
Se = Diagonal(nesr_.^2);
G = inv(K'inv(Se)K + inv(Sₐ))K'inv(Se);
# Posterior covariance matrix:
Ŝ = inv(K'inv(Se)K + inv(Sₐ));
ch4_error = sqrt(h_ch4' * Ŝ * h_ch4)*1e9
co2_error = sqrt(h_co2' * Ŝ * h_co2)*1e6
h2o_error = sqrt(h_h2o' * Ŝ * h_h2o)*1e6
hdo_error = sqrt(h_hdo' * Ŝ * h_hdo)*1e6
n2o_error = sqrt(h_n2o' * Ŝ * h_n2o)*1e9
co_error  = sqrt(h_co'  * Ŝ * h_co)*1e9
co213_error  = sqrt(h_co213'  * Ŝ * h_co213)*1e6
c2h6_error = sqrt(h_c2h6' * Ŝ * h_c2h6)*1e9

# For co-adding:
@show ch4_error/sqrt(10)
@show co2_error/sqrt(10)
@show n2o_error/sqrt(10)

@show n2o_error/sqrt(12)




albedos = 0.04:0.02:0.8
szas = 5.0:5:75.0

co_error = zeros( length(szas),length(albedos));
ch4_error = similar(co_error);
co2_error = similar(co_error);
n2o_error = similar(co_error);


for (iSza,szai) in enumerate(szas)
	@show szai
    for (iA,albedo) in enumerate(albedos)
		refl = zeros(length(wl)) .+ albedo
		@show sza, albedo
		global sza = szai
		@show sza, albedo
		F = forward_model_x_(x;sza=sza)
		ForwardDiff.jacobian!(result, forward_model_x_, x);
		K = DiffResults.jacobian(result);
		F = DiffResults.value(result);

		nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
		nesr_ = nesr./1u"mW/m^2/nm/sr"
		Se = Diagonal(nesr_.^2);
		G = inv(K'inv(Se)K + inv(Sₐ))K'inv(Se);
		Ŝ = inv(K'inv(Se)K + inv(Sₐ));
		ch4_error[iSza,iA] = sqrt(h_ch4' * Ŝ * h_ch4)*1e9
		co2_error[iSza,iA] = sqrt(h_co2' * Ŝ * h_co2)*1e6
		n2o_error[iSza,iA] = sqrt(h_n2o' * Ŝ * h_n2o)*1e9
		co_error[iSza,iA] = sqrt(h_co'  * Ŝ * h_co)*1e9
    end
end

using NCDatasets

# Define the path to save the NetCDF file
file_path = "errors_as_function_of_albedo_and_sza.nc"

# Create the NetCDF file
ds = Dataset(file_path, "c")

# Define dimensions
defDim(ds, "albedo", length(albedos))
defDim(ds, "sza", length(szas))




# Define a global attribute
ds.attrib["title"] = "Carbon-I precision error estimates (single pixel)"

# Define the variables temperature
alb = defVar(ds, "albedo", Float32, ("albedo",))
sza = defVar(ds, "sza", Float32, ("sza",))
co_error_ = defVar(ds,"CO_error",Float32,("sza","albedo"))
ch4_error_ = defVar(ds,"CH4_error",Float32,("sza","albedo"))
co2_error_ = defVar(ds,"CO2_error",Float32,("sza","albedo"))
n2o_error_ = defVar(ds,"N2O_error",Float32,("sza","albedo"))

# Define and assign coordinate variables
ds["albedo"][:] = albedos
ds["sza"][:] = szas

co_error_.attrib["units"] = "ppb"
ch4_error_.attrib["units"] = "ppb"
co2_error_.attrib["units"] = "ppm"
n2o_error_.attrib["units"] = "ppb"
alb.attrib["units"] = "unitless"
sza.attrib["units"] = "degrees"

co_error_[:] = co_error
co2_error_[:] = co2_error
ch4_error_[:] = ch4_error
n2o_error_[:] = n2o_error	

#close(ds)




# Close the NetCDF file
close(ds)

