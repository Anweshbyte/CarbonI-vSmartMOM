
using CarbonI, vSmartMOM
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials
using CairoMakie


# Load spectroscopies:
co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();

#include(joinpath(@__DIR__, "readSun.jl"))
include(joinpath(@__DIR__, "forwardModel.jl"))

# Define species in the state vector:
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);

# Define wavelength grid for forward model:
Δwl = 0.005
wl = 2035:Δwl:2385

parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")

# Generate Atmospheric Profile
profile_hr = CarbonI.generate_atmos_profile(parameters.T, 100*parameters.p, parameters.q)

# Compute cross sections for all layers:
σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);

# Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
n_layers = 10
profile, σ_matrix = CarbonI.reduce_profile(n_layers,profile_hr, σ_matrix_hr)
#profile = profile_hr
#σ_matrix = σ_matrix_hr

# Use a flat solar spectrum here (these are small in the 2 micron range anyhow):
solarIrr = ones(length(wl));

# Define an instrument:
FWHM  = 1.0  # 
SSI  = 0.68
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
lociBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));
# Generate convolution matrix:
cM = CarbonI.generate_conv_matrix(lociBox,wl, Δwl)
####

# Number of layers:
nL = length(profile.T)
	
# Define prior state vector elements
vmr_co2 = zeros(nL) .+ 400e-6
vmr_ch4 = zeros(nL) .+ 2.0e-6
#vmr_ch4[1:3] .= 1.4e-6
vmr_h2o = zeros(nL) .+0.002 # profile.vcd_h2o ./ profile.vcd_dry
vmr_co  = zeros(nL) .+ 100e-9
vmr_n2o = zeros(nL) .+ 320e-9
#vmr_n2o[1:3] .= 100e-9
vmr_hdo = vmr_h2o * 0.9
vmr_c2h6 = zeros(nL) .+ 1.0e-9
vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, 0.1e-9*vmr_hdo, 0.1e-9*vmr_co2, 0.1e-9vmr_c2h6];

# Define a polynomial scaling for the surface polynomial
p = Legendre([2.0,0.000001,0.00000001,0.0000001,0.000000001]);

# Define our state vector:
x = [vmr_co2; vmr_h2o; vmr_ch4;vmr_co; vmr_n2o; 0.1e-9*vmr_hdo; 0.1e-9*vmr_co2; 0.1e-9*vmr_c2h6;   p[:] ];

# Get SZA from Parameter file:
sza = parameters.sza

# Get prior covariance matrix:
n_state = length(x);
Sₐ = zeros(n_state,n_state);
rel_error = 0.01;
# vcd_ratio = profile_caltech.vcd_dry ./ mean(profile_caltech.vcd_dry)
	
# Fill the diagonal for the trace gases (hard-coded indices, so we have to be careful here):
for i=1:80
	Sₐ[i,i] = (rel_error*x[i])^2   
end
# higher for H2O
for i=11:19
	Sₐ[i,i] = (0.15*x[i])^2   
end
#for i=41:49
#	Sₐ[i,i] = (0.5*x[i])^2   
#end

# Enlarge prior errors near the surface dramatically:
#Sₐ[9,9] = (0.02*x[9])^2
Sₐ[10,10] = (20*x[10])^2
#Sₐ[19,19] = (0.02*x[19])^2
Sₐ[20,20] = (20*x[20])^2
#Sₐ[29,29] = (0.02*x[29])^2
Sₐ[30,30] = (20*x[30])^2
#Sₐ[39,39] = (0.02*x[39])^2
Sₐ[40,40] = (20*x[40])^2
#Sₐ[49,49] = (0.02*x[49])^2
Sₐ[50,50] = (20*x[50])^2
#Sₐ[60,60] = (20*x[60])^2
#Sₐ[70,70] = (20*x[70])^2
#Sₐ[80,80] = (20*x[80])^2
# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
for i=81:n_state
	Sₐ[i,i] = 1e5;
end

# Set the prior state vector:
xa = x;

# Define measurement error covariance matrix:
errors = 1e10*ones(length(lociBox.ν_out));
# Only use a subset:
#errors[1:200] .=1e-5
errors[200:400] .=1e-5
#errors .=1e-5
Se = Diagonal(errors.^2);

#y = CarbonI.conv_spectra(lociBox, wl, spec);


# Start defining column averaging kernels:
h_co2 = zeros(length(x));
h_co2_ = zeros(length(x));
h_ch4 = zeros(length(x));
h_h2o = zeros(length(x));
h_co  = zeros(length(x));
h_hdo = zeros(length(x));
h_n2o = zeros(length(x));
h_c2h6 = zeros(length(x));
ratio = profile.vcd_dry/sum(profile.vcd_dry);
h_co2[1:10] .= ratio;
h_h2o[11:20] .= ratio;
h_ch4[21:30] .= ratio;
h_co[31:40] .= ratio;
h_n2o[41:50] .= ratio;
h_hdo[51:60] .= ratio;
h_co2_[61:70] .= ratio;
h_c2h6[71:80] .= ratio;

N = length(lociBox.ν_out)

# Maximum number of iterations:
max_no_of_iter = 6



# Define measurement here
y = R_conv_carbonI

# Define a polynomial scaling
p = Legendre([maximum(y),0.0000001,0.000000001,0.0000001,0.000000001,0.0000001,0.000000001]);
x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; 0.1e-10*vmr_hdo;0.1e-10*vmr_co2;0.1e-10*vmr_c2h6;   p[:] ];



# Define the state vector for each iteration:
x_all   = zeros((length(x),max_no_of_iter+1))
F_all   = zeros((N,max_no_of_iter))
K = zeros(N,n_state)
A = zeros(n_state,n_state)
result = DiffResults.JacobianResult(zeros(length(lociBox.ν_out)),x);
x_all[:,1]=x
for i=1:max_no_of_iter
    @show i
    #print('Iteration #',i)
	ForwardDiff.jacobian!(result, forward_model_sat_x, x_all[:,i]);
	Kᵢ = DiffResults.jacobian(result);
    Fᵢ = DiffResults.value(result);
    Gain = inv(Kᵢ'inv(Se)Kᵢ + inv(Sₐ))Kᵢ'inv(Se);
    A.=  Gain*Kᵢ
    K .= Kᵢ
    x_all[:,i+1] = xa + Gain * (y - Fᵢ + Kᵢ *(x_all[:,i]-xa))
    println("Column averaged CO2: ", (h_co2' * x_all[:,i+1] * 1e6)/400)
    println("Column averaged N2O: ", (h_n2o' * x_all[:,i+1] * 1e6)/0.32)
    println("Column averaged CH4: ", (h_ch4' * x_all[:,i+1] * 1e9)/2000)
    @show x_all[41:50,i+1]
    #@show h_co2'*x_all[:,end,end]*1e6
    F_all[:,i] = Fᵢ
end

x_final = x_all[:,end]
println("Column averaged CO2: ", h_co2' * x_final * 1e6,"ppm")
println("Column averaged CH4: ", h_ch4' * x_final * 1e9,"ppb")
println("Column averaged N2O: ", h_n2o' * x_final * 1e9,"ppb")
println("Column averaged H2O: ", h_h2o' * x_final * 1e6,"ppm")
println("Column averaged CO: ",   h_co' * x_final * 1e9,"ppb")


plot(lociBox.ν_out, F_all[:,end],label="Carbon-I fit")
plot!(lociBox.ν_out, y, label="Carbon-I spectrum (synthetically generated with vSmartMOM)")

# CO2 column kernel:
Plots.plot(((h_co2'*A)[1,:]./h_co2)[1:10], profile.p, yflip=true, label="CO2 kernel")