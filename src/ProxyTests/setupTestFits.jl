using CarbonI, vSmartMOM
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials, Plots
#using CairoMakie

# Load spectroscopies:
co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();

# Read Solar Spectra (not used for the simple tests)
#include(joinpath(@__DIR__, "readSun.jl"))
indHR = 1:70000
indLR = 1:150
# Read non-scattering forward model
include("/home/cfranken/code/gitHub/CarbonI/src/forwardModel.jl")

# Define species in the state vector:
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);

# Define wavelength grid for forward model:
Δwl = 0.005
wl = 2031:Δwl:2385

parameters = parameters_from_yaml("/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml")
ν = parameters.spec_bands[1];
gridi = LinRange(ν[1], ν[end], length(ν))

# Generate Atmospheric Profile
profile_hr = CarbonI.generate_atmos_profile(parameters.T, 100*parameters.p, parameters.q)

# Compute cross sections for all layers:
σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);

# Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
n_layers = 10
profile, σ_matrix, indis = CarbonI.reduce_profile(n_layers, profile_hr, σ_matrix_hr)
n_layers = length(indis)

# Use a flat solar spectrum here (these are small in the 2 micron range anyhow):
solarIrr = ones(length(wl));

# Define an instrument:
cM, wl_ci = CarbonI.create_carbonI_conv_matrix(wl)

# Number of layers:
nL = length(profile.T)
h2o = [  5e-6,   5e-6,   5e-6,   5e-6,   5e-61,   5e-6,   
5e-6,   5e-6,   5e-6,   5e-6,  5e-6,   5e-6,  
5e-6,   5e-6,   5e-6,  5e-6,  5e-6,  5e-6,  
5e-6,  5e-6,  5e-6,  5e-6,  5e-6, 5e-6,
5e-6, 5e-6, 0.001, 0.0025, 0.005, 0.007, 
0.01, 0.015, 0.02, 0.03]
# Define prior state vector elements
vmr_co2 = zeros(nL) .+ 400e-6
vmr_ch4 = zeros(nL) .+ 2.0e-6
#vmr_ch4[1:3] .= 1.4e-6
vmr_h2o = [(h2o[ind]' * profile_hr.vcd_dry[ind])/sum(profile_hr.vcd_dry[ind]) for ind in indis] # profile.vcd_h2o ./ profile.vcd_dry
vmr_co  = zeros(nL) .+ 100e-9
vmr_n2o = zeros(nL) .+ 320e-9
#vmr_n2o[1:3] .= 100e-9
vmr_hdo = vmr_h2o
vmr_c2h6 = zeros(nL) .+ 1.0e-9
vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, 0.1e-9*vmr_hdo, 0.1e-9*vmr_co2, 0.1e-9vmr_c2h6];

# Define a polynomial scaling for the surface polynomial
p = Legendre([2.0,0.000001,0.00000001,0.0000001,0.000000001,0.0000001,0.000000001]);

# Define our state vector:
x = [vmr_co2; vmr_h2o; vmr_ch4;vmr_co; vmr_n2o; vmr_hdo; vmr_co2; vmr_c2h6;   p[:] ];

# Get SZA from Parameter file:
sza = parameters.sza

# Get prior covariance matrix:
n_state = length(x);
Sₐ = zeros(n_state,n_state);

rel_error = 0.25;
correlationMatrix = CarbonI.createCorrelationMatrix(profile.p/1e2,0.0,1550)
correlationMatrix_h2o = CarbonI.createCorrelationMatrix(profile.p/1e2,0.0,550)

# Create Diagonal errors for the trace gases:
e_ch4  = CarbonI.createErrorVector(profile.p/1e2,250.0,850.0, rel_error, vmr_ch4)
e_co2  = CarbonI.createErrorVector(profile.p/1e2,250.0,850.0, rel_error, vmr_co2)
#e_co2[end] = vmr_co2[end]*20
e_h2o  = CarbonI.createErrorVector(profile.p/1e2,250.0,850.0, 2*rel_error, vmr_h2o)
e_n2o  = CarbonI.createErrorVector(profile.p/1e2,250.0,850.0, rel_error, vmr_n2o)
e_co   = CarbonI.createErrorVector(profile.p/1e2,250.0,850.0, rel_error, vmr_co)
e_hdo  = CarbonI.createErrorVector(profile.p/1e2,250.0,850.0, 2*rel_error, vmr_hdo)
e_c2h6 = CarbonI.createErrorVector(profile.p/1e2,250.0,850.0, rel_error, vmr_c2h6)

dims = size(σ_matrix)	
# Fill the diagonal for the trace gases (hard-coded indices, so we have to be careful here):
e_coll = (e_co2, e_h2o, e_ch4, e_co, e_n2o, e_hdo, e_co2, e_c2h6)
index = 1
for i=1:n_layers:n_layers*8
    #@show i,i+n_layers
    if i==n_layers+1
        Sₐ[i:i+n_layers-1,i:i+n_layers-1] .= CarbonI.createPriorCovarianceMatrix(e_coll[index], correlationMatrix_h2o)
    else
        Sₐ[i:i+n_layers-1,i:i+n_layers-1] .= CarbonI.createPriorCovarianceMatrix(e_coll[index], correlationMatrix)
    end
    index += 1
end


# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
for i=8n_layers+1:n_state
	Sₐ[i,i] = 1e5;
end

# Set the prior state vector:
xa = x;

# Define measurement error covariance matrix:
errors = 1e10*ones(length(wl_ci));
# Only use a subset:
#errors[88:151] .= 0.00002
#errors[300:400] .= 0.0001
#errors[1:480] .= 0.00001
ind = 1:162
#errors[ind] .= 0.00002
errors[500:502] .= 0.00002
errors[257:410] .= 0.00002
errors[1:2] .= 0.00002
#errors[304:380] .= 0.00002
#errors[1:480] .= 0.0001
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
h_co2[1:n_layers] .= ratio;
h_h2o[1n_layers+1:2n_layers] .= ratio;
h_ch4[2n_layers+1:3n_layers] .= ratio;
h_co[3n_layers+1:4n_layers] .= ratio;
h_n2o[4n_layers+1:5n_layers] .= ratio;
h_hdo[5n_layers+1:6n_layers] .= ratio;
h_co2_[6n_layers+1:7n_layers] .= ratio;
h_c2h6[7n_layers+1:8n_layers] .= ratio;

N = length(wl_ci)



# Define a polynomial scaling
#p = Legendre([maximum(y),0.0000001,0.000000001,0.0000001,0.000000001,0.0000001,0.000000001]);
p = Legendre([0.05,0.0000001,0.000000001,0.0000001,0.000000001,0.0000001,0.000000001]);
x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; vmr_hdo;vmr_co2;vmr_c2h6;   p[:] ];

n_state = length(x);

# Define the state vector for each iteration:
max_no_of_iter = 4
x_all   = zeros((length(x),max_no_of_iter+1)) #State vector and max iters
F_all   = zeros((N,max_no_of_iter)) #Wavelengths and max iters
K = zeros(N,n_state)
A = zeros(n_state,n_state)
result = DiffResults.JacobianResult(zeros(length(wl_ci)),x);

s_aod = model.τ_aer[1][1,:]
s_aod = s_aod/sum(s_aod)


x_all[:,1]=x
n2o = []
co2 = []
co2_13 = []
ch4 = []
h2o = []
resi = []
resi2 = []
ys = []


#for i=1:200
for aodd = 0.00:0.003:0.1   
    for albedo = 0.03:0.01:0.35
 
    ##model.τ_aer[1][1,:] .= aodd .* s_aod
    #@show sum(model.τ_aer[1][1,:])
    model.τ_aer[1][1,:] .= aodd .* s_aod
    model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(albedo)
    #model.params.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(0.2)
    a = rt_run(model)
    # Get reflected radiance:   
    R = a[1][1,1,:];
    # Interpolate to the grid:
    rad_inter = CubicSplineInterpolation(gridi, R);
    R_conv_carbonI = cM*rad_inter(1e7./wl);
    x[73] = maximum(R_conv_carbonI)
    x_all[:,1]=x
# Define measurement here
y = R_conv_carbonI;# + 0.0002*randn(N)
push!(ys, y)



res = similar(y)
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
    println("Column averaged CO₂: ", (h_co2' * x_all[:,i+1] * 1e6)/400*100)
    println("Column averaged 13CO₂: ", (h_co2_' * x_all[:,i+1] * 1e6)/400*100)
    println("Column averaged N₂O: ", (h_n2o' * x_all[:,i+1] * 1e6)/0.32*100)
    println("Column averaged CH₄: ", (h_ch4' * x_all[:,i+1] * 1e9)/2000*100)
    println("Column averaged C2H6: ", (h_c2h6' * x_all[:,i+1] * 1e9)/1.0*100)
    println("Column averaged H2O: ", (h_h2o' * x_all[:,i+1] * 1e6)/5313.10461*100)
    println("Column averaged CO: ", (h_co' * x_all[:,i+1] * 1e6)/0.1*100)
    
    #@show x_all[41:50,i+1]
    res .= y-Fᵢ
    @show 1e6*x_all[1:9,i+1]
    @show 1e9*x_all[4n_layers+1:5n_layers,i+1]
    @show x_all[8n_layers+1:n_state,i+1]
    #@show res[ind]' * res[ind]
    #@show h_co2'*x_all[:,end,end]*1e6
    F_all[:,i] = Fᵢ
end
    push!(n2o,(h_n2o' * x_all[:,end] * 1e6)/0.32)
    push!(co2,(h_co2' * x_all[:,end] * 1e6)/400)
    push!(ch4,(h_ch4' * x_all[:,end] * 1e9)/2000)
    push!(h2o,(h_h2o' * x_all[:,end] * 1e6)/5313.104611587684)
    push!(co2_13, (h_co2_' * x_all[:,end] * 1e6)/400)
    push!(resi, sqrt(res[ind]' * res[ind]))
    push!(resi2, res[ind])
end
end
plot(albedo, n2o, label="N2O")
plot!(albedo, ch4, label="CH4")
plot!(albedo, co2, label="CO2")
xlabel!("Albedo")
ylabel!("Column averaged VMR (retrieved/true)")
title!("Non-scattering retrieval, AOD=0.0 @ 2.3µm")

x_final = x_all[:,end]
println("Column averaged CO2: ", h_co2' * x_final * 1e6,"ppm")
println("Column averaged CH4: ", h_ch4' * x_final * 1e9,"ppb")
println("Column averaged N2O: ", h_n2o' * x_final * 1e9,"ppb")
println("Column averaged H2O: ", h_h2o' * x_final * 1e6,"ppm")
println("Column averaged CO: ",   h_co' * x_final * 1e9,"ppb")


plot(lociBox.ν_out, F_all[:,end],label="Carbon-I fit")
plot!(lociBox.ν_out, y, label="Carbon-I spectrum (synthetically generated with vSmartMOM)")

# CO2 column kernel:
Plots.plot(((h_co2'*A)[1,:]./h_co2)[1:9], profile.p, yflip=true, label="CO2 kernel")
Plots.plot!(((h_ch4'*A)[1,:]./h_ch4)[19:19+8], profile.p, yflip=true, label="CH4 kernel")
Plots.plot!(((h_n2o'*A)[1,:]./h_n2o)[4n_layers+1:5n_layers], profile.p, yflip=true, label="N2O kernel")