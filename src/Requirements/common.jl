using CarbonI
using Unitful
using DiffResults
using DelimitedFiles
using ForwardDiff
using InstrumentOperator
using Interpolations

include(joinpath(dirname(pathof(CarbonI)), "forwardModel.jl"))

function calc_rel_error(specs, x, solarIrr, refl, sza, σ_matrix, profile, h, ins, Sₐ; return_F=false)

	# Run Forward model and Jacobian generation together:
	result = DiffResults.JacobianResult(zeros(length(specs.instrument_wl)),x);

	ForwardDiff.jacobian!(result, (x) -> 
	forward_model_x_(x, 
				sun=solarIrr, 
				instrument=specs.instrument_kernel,
				reflectance=refl, 
				sza=sza, 
				σ_matrix=σ_matrix, 
				profile=profile,
				wl=specs.modelling_wl), x);
	K = DiffResults.jacobian(result);
	F = DiffResults.value(result);  

	# Compute noise equivalent radiance:
	nesr = InstrumentOperator.noise_equivalent_radiance(ins,(specs.instrument_wl)u"nm", (F)u"mW/m^2/nm/sr");
	nesr_ = nesr./1u"mW/m^2/nm/sr"

	# Generate S\_epsilon matrix
	Se = Diagonal(nesr_.^2);

	# Compute the Gain Matrix:
	#G = inv(K'inv(Se)K + inv(Sₐ))K'inv(Se);

	# Posterior covariance matrix (Profile retrievals, 10 layers):
	Ŝ = inv(K'inv(Se)K + inv(Sₐ));

	# Apply column operator to get to total column uncertainties:
	error_ppb = Dict()
	for spec in ["co2","co213", "ch4", "h2o", "co", "hdo", "n2o", "c2h6"]
		error_ppb[spec] = sqrt(h[spec]' * Ŝ * h[spec])*1e9
	end
    if return_F
        return error_ppb, F
    else
        return error_ppb
    end
	
	
end


function setup_data(scenario, specs)

	# Load spectroscopies:
	co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();


	DS = Dataset(CarbonI.solar_file);
	wlSol = 1e3*DS["wl"][:]
	solar_irr = 1e3*DS["solar_irr"][:] # convert to mW/m2/nm
	close(DS)

	# optionally load a custom reflectance
	#DS = Dataset("data/reflectance_cube_all_1nm_from450nm.h5")
	#r = DS["reflectance_cube"][:]
	#close(DS)
	
	hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);
	# Precompute the cross sections:
	σ_matrix_hr = CarbonI.compute_profile_crossSections(scenario.profile_hr, hitran_array , specs.modelling_wl);

	nL = length(scenario.profile_hr.T)

	vmr_co2 = zeros(nL) .+ 407e-6
	vmr_ch4 = zeros(nL) .+ 1.8e-6
	vmr_ch4[1:3] .= 1.4e-6
	vmr_h2o = scenario.profile_hr.vcd_h2o ./ scenario.profile_hr.vcd_dry
	vmr_co  = zeros(nL) .+ 100e-9
	vmr_n2o = zeros(nL) .+ 337e-9
	vmr_n2o[1:3] .= 100e-9
	vmr_hdo = vmr_h2o * 0.9
	vmr_c2h6 = zeros(nL) .+ 1.0e-9
	vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, vmr_hdo, vmr_co2, vmr_c2h6];

	sol  = CubicSplineInterpolation(range(wlSol[1],wlSol[end], length=length(wlSol)),solar_irr, extrapolation_bc=Interpolations.Flat());

	# Reduce dimensions, group layers together to get  layers of equal pressure difference:
	n_layers = 10;

	# Reduce profile and cross sections to fewer dimensions:
	profile, σ_matrix, indis, gasProfiles = CarbonI.reduce_profile(n_layers, scenario.profile_hr, σ_matrix_hr,vmrs);

	# Number of Legendre Polynomials for the spectrally resolved albedo term
	nLeg = 10

	# Define state vector elements for legendre polynomials
	xPoly = zeros(nLeg).+eps()

	# First element has to be unity here
	xPoly[1] = 1.0

	# Define full state vector (includes trace gas VMRs and Legendre Polynomial)
	x = [reduce(vcat,gasProfiles) ; xPoly ];

	# Define Covariance matrices and column kernel operator
	# Get prior covariance matrix:
	n_state = length(x);
	Sₐ = zeros(n_state,n_state);
	rel_error = 0.0001;
	rel_error = 0.02;
	# could try modifying this to something like 0.02 (TBD)
	# currently good for point surfaces and surface components

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

	h = Dict()
	for spec in ["co2","co213", "ch4", "h2o", "co", "hdo", "n2o", "c2h6"]
		h[spec] = zeros(length(x));
	end

	h["co2"][1:10] .= ratio;
	h["h2o"][11:20] .= ratio;
	h["ch4"][21:30] .= ratio;
	h["co"][31:40] .= ratio;
	h["n2o"][41:50] .= ratio;
	h["hdo"][51:60] .= ratio;
	h["co213"][61:70] .= ratio;
	h["c2h6"][71:80] .= ratio;

	# Load tropical albedo
	clima_alb = readdlm(CarbonI.albedo_file,',', skipstart=1);

	# Generate albedo function generator
	soil = CubicSplineInterpolation(300:2400,clima_alb[:,2]/1.16, extrapolation_bc=Interpolations.Flat());

	# Create solar spectrum at Forward model resolution
	solarIrr = sol(specs.modelling_wl);


	return soil, x, solarIrr, σ_matrix, profile, h, Sₐ

end

