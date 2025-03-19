

# Example Gas Array (we merge CO2_13 and CO2 now) ########################
gas_array1 = ["co2", "ch4", "n2o", "h2o","hdo"];
#gas_array = ["ch4", "co2", "n2o", "h2o", "co", "hdo", "c2h6"];
#gas_array = ["ch4", "co2", "n2o",  "co",  "c2h6"];
#gas_array = ["ch4", "n2o", "h2o", "co", "hdo", "c2h6"];
setupFile = "/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml"
n_layers = 10
indLR1 = 8:260
#indLR = 287:410
#indLR = 7:410
# CLS = correlation length scale
cls        = Dict(gas => 1550.0 for gas in gas_array1)
cls["h2o"] = 150.0
cls["co2"] = 150.0
cls["hdo"] = 150.0
# Relative error of trace gases (0.15=15%, used everywhere)
rel_errors = Dict(gas => 0.15 for gas in gas_array1)
rel_errors["h2o"] = 0.15
rel_errors["hdo"] = 0.15
rel_errors["n2o"] = 0.15
# PBL error (100=amplify previous noise by this factor for the PBL (boundary layer)
pbl_error = 100.0

############################################################################



n_poly = 4
a1 = createFitParams(indLR1, setupFile, n_layers, gas_array1)
xa1, Sa1, h_column1 = createBayesianConstraints(gas_array1, a1.gasProfiles,n_poly, a1.profile, 150.0, 800.0, cls, rel_errors, pbl_error)
ff = define_forward_model(a1)
errors = 1e10*ones(length(indLR1));
errors .= 0.00002
Se1 = Diagonal(errors.^2);
#Sa[15,15] = 100^2
iP = length(xa1)-n_poly+1
ii_mw1 = define_inverse_model(ff,xa1,Sa1,Se1,length(indLR1),4,iP);

# Example Gas Array (we merge CO2_13 and CO2 now) ########################
#gas_array = ["co2", "ch4", "n2o", "h2o","hdo"];
gas_array = ["ch4", "co2", "n2o", "h2o", "co", "hdo", "c2h6"];
#gas_array = ["ch4", "co2", "n2o",  "co",  "c2h6"];
#gas_array = ["ch4", "n2o", "h2o", "co", "hdo", "c2h6"];
setupFile = "/home/cfranken/code/gitHub/CarbonI/src/yaml/carbon-i.yaml"
n_layers = 10
#indLR = 8:260
indLR2 = 260:493
#indLR = 7:410
cls        = Dict(gas => 1550.0 for gas in gas_array)
cls["h2o"] = 150.0
cls["co2"] = 150.0
cls["hdo"] = 150.0
rel_errors = Dict(gas => 0.15 for gas in gas_array)
rel_errors["h2o"] = 0.15
rel_errors["hdo"] = 0.15
rel_errors["n2o"] = 0.15
pbl_error = 100.0

############################################################################



n_poly = 6
a = createFitParams(indLR2, setupFile, n_layers, gas_array)
xa, Sa, h_column = createBayesianConstraints(gas_array, a.gasProfiles,n_poly, a.profile, 150.0, 800.0, cls, rel_errors, pbl_error)
ff = define_forward_model(a)
errors = 1e10*ones(length(indLR2));
errors .= 0.00002
Se = Diagonal(errors.^2);
#Sa[15,15] = 100^2
iP = length(xa)-n_poly+1
ii_mw2 = define_inverse_model(ff,xa,Sa,Se,length(indLR2),4,iP);



