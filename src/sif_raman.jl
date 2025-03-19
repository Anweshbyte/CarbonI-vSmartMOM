using Plots, LegendrePolynomials, DelimitedFiles, Polynomials, Statistics

fit_window = [758.0, 759.2]
sza_str = "50"
psurf   = "1000"
albedo  = "0p2"

fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*psurf*"hpa_nors_ABO2.dat"
fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*psurf*"hpa_rrs_ABO2.dat"
#fnameSIF = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*psurf*"hpa_rrs_ABO2.dat"
fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*albedo*"_psurf"*psurf*"hpa_nors_ABO2.dat"
fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*albedo*"_psurf"*psurf*"hpa_rrs_ABO2.dat"
#fnameSIF = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*albedo*"_psurf"*psurf*"hpa_rrs_ABO2.dat"

specNoRS = readdlm(fname0)
specRRS  = readdlm(fname1)

specNoRS0 = readdlm(fname0_0)
specRRS0  = readdlm(fname1_0)

wl = 1e7./specNoRS[:,1]
radRRS_I = specRRS[:,2] + specRRS[:,5] 
radnoRS_I = specNoRS[:,2] 

radRRS_I0 = specRRS0[:,2] + specRRS0[:,5] 
radnoRS_I0 = specNoRS0[:,2] 

F0 = specNoRS[:,5] 
# Find indices for wavelength window:
ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);

# Fit SIF
# Grid for Legendre Polynomials:
iLeg = range(-1,1,length(ind))
# Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
Io = specNoRS0[ind,2] #radnoRS_I[ind]
# Define polynomial terms
poly = Pl.(iLeg,(0:3)');
# Multiply polynomial with solar spectrum
K_ = Io .* poly 
# Add a constant SIF term to the Jacobian (can add Shape later)
# Crude, danach direkt nehmen vom SIF spectrum:
p_sif = fit(specNoRS[ind,1], specNoRS[ind,2]-specNoRS0[ind,2], 1)
shape = p_sif.(specNoRS[ind,1])
shape = shape ./ mean(shape)
K = [K_ shape];

# Fit with Least squares:
x = K \ radRRS_I[ind];

# COmpute fitted spectrum:
radRRS_I_fit = K * x;

l = @layout [a ; b]
p1 = plot(wl[ind], radRRS_I[ind], label="Measured")
plot!(wl[ind], radRRS_I_fit, label="Fitted")
p2 = plot(wl[ind], (radRRS_I[ind] .- radRRS_I_fit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)

# Fitted SIF:
SIF_fit = x[end]
relative_SIF = SIF_fit / maximum(radRRS_I[ind])


# Fit without SIF:
#x2 = K_ \ radRRS_I[ind];
#radStupidFit = K_ * x2;
#l = @layout [a ; b]
# Fit with Least squares:
x0 = K \ radRRS_I0[ind];

# COmpute fitted spectrum:
radRRS_I_fit0 = K * x0;


p1 = plot(wl[ind], radRRS_I0[ind], label="Measured")
plot!(wl[ind], radRRS_I_fit0, label="Fit without SIF")
p2 = plot(wl[ind], (radRRS_I0[ind] .- radRRS_I_fit0)./radRRS_I0[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)

#plot!(wl[ind], radStupidFit, label="Fit without SIF")
#p2 = plot(wl[ind], (radRRS_I[ind] .- radStupidFit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
#plot(p1, p2,  layout = l)
 

# Fit without Raman:
#x2 = K_ \ radRRS_I[ind];
#radStupidFit = K_ * x2;
#l = @layout [a ; b]
# Fit with Least squares:
x_noRS = K \ radnoRS_I[ind];

# COmpute fitted spectrum:
radnoRS_I_fit = K * x_noRS;


p1 = plot(wl[ind], radnoRS_I[ind], label="Measured")
plot!(wl[ind], radnoRS_I_fit, label="Fit without SIF")
p2 = plot(wl[ind], (radnoRS_I[ind] .- radnoRS_I_fit)./radnoRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)

#plot!(wl[ind], radStupidFit, label="Fit without SIF")
#p2 = plot(wl[ind], (radRRS_I[ind] .- radStupidFit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
#plot(p1, p2,  layout = l)
 
x_noRS0 = K \ radnoRS_I0[ind];

# COmpute fitted spectrum:
radnoRS_I_fit0 = K * x_noRS0;


p1 = plot(wl[ind], radnoRS_I0[ind], label="Measured")
plot!(wl[ind], radnoRS_I_fit0, label="Fit without SIF")
p2 = plot(wl[ind], (radnoRS_I0[ind] .- radnoRS_I_fit0)./radnoRS_I0[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)