using Glob, Dates
spec,sza,psurf, timer = loadTCCON("/home/cfranken/for-christian/oc20170601sfddaa.1001.nc",wl)
solarIrr = sol(wl);
p = Legendre([2.0,0.0001,0.000001,0.0001,0.000001,0.0001,0.000001,0.0001,0.000001,0.0001,0.000001]);
x2 = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; vmr_hdo;vmr_co2; 1.0;  p[:] ];
result = DiffResults.JacobianResult(zeros(length(lociBox.ν_out)),x2);
F = DiffResults.value(result);
sza = 73.01;

# Get prior covariance matrix:
n_state = length(x2);
Sₐ = zeros(n_state,n_state);
rel_error = 0.1;
# vcd_ratio = profile_caltech.vcd_dry ./ mean(profile_caltech.vcd_dry)
	
# Fill the diagonal for the trace gases:
for i=1:70
	Sₐ[i,i] = (rel_error*x2[i])^2   
end
# For H2O
for i=11:19
	Sₐ[i,i] = (0.5*x2[i])^2   
end
for i=41:49
	Sₐ[i,i] = (0.5*x2[i])^2   
end

# CO2 at surface, 100% error
Sₐ[10,10] = (20*x2[10])^2
Sₐ[20,20] = (20*x2[20])^2
Sₐ[30,30] = (20*x2[30])^2
Sₐ[40,40] = (20*x2[40])^2
Sₐ[50,50] = (20*x2[50])^2
Sₐ[60,60] = (20*x2[60])^2
Sₐ[70,70] = (20*x2[70])^2
# Put in arbitrarily high numbers for the polynomial term, so these won't be constrained at all! 
for i=71:n_state
	Sₐ[i,i] = 1e5;
end
xa = x2;

ForwardDiff.jacobian!(result, forward_model_x, x2);
#Se = Diagonal(sqrt.(^2);
K = DiffResults.jacobian(result);
F = DiffResults.value(result);
k = 1000000.0
Se = Diagonal((sqrt.(F*k)/k).^2);

y = CarbonI.conv_spectra(lociBox, wl, spec);

G = inv(K'inv(Se)K + inv(Sₐ))K'inv(Se);
y_prime = y - F;
x_primeBayes = G * y_prime;
x̂_bayes = x_primeBayes + x2;
Ŝ = inv(K'inv(Se)K + inv(Sₐ));

ratio = profile.vcd_dry/sum(profile.vcd_dry);

h_co2 = zeros(length(x2));
h_co2_ = zeros(length(x2));
h_ch4 = zeros(length(x2));
h_h2o = zeros(length(x2));
h_co  = zeros(length(x2));
h_hdo = zeros(length(x2));
h_n2o = zeros(length(x2));

h_co2[1:10] .= ratio;
h_h2o[11:20] .= ratio;
h_ch4[21:30] .= ratio;
h_co[31:40] .= ratio;
h_n2o[41:50] .= ratio;
h_hdo[51:60] .= ratio;
h_co2_[61:70] .= ratio;

plot(lociBox.ν_out, K[:,50]/1e6,linewidth=2, alpha=0.5, label="N₂O")
plot!(lociBox.ν_out, K[:,40]/1e6,linewidth=2, alpha=0.5, label="CO")
plot!(lociBox.ν_out, K[:,30]/1e6,linewidth=2, alpha=0.5, label="CH₄")
plot!(lociBox.ν_out, K[:,10]/1e4,linewidth=2, alpha=0.5, label="CO₂")
plot!(lociBox.ν_out, K[:,70]/1e4,linewidth=2, alpha=0.5, label="¹³CO₂")
plot!(lociBox.ν_out, K[:,60]/1e3,linewidth=2, alpha=0.5, label="HDO")
plot!(lociBox.ν_out, K[:,20]/1e3,linewidth=2, alpha=0.5, label="H₂O")

for h in (h_co2, h_h2o, h_ch4, h_co, h_n2o, h_hdo, h_co2_)
    @show sqrt(h' * Ŝ * h)*1e6, h'*x_all[:,end]*1e6,h'*xa*1e6
end 





# Define an instrument:
FWHM  = 7.0  # 
SSI  = 7.1
kern1 = CarbonI.box_kernel(SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
lociBox = CarbonI.KernelInstrument(kernf, collect(2035.0:SSI:2380));

N = length(lociBox.ν_out)
result = DiffResults.JacobianResult(zeros(N),x2);
max_no_of_iter = 6
x_all   = zeros((length(x2),max_no_of_iter+1))
F_all   = zeros((N,max_no_of_iter))
x_all[:,1]=x2

files = glob("ci20*.nc", "/home/cfranken/for-christian-20200610/")
iFile = 10
spec,sza,psurf, time = loadTCCON(files[iFile],wl)
y = CarbonI.conv_spectra(lociBox, wl, spec);

Se = Diagonal((sqrt.(y*k)/k).^2)
# ╔═╡ 3865ad63-2cb2-442b-86c1-70916a9dcab4
for i=1:max_no_of_iter
    @show i
    #print('Iteration #',i)
	ForwardDiff.jacobian!(result, forward_model_x, x_all[:,i]);
	Kᵢ = DiffResults.jacobian(result);
    Fᵢ = DiffResults.value(result);
    G = inv(Kᵢ'inv(Se)Kᵢ + inv(Sₐ))Kᵢ'inv(Se);
    x_all[:,i+1] = xa + G * (y - Fᵢ + Kᵢ *(x_all[:,i]-xa))
    #@show h_co2'*x_all[:,end,end]*1e6
    F_all[:,i] = Fᵢ
end
f = Figure(resolution=(600,400))
ax = Axis(f[1,1], title="Carbon-I Fit to TCCON Spectrum", ylabel="Radiance [AU]")
ax_ = Axis(f[2,1],xlabel="Wavenlength [nm]",ylabel="Modeled-Measured" )
ax2 = lines!(ax,lociBox.ν_out, F_all[:,end],label="Carbon-I fit")
lines!(ax,lociBox.ν_out, y, label="Carbon-I spectrum (convolved TCCON data)")
lines!(ax_,lociBox.ν_out, (F_all[:,end]-y)/0.25*100, label="Residual (%)")
ylims!(ax,0.08,0.28)
ylims!(ax_,-2,2)
xlims!(ax_,2040,2380)
xlims!(ax_,2040,2380)
axislegend(ax, position=:cb)
#axislegend(ax_,position=:rt)
hidexdecorations!(ax, grid=false)
rowsize!(f.layout,1,Relative(0.7))
rowgap!(f.layout,5)
f

#files = glob("oc20*.nc", "/home/cfranken/for-christian/20200610/");
files = glob("ci20*.????.nc", "/home/cfranken/for-christian-20200610/")
time_all = a = zeros(Dates.DateTime, length(files))
x_all   = zeros((length(x2),max_no_of_iter+1,length(files)))
F_all   = zeros((N,max_no_of_iter,length(files)))
y_all   = zeros((N,length(files)))
sza_all = zeros(length(files))
p_all = zeros(length(files))
for iFile in eachindex(files)
    spec,sza,psurf, time = loadTCCON(files[iFile],wl)
    N = length(lociBox.ν_out)
    max_no_of_iter = 6
    y = CarbonI.conv_spectra(lociBox, wl, spec);
    y_all[:,iFile] = y;
    x_all[:,1,iFile]=x2
    sza_all[iFile] = sza
    p_all[iFile] = psurf
    for i=1:max_no_of_iter
        @show i, @show h_co2'*x_all[:,i,iFile]*1e6 @show h_ch4'*x_all[:,i,iFile]*1e6
        #print('Iteration #',i)
        ForwardDiff.jacobian!(result, forward_model_x, x_all[:,i,iFile]);
        Kᵢ = DiffResults.jacobian(result);
        Fᵢ = DiffResults.value(result);
        Gain = inv(Kᵢ'inv(Se)Kᵢ + inv(Sₐ))Kᵢ'inv(Se);
        x_all[:,i+1,iFile] = xa + Gain * (y - Fᵢ + Kᵢ *(x_all[:,i,iFile]-xa))
        
        F_all[:,i,iFile] = Fᵢ
    end
end
