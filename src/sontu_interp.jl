using LegendrePolynomials
using DelimitedFiles
using Interpolations, InstrumentOperator

fit_window = [758.0, 759.2]
#sza_str   = "50"
#psurf_str = "1000"
#alb_str   = "0p2"
FT = Float32
a=[]
for i=0:20
    push!(a, acosd(i/20))
end
sza=reverse(Int.(ceil.(a[8:21])))
ρ = zeros(FT,15)
ρ_str = []
for iρ = 1:15
    ρ[iρ] = (iρ-1)*0.05
    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
end
psurf=[1000 750 500]
xSIF = zeros(length(psurf), length(ρ_str), length(sza))
xSIF_noRS = zeros(length(psurf), length(ρ_str), length(sza))
addRaman = zeros(length(psurf), length(ρ_str), length(sza))
contRad = similar(addRaman)
isurf = 1
iρ = 9
iA = 5
for isurf = 1:3 # 1:1 #2:2 #
    for iρ = 1:15 #3 #1:15
        for iA = 1:14
            sza_str  = string(sza[iA])
            alb_str  = ρ_str[iρ]
            psurf_str= string(psurf[isurf])
            fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"
            fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"
            specNoRS = readdlm(fname0)
            specRRS  = readdlm(fname1)
            #specNoRS0 = readdlm(fname0_0)
            #specRRS0  = readdlm(fname1_0)
            wl = 1e7./specNoRS[:,1]
            ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);
            radRRS_I = specRRS[ind,2] + specRRS[ind,5]
            radnoRS_I = specNoRS[ind,2]
            #radRRS_I0 = specRRS0[ind,2] + specRRS0[ind,5]
            #radnoRS_I0 = specNoRS0[ind,2]
            F0 = specNoRS[ind,5]
            # Find indices for wavelength window:
            # Fit SIF
            # Grid for Legendre Polynomials:
            iLeg = range(-1,1,length(ind))
            # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn’t matter)
            #Io = specNoRS[ind,5] #radnoRS_I[ind]
            # Define polynomial terms
            poly = Pl.(iLeg,(0:3)');
            # Multiply polynomial with solar spectrum
            K_ = F0 .* poly
            # Add a constant SIF term to the Jacobian (can add Shape later)
            K = [K_ ones(length(ind))];
            # Fit with Least squares:
            x = K \ radRRS_I;
            x_noRS = K \ radnoRS_I;
            #x0 = K \ radRRS_I0[ind];
            #x0_noRS = K \ radnoRS_I0[ind];
            xSIF[isurf, iρ, iA] = x[end]
            xSIF_noRS[isurf, iρ, iA] = x_noRS[end]
            addRaman[isurf, iρ, iA] =mean(specRRS[ind,5])
            contRad[isurf, iρ, iA] =mean(specRRS[ind,2])
            #xSIF0[isurf, iρ, iA] = x0[end]
            #xSIF0_noRS[isurf, iρ, iA] = x0_noRS[end]
        end
    end
end

itp = interpolate(xSIF[end:-1:1,:,end:-1:1], BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp, 500.0:250:1000, 0:0.05:0.7, 0.35:0.05:1.0)

x = 800.0
y = 0:0.01:0.7
z = 0.35:0.01:1.0
contourf(y,z,preFac*sitp.(x,y,z'),levels=10)


L1File   = "/net/fluo/data3/data/FluoData1/group/oco2/L1bSc/oco2_L1bScND_26780a_190715_B10003r_200429212407.h5"
metFile  = "/net/fluo/data3/data/FluoData1/group/oco2/L2Met/oco2_L2MetND_26780a_190715_B10003r_200429212406.h5"
dictFile = "/home/cfranken/code/gitHub/InstrumentOperator.jl/json/oco2.yaml"

oco = InstrumentOperator.load_L1(dictFile,L1File, metFile);

# Pick some bands as tuple (or just one)
bands = (1,);
#bands = (1,3);
# Indices within that band:
indices = (10:925,);
#indices = (92:885,50:916);
# Geo Index (footprint,sounding):
GeoInd = [5,5000];

# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd);

wl = 1e7./specNoRS[:,1]
i = sortperm(wl)
λ_grid = wl[i]  ./1e3
spec = specNoRS[end:-1:1,2] 

interp_I = LinearInterpolation(λ_grid, reverse(spec[i]));
res = 0.001e-3;
off = 0.5e-3
wl = oco_sounding.ils[1].ν_out[1]-off:res:oco_sounding.ils[1].ν_out[end]+off;
@show wl[1],wl[end], λ_grid[1],λ_grid[end]
I_wl = interp_I(wl);

# Convolve input spectrum with variable kernel
@time I_conv = InstrumentOperator.conv_spectra(oco_sounding.ils[1], wl, I_wl)


wl = 1e7./specRRS[:,1]
i = sortperm(wl)
λ_grid = wl[i]  ./1e3
spec2 = specRRS[end:-1:1,2] +  specRRS[end:-1:1,5]
wl = oco_sounding.ils[1].ν_out[1]-off:res:oco_sounding.ils[1].ν_out[end]+off;
interp_I2 = LinearInterpolation(λ_grid, reverse(spec2[i]));
I_wl2 = interp_I2(wl);

# Convolve input spectrum with variable kernel
@time I_convRRS = InstrumentOperator.conv_spectra(oco_sounding.ils[1], wl, I_wl2)


# Fit:

ind = findall(x->x>fit_window[1] && x<fit_window[2], oco_sounding.SpectralGrid*1e3);
iLeg = range(-1,1,length(ind))
poly = Pl.(iLeg,(0:3)');
# Multiply polynomial with solar spectrum
K_ = I_conv[ind] .* poly
# Add a constant SIF term to the Jacobian (can add Shape later)
K = [K_ ones(length(ind))];
# Fit with Least squares:
x = K \ I_convRRS[ind];


Makie.set_theme!(ggthemr(:fresh))
f = Figure(resolution=(800,430))
ax2 = Axis(f[1,1], halign=:left, title="SIF vs. Raman impact", ylabel="Radiance", xlabel="Wavelength (nm)")
ax1 = Axis(f[1,2], halign=:right, title="SIF vs. Raman impact", ylabel="Radiance", xlabel="Wavelength (nm)", yaxisposition=:right)

preFac = (specRRS[:,1].^2/1e7)
SIF = π*(specRRS[:,1].^2/1e7).*(specRRS[:,2]+specRRS[:,5]-(specRRS0[:,2]+specRRS0[:,5]))

lines!(ax1,1e7./specRRS[:,1],SIF,color=:gray, linewidth=2, label="SIF")   
lines!(ax1,1e7./specRRS[:,1],preFac .* specRRS[:,5],color=:red, alpha=0.5,linewidth=1, label="Raman")

preFacN = (specRRS0_B[:,1].^2/1e7)
lines!(ax2,wl_B,preFacN .* specRRS0_B[:,5],color=:red, alpha=0.5,linewidth=1, label="Raman")
#axislegend(ax1, position=:cb, framevisible=false)
axislegend(ax1, position=:lb, framevisible=false)
#ax1.yticklabelsvisible = false

CairoMakie.xlims!(ax1,746.0,780)
CairoMakie.xlims!(ax2,670.0,700)
CairoMakie.ylims!(ax1,0.0,0.6)
CairoMakie.ylims!(ax2,0.0,0.6)

fname0 = "/home/sanghavi/RamanSIFgrid/raylSIF_sza37_alb0p4_psurf1000hpa_nors_ABO2.dat"
fname1 = "/home/sanghavi/RamanSIFgrid/raylSIF_sza37_alb0p4_psurf1000hpa_rrs_ABO2.dat"
fname0_0 = "/home/sanghavi/RamanSIFgrid/rayl_sza37_alb0p4_psurf1000hpa_nors_ABO2.dat"
fname0_1 = "/home/sanghavi/RamanSIFgrid/rayl_sza37_alb0p4_psurf1000hpa_rrs_ABO2.dat"

specNoRS = readdlm(fname0)
specRRS  = readdlm(fname1)
specNoRS0 = readdlm(fname0_0)
specRRS0  = readdlm(fname1_0)
wl = 1e7./specNoRS[:,1]

fname0_B   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza37_alb0p4_psurf1000hpa_nors_BBO2.dat"
fname1_B   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza37_alb0p4_psurf1000hpa_rrs_BBO2.dat"
fname0_0_B = "/home/sanghavi/RamanSIFgrid/rayl_sza37_alb0p4_psurf1000hpa_nors_BBO2.dat"
fname0_1_B = "/home/sanghavi/RamanSIFgrid/rayl_sza37_alb0p4_psurf1000hpa_rrs_BBO2.dat"

specNoRS_B = readdlm(fname0_B)
specRRS_B  = readdlm(fname1_B)
specNoRS0_B = readdlm(fname0_0_B)
specRRS0_B  = readdlm(fname0_1_B)
wl_B = 1e7./specRRS0_B[:,1]