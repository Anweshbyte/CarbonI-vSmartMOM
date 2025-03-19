
using CarbonI
using ImageFiltering, DiffResults, ForwardDiff, InstrumentOperator, Unitful, Interpolations
using NCDatasets, Polynomials, LinearAlgebra, SpecialPolynomials
# Load spectroscopies:
co2, ch4, h2o, hdo, n2o, co, co2_iso2, c2h6 = CarbonI.loadXSModels();

#include(joinpath(@__DIR__, "readSun_DC.jl"))
include(joinpath(@__DIR__, "readSun.jl"))
include(joinpath(@__DIR__, "forwardModel.jl"))

# Load some profile:
#MD = "/net/fluo/data1/ftp/XYZT_ESE156/Data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"
#MD = "/home/cfranken/for-christian/MERRA2_400.tavg3_3d_asm_Nv.20170601.nc4"
MD = CarbonI.default_merra_file
hitran_array = (co2, h2o, ch4, co, n2o, hdo, co2_iso2, c2h6);


# What latitude do we want? Take Caltech as example
myLat = 36.604
myLon = -97.486

myLat = 34.1478
myLon = -118.1445
profile_hr = CarbonI.read_atmos_profile_MERRA2(MD, myLat, myLon, 7);

# Reduce dimensions, group layers together to get roughly layers of equal pressure difference:
n_layers = 10

Δwl = 0.01
wl = 2000:Δwl:2400
σ_matrix_hr = CarbonI.compute_profile_crossSections(profile_hr, hitran_array , wl);
sol  = CubicSplineInterpolation(range(wlSol[1],wlSol[end], length=length(wlSol)),solar_irr, extrapolation_bc=Interpolations.Flat());
# Reduce to fewer dimensions:
profile, σ_matrix = CarbonI.reduce_profile(n_layers, profile_hr, σ_matrix_hr);

nL = length(profile.T)
	
vmr_co2 = zeros(nL) .+ 407e-6
vmr_ch4 = zeros(nL) .+ 1.8e-6
vmr_ch4[1:3] .= 1.4e-6
vmr_h2o = profile.vcd_h2o ./ profile.vcd_dry
vmr_co  = zeros(nL) .+ 100e-9
vmr_n2o = zeros(nL) .+ 337e-9
vmr_n2o[1:3] .= 100e-9
vmr_hdo = vmr_h2o * 0.9
vmr_c2h6 = zeros(nL) .+ 1.0e-9
vmrs = [vmr_co2, vmr_h2o, vmr_ch4,vmr_co, vmr_n2o, vmr_hdo, vmr_co2, vmr_c2h6];

# Define a polynomial scaling
p = Polynomial([0.2,0.0001,0.000001]);

# Define an instrument:
FWHM  = 1.0  # 
SSI  = 0.68
kern1 = CarbonI.box_kernel(2*SSI, Δwl)
kern2 = CarbonI.gaussian_kernel(FWHM, Δwl)
kernf = imfilter(kern1, kern2)
lociBox = CarbonI.KernelInstrument(kernf, collect(2040:SSI:2380));

# Define state vector:
x = [vmr_co2; vmr_h2o; vmr_ch4; vmr_co; vmr_n2o; vmr_hdo; vmr_co2 ; vmr_c2h6 ;zeros(10) ];

@show size(x)
sza = 45.0

result = DiffResults.JacobianResult(zeros(length(lociBox.ν_out)),x);


# Define the instrument specs:
ET  = 44.0u"ms" # Exposure time
SSI = (2*0.7)u"nm" # Spectral resolution
Pitch = 18.0u"μm" # Pixel pitch
FPA_QE = 0.9 # FPA quantum efficiency
Bench_efficiency = 0.65 # Bench efficiency
Fnumber = 2.2   # F-number
readout_noise = 70.0 # Readout noise
dark_current = 100.0u"1/s" # Dark current

ins = InstrumentOperator.createGratingNoiseModel(ET, Pitch,FPA_QE, Bench_efficiency, Fnumber, SSI, (readout_noise), dark_current);

soil = CubicSplineInterpolation(450:2500,r[:,156], extrapolation_bc=Interpolations.Flat());

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
rel_error = 0.02;
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
Sₐ[80,80] = (122x[80])^2
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

# Run model:
szas = collect(0:5:80)
alb = collect(0.02:0.02:0.5)
resis = zeros(length(szas),length(alb),8)
iSZA = 10
iAlb = 10
for iSZA in eachindex(szas), iAlb in eachindex(alb)
    @show szas[iSZA], alb[iAlb]
    sza = szas[iSZA]
    refl = alb[iAlb]/0.1*soil(wl)
    refl = alb[iAlb].+0.0*soil(wl)
    ForwardDiff.jacobian!(result, forward_model_x_, x);
    K = DiffResults.jacobian(result);
    F = DiffResults.value(result);
    # Adapt K for the legendre polynomials:
    ranger = range(-1,1,length(F))
    pp = zeros(10)
    for ii=1:10
        pp = zeros(10)
        pp[ii] = 1.0
        p = Legendre(pp)
        K[:,80+ii] .= p.(ranger)
    end
    
    
    nesr = InstrumentOperator.noise_equivalent_radiance(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
    nesr_ = nesr./1u"mW/m^2/nm/sr"
    e = InstrumentOperator.photons_at_fpa(ins, (lociBox.ν_out)u"nm", (F)u"mW/m^2/nm/sr");
    Se = Diagonal(nesr_.^2);
    # Gain matrix:
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
    resis[iSZA,iAlb,:] = [ch4_error, co2_error, h2o_error, hdo_error, n2o_error, co_error, c2h6_error, maximum(e)]

end

iAlb = argmin(abs.(alb.-0.2))
iSZA = argmin(abs.(szas.-40))

f = Figure(resolution=(600,500))

ax1 = f[1,1] = GridLayout()
ax2 = f[1, 2]= GridLayout()#,xlabel="Solar Zenith Angle (degrees)", ylabel="Surface Albedo")
ax3 = f[2, 1]= GridLayout()#,xlabel="Solar Zenith Angle (degrees)", ylabel="Surface Albedo")
ax4 = f[2, 2]= GridLayout()#xlabel="Solar Zenith Angle (degrees)", ylabel="Surface Albedo")

axtop = Axis(ax1[1,1])
axtop.title = "CH₄ error (ppb)"
axtop2 = Axis(ax2[1,1])
axtop2.title = "CO₂ error (ppm)"
axtop3 = Axis(ax3[1,1])
axtop3.title = "CO error (ppb)"
axtop4 = Axis(ax4[1,1])
axtop4.title = "N₂O error (ppb)"
c = to_color(:red)
c2 = RGBAf(c.r, c.g, c.b,0.5)
scale =sqrt(10)

ch4 = CairoMakie.contourf!(ax1[1,1],szas, alb, resis[:,:,1]/scale, levels=range(1,8,20), extendhigh=c2)
CairoMakie.scatter!(ax1[1,1],szas[iSZA], alb[iAlb], resis[iSZA,iAlb,1]/scale, color=:red)
tt = @sprintf("%.2g", resis[iSZA,iAlb,1]/scale) * "ppb"
CairoMakie.text!(ax1[1,1],szas[iSZA]+1, alb[iAlb]+0.005, text = tt, color=:white, fontsize=18, font = :bold)


co  = CairoMakie.contourf!(ax3[1,1],szas, alb, resis[:,:,6]/scale, levels=range(5,25,20), extendhigh=c2)
CairoMakie.scatter!(ax3[1,1],szas[iSZA], alb[iAlb], resis[iSZA,iAlb,1]/scale, color=:red)
tt = @sprintf("%.2g", resis[iSZA,iAlb,6]/scale) * "ppb"
CairoMakie.text!(ax3[1,1],szas[iSZA]+1, alb[iAlb]+0.005, text = tt, color=:white, fontsize=18, font = :bold)

co2 = CairoMakie.contourf!(ax2[1,1],szas, alb, resis[:,:,2]/scale, levels=range(0.4,2.0,20), extendhigh=c2)
CairoMakie.scatter!(ax2[1,1],szas[iSZA], alb[iAlb], resis[iSZA,iAlb,1]/scale, color=:red)
tt = @sprintf("%.2g", resis[iSZA,iAlb,2]/scale) * "ppm"
CairoMakie.text!(ax2[1,1],szas[iSZA]+1, alb[iAlb]+0.005, text = tt, color=:white, fontsize=18, font = :bold)
n2o = CairoMakie.contourf!(ax4[1,1],szas, alb, resis[:,:,5]/scale, levels=range(2,9,20), extendhigh=c2)
CairoMakie.scatter!(ax4[1,1],szas[iSZA], alb[iAlb],color=:red)
tt = @sprintf("%.2g", resis[iSZA,iAlb,5]/scale) * "ppb"
CairoMakie.text!(ax4[1,1],szas[iSZA]+1, alb[iAlb]+0.005, text = tt, color=:white, fontsize=18, font = :bold)
for a in (ax1, ax2, ax3, ax4)
    CairoMakie.lines!(a[1,1], [0,80],[0.05,0.05], color=:black, alpha=0.5)
end
Colorbar(ax1[1, 2], ch4)#, label="CH₄ error (ppb)")
Colorbar(ax3[1, 2], co)#, label="CO error (ppb)")
Colorbar(ax2[1, 2], co2)#, label="CO₂ error (ppm)")
Colorbar(ax4[1, 2], n2o)#, label="N₂O error (ppb)")

for a in (axtop,axtop2,axtop3,axtop4)
    CairoMakie.xlims!(a,(0,70))
    CairoMakie.ylims!(a,(0.02,0.5))
    #colgap!(a,10)
    #rowgap!(a,10)
end
for a in (ax1,ax2,ax3,ax4)
    colgap!(ax,2)
end

rowgap!(f.layout,3)
colgap!(f.layout,3)
hidexdecorations!(axtop)
hidexdecorations!(axtop2)
hideydecorations!(axtop2)
hideydecorations!(axtop4)
f
CairoMakie.save("error_budget_CarbonI.pdf", f)
