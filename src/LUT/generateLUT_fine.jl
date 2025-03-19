using vSmartMOM.Absorption
using vSmartMOM
using DelimitedFiles
using Interpolations

# Minimum wavenumber
ν_min  = 4000.0
# Maximum wavenumber
ν_xmax = 5050.0

pressures = 0.0001:25:1150.0
temperatures = 160:10:360.0
grid = 2030:0.01:2390
ν_grid = ν_min:0.005:ν_xmax
lineModel = Voigt()

# http://bernath.uwaterloo.ca/publicationfiles/2019/C2H6-overtones-2020.pdf

co2_par = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=ν_min, ν_max=ν_xmax);

ch4_par = read_hitran(artifact("CH4"), mol=6, iso=1, ν_min=ν_min, ν_max=ν_xmax);
ch4_13par = read_hitran(artifact("CH4"), mol=6, iso=2, ν_min=ν_min, ν_max=ν_xmax);
h2o_par = read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=ν_min, ν_max=ν_xmax);
hdo_par = read_hitran(artifact("H2O"), mol=1, iso=4, ν_min=ν_min, ν_max=ν_xmax);
co_par  = read_hitran(artifact("CO"),  ν_min=ν_min, ν_max=ν_xmax);
n2o_par = read_hitran(artifact("N2O"),  ν_min=ν_min, ν_max=ν_xmax);

ethane = readdlm("data/C2H6_overtone/N2/295/ethanepressure697torr_n2pressure101torr_temperature295k.txt",',')
egrid = range(ethane[1,1], ethane[end,1], length(ethane[:,1]))
interp_linear = LinearInterpolation(egrid, ethane[:,2])
ethane_cs = interp_linear(1e7./wl) * 1e-21

# 1.6 band for CH4
w_ch4_par = read_hitran(artifact("CH4"), mol=6, iso=1, ν_min=1e7/1800, ν_max=1e7/1600);
w_ch4_13par = read_hitran(artifact("CH4"), mol=6, iso=2, ν_min=1e7/1800, ν_max=1e7/1600);

model_interp_co2     = make_interpolation_model(co2_par,   lineModel,ν_grid, pressures, temperatures)
model_interp_ch4     = make_interpolation_model(ch4_par,   lineModel,ν_grid, pressures, temperatures)
model_interp_ch4_13  = make_interpolation_model(ch4_13par, lineModel,ν_grid, pressures, temperatures)
model_interp_h2o     = make_interpolation_model(h2o_par,   lineModel,ν_grid, pressures, temperatures)
model_interp_hdo     = make_interpolation_model(hdo_par,   lineModel,ν_grid, pressures, temperatures)
model_interp_co      = make_interpolation_model(co_par,    lineModel,ν_grid, pressures, temperatures)
model_interp_n2o     = make_interpolation_model(n2o_par,   lineModel,ν_grid, pressures, temperatures)

save_interpolation_model(model_interp_co2,    "data/co2_model.jld2")
save_interpolation_model(model_interp_ch4,    "data/ch4_model.jld2")
save_interpolation_model(model_interp_ch4_13, "data/ch413_model.jld2")
save_interpolation_model(model_interp_h2o,    "data/h2o_model.jld2")
save_interpolation_model(model_interp_hdo,    "data/hdo_model.jld2")
save_interpolation_model(model_interp_co,     "data/co_model.jld2")
save_interpolation_model(model_interp_n2o,    "data/n2o_model.jld2")



wl = 1990:0.01:2480
offset = 1e-32
plot(wl, offset.+absorption_cross_section(model_interp_co2, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="CO₂")
plot!(wl, offset.+absorption_cross_section(model_interp_ch4, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="CH₄")
plot!(wl, offset.+absorption_cross_section(model_interp_ch4_13, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="¹³CH₄")
plot!(wl, offset.+absorption_cross_section(model_interp_h2o, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="H₂O")
plot!(wl, offset.+absorption_cross_section(model_interp_hdo, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="HDO")
plot!(wl, offset.+absorption_cross_section(model_interp_co, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="CO")
plot!(wl, offset.+absorption_cross_section(model_interp_n2o, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="N₂O")
plot!(wl, ethane_cs,  yscale=:log10, label="C2H6")

wl = 1990:0.01:2480
offset = 1e-32
plot(wl, offset.+absorption_cross_section(co2_par, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="CO₂")
plot!(wl, offset.+absorption_cross_section(ch4_par, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="CH₄")
plot!(wl, offset.+absorption_cross_section(ch4_13par, reverse(wl), 1010.0, 298; wavelength_flag=true), yscale=:log10, label="¹³CH₄")