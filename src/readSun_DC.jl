#DS = Dataset("data/solar_irr.nc")
using DelimitedFiles
sol = readdlm("data/solar_center_carbonI.dat")
wlSol = sol[:,1]
solar_irr = sol[:,2] # convert to mW/m2/nm
#close(DS)

DS = Dataset("data/reflectance_cube_all_1nm_from450nm.h5")
r = DS["reflectance_cube"][:]
close(DS)