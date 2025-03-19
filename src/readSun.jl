DS = Dataset("data/solar_irr.nc")
wlSol = 1e3*DS["wl"][:]
solar_irr = 1e3*DS["solar_irr"][:] # convert to mW/m2/nm
close(DS)

#DS = Dataset("data/reflectance_cube_all_1nm_from450nm.h5")
#r = DS["reflectance_cube"][:]
#close(DS)