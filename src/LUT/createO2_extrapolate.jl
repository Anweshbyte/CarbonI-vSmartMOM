o2_i = load_interpolation_model("/net/fluo/data3/data/Databases/ABSCO_CS_Database/v5.2_final/o2_v52.jld2")
o2_e = extrapolate(o2_i.itp, 4e-28)
o2_extra = InterpolationModel(o2_e, o2_i.mol, o2_i.iso, o2_i.ν_grid, o2_i.p_grid, o2_i.t_grid);
save_interpolation_model(o2_extra, "/net/fluo/data3/data/Databases/ABSCO_CS_Database/v5.2_final/o2_v52_extra.jld2")