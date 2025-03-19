pa = Dataset("/home/cfranken/pa20040602_20221205.public.qc.nc")
timePA = pa["time"][:]
woP = findall(timePA .<= time_all[end])