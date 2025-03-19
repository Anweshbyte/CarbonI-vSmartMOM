function loadTCCON(file, wl_grid)
    fts = Dataset(file);
    wn = fts["frequency"][:]
    wo = findall(1e7/2500 .< wn .< 1e7/1980)
    wn = 0.999999*wn[wo]
    sp = fts["intensity"][wo]
    wnn = range(wn[1], wn[end], length(wn));
    interr = CubicSplineInterpolation(wnn, sp)
    psurf = fts["outside_pressure"][:]
   # @show psurf, fts["sza"][:]
    return interr(1e7./wl_grid), fts["sza"][:],psurf, fts["time"][:]
end