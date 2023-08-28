function rep = reach_interp(ind, start_h, end_h, reach_init, reach_fin)
    rep = round(reach_init + (reach_fin-reach_init)*(ind - start_h)/(end_h-start_h));
end