g = mc.groupby('mc_cat')

for mcname, group in g:
    sel = get_samples(group)

    print mcname, sum(group[sel['sig_sf']].weight), sum(group[sel['sig_mct_low_sf']].weight), sum(group[sel['sig_mct_high_sf']].weight)
