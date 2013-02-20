
def categorize(mctype):
    if mctype.lower() in ["tw", 'ttw', 'ttbar', 'ttg', 'ttww', 'ttz']:
        return "top"
    elif mctype.lower() in ['wwto2l2nu', 'wwg', ]:
        return 'WW'
    elif mctype.lower() in ['wzto3lnu', 'wgstartolnu2e', 'wgstartolnu2mu']:
        return 'WZ'
    elif mctype.lower() in ['zzto2l2nu']:
        return 'ZZ'
    elif mctype.lower() in ['www', 'wwznogstar']:
        return "VVV"
    elif mctype.lower() in ['dy_m10', 'dyjets_m-50', 'wzto2l2q', 'zzto2l2q', 'zzznogstar', 'zzto4l']:
        return 'DY'
    elif mctype.lower() in ['wjets']:
        return 'fake'
    elif mctype.lower() in ['vbftohww', 'ggtohww', 'rarehww']:
        return "HWW"
    else :
        raise RuntimeError("Unrecognied mctype "+mctype)

# cats = s['all'].mctype.apply(categorize)
# d = s['all'].insert(len(s['all'].columns), 'mc_cat', cats)
