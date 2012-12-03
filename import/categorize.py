
def categorize(mctype):
    if mctype.lower() in ["tw", 'ttw', 'ttbar', 'ttg', 'ttww', 'ttz']:
        return "top"
    elif mctype.lower() in ['wwto2l2nu', 'wwg', 'www', 'wzto3lnu']:
        return 'WV'
    elif mctype.lower() in ['zzto2l2nu']:
        return 'ZZ'
    elif mctype.lower() in ['dy_m10', 'dyjets_m-50', 'wz', 'zzto2l2q', 'wwznogstar', 'zzznogstar', 'zzto4l']:
        return 'DY'
    elif mctype.lower() in ['wjets']:
        return 'fake'
    else :
        raise RuntimeError("Unrecognied mctype "+mctype)

# cats = s['all'].mctype.apply(categorize)
# d = s['all'].insert(len(s['all'].columns), 'mc_cat', cats)
