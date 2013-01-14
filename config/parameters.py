lumi = 19490. # in pbinv, MUST BE A FLOAT

ee_trigger_eff = 0.95
mumu_high_eta_trigger_eff = 0.83
mumu_low_eta_trigger_eff = 0.91
emu_trigger_eff = 0.92

ee_trigger_eff_frac_unc = 0.05
mumu_trigger_eff_frac_unc = 0.05
emu_trigger_eff_frac_unc = 0.05

ele_selection_eff_frac_unc = 0.02
mu_selection_eff_frac_unc = 0.02

bkg_colors = {"top": (.9,.6,0),
          "WW": (.35, .7, .9),
          "WZ": (0,.6,.5),
          "ZZ": (.95, .9, .25),
          "Z": (0, .45, .70),
          "DY": (0, .45, .7),
          "wjets": (.8,.4, .0)
         }

bkg_labels = {"top": "Top",
              "WW": "WW",
              "WZ": "WZ",
              "ZZ": "ZZ",
              "Z": r"Z/$\gamma^*$",
              "DY": r"Z/$\gamma^*$",
              "wjets": "Non-prompt",
              "fake": "Non-prompt"
             }
