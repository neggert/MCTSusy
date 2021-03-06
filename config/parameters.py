lumi = 19490. # in pbinv, MUST BE A FLOAT

ee_trigger_eff = 0.95
mumu_high_eta_trigger_eff = 0.81
mumu_low_eta_trigger_eff = 0.90
emu_trigger_eff = 0.93

ee_trigger_eff_frac_unc = 0.05
mumu_trigger_eff_frac_unc = 0.05
emu_trigger_eff_frac_unc = 0.05

ele_selection_eff_frac_unc = 0.02
mu_selection_eff_frac_unc = 0.02

bkg_colors = {"top": (.9,.6,0),
          "WW": (.35, .7, .9),
          "diboson": (.35, .7, .9),
          "WZ": (0,.6,.5),
          "FS": (.35, .7, .9),
          "ZZ": (.95, .9, .25),
          'Rare': (0, .45, .70),
          'Z': (.80,.40,0),
          'DY': (.80,.40,0),
          "wjets": (.80, .60, .7),
          "fake": (.80, .60, .7),
          # "fake": (0.65, 0.48, 0.36)
         }

bkg_labels = {"top": "Top",
              "WW": "WW",
              "WZ": "WZ",
              "ZZ": "ZZ",
              "Z": r"Z/$\gamma^*$",
              "DY": r"Z/$\gamma^*$",
              "HWW": r"H$\rightarrow$WW",
              "Rare": "Rare SM",
              "diboson": "Diboson & Rare SM",
              "VVV": "VVV",
              "wjets": "Non-prompt",
              "fake": "Non-prompt",
              "FS": "Flavor symmetric"
             }
