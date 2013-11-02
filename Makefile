CONFIG=config/data.py config/parameters.py

# fixme
NO_SIG_MODEL_CHI=limits/chi_templates_600_200_combined_meas_model.root
NO_SIG_MODEL_SLEP=limits/slep_templates_325_225_combined_meas_model.root

MONEY_REQS=money_take2.py diboson_fracs.json chi_templates.root config/parameters.py

MONEY2_REQS=money_take2_2.py diboson_fracs2.json slep_templates.root config/parameters.py

fit_results.json: $(NO_SIG_MODEL_CHI) $(MONEY2_REQS)
	python bkg_fit.py $(NO_SIG_MODEL_CHI)

fit_results2.json: $(NO_SIG_MODEL_SLEP) $(MONEY2_REQS)
	python bkg_fit2.py $(NO_SIG_MODEL_SLEP)

bkg_fit: fit_results.json
bkg_fit2: fit_results2.json
fit: bkg_fit bkg_fit2
moneyplots: fit