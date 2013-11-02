CONFIG=config/data.py config/parameters.py

# fixme
NO_SIG_MODEL_CHI=limits/chi_templates_600_200_combined_meas_model.root
NO_SIG_MODEL_SLEP=limits/slep_templates_325_225_combined_meas_model.root

#fixme
no_sig_models: $(NO_SIG_MODEL_SLEP) $(NO_SIG_MODEL_CHI)

money: $(CONFIG) money_take2.py money_take2_2.py no_sig_models
	python money_take2.py $(NO_SIG_MODEL_CHI)
	python money_take2_2.py $(NO_SIG_MODEL_SLEP)