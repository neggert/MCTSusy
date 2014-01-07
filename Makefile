CONFIG=config/data.py config/parameters.py

NO_SIG_MODEL_CHI=limits/chi_bkg_only_combined_meas_model.root
NO_SIG_MODEL_CHI_CONSTRAINED=limits/chi_bkg_only_constrained.root
NO_SIG_MODEL_SLEP=limits/slep_bkg_only_combined_meas_model.root

PMSSM_MODELS=$(shell python list_models.py work/pMSSM.hdf5 pMSSM_templates)
SIG_MODELS_CHI=$(shell python list_models.py work/sms/sms_chi_201311.hdf5 chi_templates)
SIG_MODELS_SLEP=$(shell python list_models.py work/sms/sms_slep_201311.hdf5 slep_templates)


AD=KS/AndersonDarlingTestStat.cc KS/AndersonDarlingTestStat.h

# models
$(NO_SIG_MODEL_CHI): histfactory.py templates.root data.root
	./histfactory.py templates.root

$(NO_SIG_MODEL_CHI_CONSTRAINED): histfactory.py templates.root data.root
	./histfactory.py templates.root --constrain

$(NO_SIG_MODEL_SLEP): histfactory2.py templates2.root data.root
	./histfactory2.py templates2.root

$(SIG_MODELS_CHI): histfactory.py templates.root data.root chi_templates.root
	./histfactory.py templates.root chi_templates.root $@ > /dev/null

$(SIG_MODELS_SLEP): histfactory2.py templates2.root data.root slep_templates.root
	./histfactory2.py templates2.root slep_templates.root $@ > /dev/null

$(PMSSM_MODELS): histfactory.py pMSSM_templates.root templates.root
	./histfactory.py templates.root pMSSM_templates.root $@ > /dev/null

# misc
diboson_fracs.json: one-offs/diboson_fracs.py $(CONFIG)
	python one-offs/diboson_fracs.py

# templates
chi_templates.root templates.root data.root: prep_hists.py $(CONFIG)
	./prep_hists.py work/sms/sms_chi_201311.hdf5 chi_templates.root limits/8TeVc1c1.xsec limits/TChipmSlepSnu_Nevents.root --low=10 --high=200 --bins=19

slep_templates.root templates2.root: prep_hists2.py $(CONFIG)
	./prep_hists2.py work/sms/sms_slep_201311.hdf5 slep_templates.root limits/8TeVeLeL.xsec limits/TSlepSlep_Nevents.root --low=10 --high=300 --bins=29 -x 2

pMSSM_templates.root: prep_hists_pMSSM.py pMSSM_points.json limits/nevents_pMSSM.json work/pMSSM.hdf5
	./prep_hists_pMSSM.py work/pMSSM.hdf5 pMSSM_templates.root --low=10 --high=200 --bins=19

# results
pMSSM_likelihoods.txt: $(PMSSM_MODELS) likelihood_pMSSM.py limits/nevents_pMSSM.json
	./likelihood_pMSSM.py $(PMSSM_MODELS)

# requirements for money plots
MONEY_REQS=money_take2.py diboson_fracs.json chi_templates.root config/parameters.py
MONEY2_REQS=money_take2_2.py diboson_fracs2.json slep_templates.root config/parameters.py

fit_results.json: $(NO_SIG_MODEL_CHI) $(MONEY_REQS) $(AD) IntegralError.C bkg_fit.py
	python bkg_fit.py -m $(NO_SIG_MODEL_CHI) $@

count_fit_results.json: $(NO_SIG_MODEL_CHI_CONSTRAINED) $(MONEY_REQS) $(AD) IntegralError.C bkg_fit.py
	python bkg_fit.py $(NO_SIG_MODEL_CHI_CONSTRAINED) $@ -m --cut=120

fit_results2.json: $(NO_SIG_MODEL_SLEP) $(MONEY2_REQS) $(AD) IntegralError.C bkg_fit2.py
	python bkg_fit2.py $(NO_SIG_MODEL_SLEP) $@

count_fit_results2.json: $(NO_SIG_MODEL_SLEP) $(MONEY2_REQS) $(AD) IntegralError.C bkg_fit2.py
	python bkg_fit2.py $(NO_SIG_MODEL_SLEP) $@ --cut=120

p-value: $(NO_SIG_MODEL_CHI) $(MONEY_REQS) $(AD)
	ipcluster start -n 8 --daemonize
	python bkg_fit.py $(NO_SIG_MODEL_CHI) tmp.json -p
	rm tmp.json
	ipcluster stop

p-value2: $(NO_SIG_MODEL_SLEP) $(MONEY2_REQS) $(AD)
	ipcluster start -n 8 --daemonize
	python bkg_fit2.py $(NO_SIG_MODEL_SLEP) tmp.json -p
	rm tmp.json
	ipcluster stop


# for documentation
mc_preds.json: $(CONFIG) mc_preds.py selection.py
	python mc_preds.py

mc_preds2.json: $(CONFIG) mc_preds2.py selection.py
	python mc_preds2.py

table: fit_results.json mc_preds.json print_results_table.py
	python print_results_table.py

table2: fit_results2.json mc_preds2.json print_results_table2.py
	python print_results_table2.py

# shortcuts
models_chi: $(SIG_MODELS_CHI)
models_slep: $(SIG_MODELS_SLEP)
bkg_fit: fit_results.json
bkg_fit2: fit_results2.json
counting: count_fit_results.json
pMSSM: pMSSM_likelihoods.txt