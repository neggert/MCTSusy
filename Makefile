CONFIG=config/data.py config/parameters.py

# fixme
NO_SIG_MODEL_CHI=limits/chi_bkg_only_combined_meas_model.root
NO_SIG_MODEL_SLEP=limits/slep_bkg_only_combined_meas_model.root

AD=KS/AndersonDarlingTestStat.cc KS/AndersonDarlingTestStat.h

$(NO_SIG_MODEL_CHI): histfactory.py templates.root data.root
	./histfactory.py templates.root

$(NO_SIG_MODEL_SLEP): histfactory2.py templates2.root data.root
	./histfactory2.py templates2.root

chi_templates.root templates.root data.root: prep_hists.py
	./prep_hists.py work/sms/sms_chi_25GeV.hdf5 chi_templates.root limits/8TeVc1c1.xsec limits/TChipmSlepSnu_Nevents.root --low=10 --high=200 --bins=19

MONEY_REQS=money_take2.py diboson_fracs.json chi_templates.root config/parameters.py

MONEY2_REQS=money_take2_2.py diboson_fracs2.json slep_templates.root config/parameters.py

fit_results.json: $(NO_SIG_MODEL_CHI) $(MONEY2_REQS) $(AD) IntegralError.C
	python bkg_fit.py $(NO_SIG_MODEL_CHI)

fit_results2.json: $(NO_SIG_MODEL_SLEP) $(MONEY2_REQS) $(AD) IntegralError.C
	python bkg_fit2.py $(NO_SIG_MODEL_SLEP)

p-value: $(NO_SIG_MODEL_CHI) $(MONEY2_REQS) $(AD)
	ipcluster start -n 8 --daemonize
	python bkg_fit.py $(NO_SIG_MODEL_CHI) -p
	ipcluster stop

mc_preds.json: $(CONFIG) mc_preds.py selection.py
	python mc_preds.py

table: fit_results.json mc_preds.json print_results_table.py
	python print_results_table.py

bkg_fit: fit_results.json
bkg_fit2: fit_results2.json
moneyplots: fit