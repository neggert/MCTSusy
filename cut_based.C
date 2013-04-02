{
    using namespace std;
    using namespace RooFit;
    using namespace RooStats;

    // Setup
	TFile f("limits/sig_chi_600_200_combined_meas_model.root");

    RooWorkspace* wspace = (RooWorkspace*)f.Get("combined");

    RooAbsData* data = wspace->data("obsData");
    ModelConfig* model = (ModelConfig*)wspace->obj("ModelConfig");

    RooRealVar* poi = (RooRealVar*)model->GetParametersOfInterest()->first();

    RooRealVar *obs_sf = wspace->obj("obs_x_sf");
    RooRealVar *obs_of = wspace->obj("obs_x_of");

    obs_sf->setRange("fit_sf", 10., 120.);
    // obs_sf->setRange("fit_of", 10., 300.);
    // obs_of->setRange("fit_sf", 10., 300.);
    obs_of->setRange("fit_of", 10., 120.);




    // Run the fit to get some reasonable parameters
    poi->setVal(0.);
    poi->setConstant(kTRUE);

    RooArgSet constr;
    constr.add(*(model->GetNuisanceParameters()));
    RemoveConstantParameters(&constr);

    RooFitResult *res = model->GetPdf()->fitTo(*data, Constrain(constr), PrintLevel(0), Save(),
                                               Range("fit"), SplitRange());

    RooArgSet* params = const_cast<RooArgSet*>(model->GetNuisanceParameters());
    params->add(*const_cast<RooArgSet*>(model->GetParametersOfInterest()));
    model->SetSnapshot(*params);

    // SF

    // top
    RooProduct *top_shape_sf = wspace->obj("top_sf_sf_overallSyst_x_StatUncert");
    RooArgSet *top_params = top_shape_sf->getParameters(*data);
    RemoveConstantParameters(top_params)

    TF1* top_tf = top_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*top_params));

    double top_sf_low = top_tf->Integral(10, 120)/10;
    double top_sf_high = top_tf->Integral(120, 300)/10;
    double top_sf_high_err = top_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*top_params)->GetMatrixArray()) /10.;

    cout << top_sf_low << "\t" << top_sf_high << " +/- " << top_sf_high_err << endl;
    cout << top_sf_high+top_sf_low << endl;
    cout << ((RooRealVar*)wspace->obj("n_top_sf"))->getVal() << endl;

    // // Set up the test statistic and sampler
    // AndersonDarlingTestStat ad(model->GetPdf());

    // ToyMCSampler* mc = new ToyMCSampler(ad, 10);
    // mc->SetPdf(*model->GetPdf());
    // mc->SetObservables(*model->GetObservables());
    // mc->SetGlobalObservables(*model->GetObservables());
    // mc->SetNuisanceParameters(*model->GetNuisanceParameters());
    // mc->SetParametersForTestStat(*model->GetParametersOfInterest());

    // TStopwatch t;
    // // params = const_cast<RooArgSet*>(model->GetSnapshot());
    // t.Start();
    // mc->GetSamplingDistribution(model->GetSnapshot());
    // t.Stop();
    // t.Print();
}