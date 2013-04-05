using namespace std;
using namespace RooFit;
using namespace RooStats;


void cut_based() {
    // Setup
	TFile f("limits/sig_chi_600_200_combined_meas_model.root");

    RooWorkspace* wspace = (RooWorkspace*)f.Get("combined");

    RooAbsData* data = wspace->data("obsData");
    ModelConfig* model = (ModelConfig*)wspace->obj("ModelConfig");

    RooRealVar* poi = (RooRealVar*)model->GetParametersOfInterest()->first();

    RooRealVar *obs_sf = wspace->obj("obs_x_sf");
    RooRealVar *obs_of = wspace->obj("obs_x_of");

    obs_sf->setRange("fitRange", 10., 120.);
    obs_sf->setRange(10., 300.);
    obs_of->setRange("fitRange", 10., 120.);
    obs_of->setRange(10., 300.);




    // Run the fit to get some reasonable parameters
    poi->setVal(0.);
    poi->setConstant(kTRUE);

    RooArgSet constr;
    constr.add(*(model->GetNuisanceParameters()));
    RemoveConstantParameters(&constr);

    RooSimultaneous *pdf = model->GetPdf();

    RooFitResult *res = model->GetPdf()->fitTo(*data, Constrain(constr), PrintLevel(0), Save(), InitialHesse(), Minos(),
                                               Range("fitRange"));

    RooArgSet* params = const_cast<RooArgSet*>(model->GetNuisanceParameters());
    params->add(*const_cast<RooArgSet*>(model->GetParametersOfInterest()));
    model->SetSnapshot(*params);

    NumEventsTestStat dummy(*model->GetPdf());

    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    ToyMCSampler* mc = new ToyMCSampler(dummy, 1000);
    mc->SetPdf(*pdf);
    mc->SetObservables(*model->GetObservables());
    // mc->SetGlobalObservables(*model->GetObservables());
    mc->SetNuisanceParameters(*model->GetNuisanceParameters());
    mc->SetParametersForTestStat(*model->GetParametersOfInterest());
    mc->SetNEventsPerToy(0); // draw from poisson
    mc->SetGenerateBinned();

    get_results(wspace, res, data);

    vector<double> vals;

    RooAbsData* toy_data;
    for (int i = 0; i < 20; i++) {
        toy_data = mc->GenerateToyData(*model->GetSnapshot());
        toy_data->Print();
        res = pdf->fitTo(*toy_data, Constrain(constr), PrintLevel(-1), Save(), InitialHesse(), Minos(),
                                               Range("fit"));
        vals.push_back(get_results(wspace, res, data));

    }

    for (vector<double>::iterator v = vals.begin(); v != vals.end(); v++) {
        cout << *v << endl;
    }
}



double get_results(RooWorkspace* wspace, RooFitResult* res, RooAbsData* data) {

    RooRealVar *obs_sf = wspace->obj("obs_x_sf");
    RooRealVar *obs_of = wspace->obj("obs_x_of");

    // ModelConfig *model = (ModelConfig*)wspace->obj("ModelConfig");
    // model->LoadSnapshot();

    // SF
    cout << "Same-Flavor" << endl;
    cout << "===========" << endl;

    // sum

    RooProdPdf *sum_shape_sf = wspace->obj("sf_model");
    RooArgSet *sum_params_sf = sum_shape_sf->getParameters(*data);
    RemoveConstantParameters(sum_params_sf);

    TF1* sum_tf = sum_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*sum_params_sf));

    double sum_sf_low = sum_tf->Integral(10, 120);
    double sum_sf_low_err = sum_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*sum_params_sf).GetMatrixArray());
    double sum_sf_high = sum_tf->Integral(120, 300);
    double sum_sf_high_err = sum_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*sum_params_sf).GetMatrixArray());

    cout << "Total: " << sum_sf_low << " +/- " << sum_sf_low_err << "\t\t" << sum_sf_high << " +/- " << sum_sf_high_err << endl;

    // top
    RooProduct *top_shape_sf = wspace->obj("top_sf_sf_overallSyst_x_StatUncert");
    RooArgSet *top_params = top_shape_sf->getParameters(*data);
    RemoveConstantParameters(top_params);

    TF1* top_tf = top_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*top_params));

    double top_sf_low = top_tf->Integral(10, 120)/10;
    double top_sf_low_err = top_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*top_params).GetMatrixArray()) /10.;
    double top_sf_high = top_tf->Integral(120, 300)/10;
    double top_sf_high_err = top_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*top_params).GetMatrixArray()) /10.;

    cout << "Top: " << top_sf_low << " +/- " << top_sf_low_err << "\t\t" << top_sf_high << " +/- " << top_sf_high_err << endl;


    // Diboson
    RooProduct *vv_shape_sf = wspace->obj("vv_sf_sf_overallSyst_x_StatUncert_x_sf_ww_syst_sf_ShapeSys");
    RooArgSet *vv_params = vv_shape_sf->getParameters(*data);
    RemoveConstantParameters(vv_params);
    TF1* vv_tf = vv_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*vv_params));

    double vv_sf_low = vv_tf->Integral(10, 120)/10;
    double vv_sf_low_err = vv_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*vv_params).GetMatrixArray()) /10.;
    double vv_sf_high = vv_tf->Integral(120, 300)/10;
    double vv_sf_high_err = vv_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*vv_params).GetMatrixArray()) /10.;

    cout << "VV: " << vv_sf_low << " +/- " << vv_sf_low_err << "\t\t" << vv_sf_high << " +/- " << vv_sf_high_err << endl;

    // Z
    RooProduct *z_shape_sf = wspace->obj("z_sf_sf_overallSyst_x_StatUncert_x_sf_z_syst_sf_ShapeSys");
    RooArgSet *z_params = z_shape_sf->getParameters(*data);
    RemoveConstantParameters(z_params);
    TF1* z_tf = z_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*z_params));

    double z_sf_low = z_tf->Integral(10, 120)/10;
    double z_sf_low_err = z_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*z_params).GetMatrixArray()) /10.;
    double z_sf_high = z_tf->Integral(120, 300)/10;
    double z_sf_high_err = z_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*z_params).GetMatrixArray()) /10.;

    cout << "z: " << z_sf_low << " +/- " << z_sf_low_err << "\t\t" << z_sf_high << " +/- " << z_sf_high_err << endl;

    // W+Jets
    RooProduct *fake_shape_sf = wspace->obj("wjets_sf_sf_overallSyst_x_StatUncert_x_sf_wjets_syst_sf_ShapeSys");
    RooArgSet *fake_params = fake_shape_sf->getParameters(*data);
    RemoveConstantParameters(fake_params);
    TF1* fake_tf = fake_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*fake_params));

    double fake_sf_low = fake_tf->Integral(10, 120)/10;
    double fake_sf_low_err = fake_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*fake_params).GetMatrixArray()) /10.;
    double fake_sf_high = fake_tf->Integral(120, 300)/10;
    double fake_sf_high_err = fake_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*fake_params).GetMatrixArray()) /10.;

    cout << "fake: " << fake_sf_low << " +/- " << fake_sf_low_err << "\t\t" << fake_sf_high << " +/- " << fake_sf_high_err << endl;

    // OF
    cout << "Opposite-Flavor" << endl;
    cout << "===========" << endl;

    // sum

    RooProdPdf *sum_shape_of = wspace->obj("of_model");
    RooArgSet *sum_params_of = sum_shape_of->getParameters(*obs_of);
    RemoveConstantParameters(sum_params_of);

    TF1* sum_tf = sum_shape_of->asTF(RooArgList(*obs_of), *sum_params_of);

    double sum_of_low = sum_tf->Integral(10, 120);
    double sum_of_low_err = sum_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*sum_params_of).GetMatrixArray());
    double sum_of_high = sum_tf->Integral(120, 300);
    double sum_of_high_err = sum_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*sum_params_of).GetMatrixArray());

    cout << "Total: " << sum_of_low << " +/- " << sum_of_low_err << "\t\t" << sum_of_high << " +/- " << sum_of_high_err << endl;

    // top
    RooProduct *top_shape_of = wspace->obj("top_of_of_overallSyst_x_StatUncert");

    // RooRealIntegral *i = (RooRealIntegral*) top_shape_of->createIntegral(*obs_of, "fitRange");
    // cout << i->getVal() << endl;

    RooArgSet *top_params = top_shape_of->getParameters(*data);
    RemoveConstantParameters(top_params);

    TF1* top_tf = top_shape_of->asTF(RooArgList(*obs_of), RooArgList(*top_params));

    double top_of_low = top_tf->Integral(10, 120)/10;
    double top_of_low_err = top_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*top_params).GetMatrixArray()) /10.;
    double top_of_high = top_tf->Integral(120, 300)/10;
    double top_of_high_err = top_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*top_params).GetMatrixArray()) /10.;

    cout << "Top: " << top_of_low << " +/- " << top_of_low_err << "\t\t" << top_of_high << " +/- " << top_of_high_err << endl;


    // Diboson
    RooProduct *vv_shape_of = wspace->obj("vv_of_of_overallSyst_x_StatUncert_x_of_ww_syst_of_ShapeSys");
    RooArgSet *vv_params = vv_shape_of->getParameters(*data);
    RemoveConstantParameters(vv_params);
    TF1* vv_tf = vv_shape_of->asTF(RooArgList(*obs_of), RooArgList(*vv_params));

    double vv_of_low = vv_tf->Integral(10, 120)/10;
    double vv_of_low_err = vv_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*vv_params).GetMatrixArray()) /10.;
    double vv_of_high = vv_tf->Integral(120, 300)/10;
    double vv_of_high_err = vv_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*vv_params).GetMatrixArray()) /10.;

    cout << "VV: " << vv_of_low << " +/- " << vv_of_low_err << "\t\t" << vv_of_high << " +/- " << vv_of_high_err << endl;

    // Z
    RooProduct *z_shape_of = wspace->obj("z_of_of_overallSyst_x_StatUncert_x_of_z_syst_of_ShapeSys");
    RooArgSet *z_params = z_shape_of->getParameters(*data);
    RemoveConstantParameters(z_params);
    TF1* z_tf = z_shape_of->asTF(RooArgList(*obs_of), RooArgList(*z_params));

    double z_of_low = z_tf->Integral(10, 120)/10;
    double z_of_low_err = z_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*z_params).GetMatrixArray()) /10.;
    double z_of_high = z_tf->Integral(120, 300)/10;
    double z_of_high_err = z_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*z_params).GetMatrixArray()) /10.;

    cout << "z: " << z_of_low << " +/- " << z_of_low_err << "\t\t" << z_of_high << " +/- " << z_of_high_err << endl;

    // W+Jets
    RooProduct *fake_shape_of = wspace->obj("wjets_of_of_overallSyst_x_StatUncert_x_of_wjets_syst_of_ShapeSys");
    RooArgSet *fake_params = fake_shape_of->getParameters(*data);
    RemoveConstantParameters(fake_params);
    TF1* fake_tf = fake_shape_of->asTF(RooArgList(*obs_of), RooArgList(*fake_params));

    double fake_of_low = fake_tf->Integral(10, 120)/10;
    double fake_of_low_err = fake_tf->IntegralError(10, 120, 0, res->reducedCovarianceMatrix(*fake_params).GetMatrixArray()) /10.;
    double fake_of_high = fake_tf->Integral(120, 300)/10;
    double fake_of_high_err = fake_tf->IntegralError(120, 300, 0, res->reducedCovarianceMatrix(*fake_params).GetMatrixArray()) /10.;

    cout << "fake: " << fake_of_low << " +/- " << fake_of_low_err << "\t\t" << fake_of_high << " +/- " << fake_of_high_err << endl;


    return sum_sf_high;

}


