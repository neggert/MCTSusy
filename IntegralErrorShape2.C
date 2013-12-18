#include "TF1.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooRandom.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooNumber.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

Double_t get_integral_error(TF1* f, const Double_t low, const Double_t high, TMatrixTSym<double>* cov) {
	return f->IntegralError(low, high, 0, cov->GetMatrixArray());
}

struct val_with_error {
	double val;
	double error;
};

struct background {
	val_with_error low;
	val_with_error high;
};

struct channel {
        val_with_error generated_sum;
	background sum;
	background top;
	background vv;
	background z;
	background fake;
};

struct result {
	channel of;
	channel sf;
};

channel get_results2(RooWorkspace* wspace, RooFitResult* res, float cut, float upper) {

    RooRealVar *obs_sf = (RooRealVar*)wspace->obj("obs_x_sf");

    // ModelConfig *model = (ModelConfig*)wspace->obj("ModelConfig");
    // model->LoadSnapshot();

    // SF
    cout << "Same-Flavor" << endl;
    cout << "===========" << endl;

    // sum

    RooProdPdf *sum_shape_sf = (RooProdPdf*)wspace->obj("sf_model");
    RooArgSet *sum_params_sf = sum_shape_sf->getParameters(*obs_sf);
    RemoveConstantParameters(sum_params_sf);

    TF1* sum_tf = sum_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*sum_params_sf));

    double sum_sf_low = sum_tf->Integral(10, cut);
    double sum_sf_low_err = sum_tf->IntegralError(10, cut, 0, res->reducedCovarianceMatrix(*sum_params_sf).GetMatrixArray());
    double sum_sf_high = sum_tf->Integral(cut, upper);
    double sum_sf_high_err = sum_tf->IntegralError(cut, upper, 0, res->reducedCovarianceMatrix(*sum_params_sf).GetMatrixArray());

    cout << "Total: " << sum_sf_low << " +/- " << sum_sf_low_err << "\t\t" << sum_sf_high << " +/- " << sum_sf_high_err << endl;
    val_with_error low = {sum_sf_low, sum_sf_low_err};
    val_with_error high = {sum_sf_high, sum_sf_high_err};
    background sum = {low, high};

    // top
    RooProduct *top_shape_sf = (RooProduct*)wspace->obj("of_sf_overallSyst_x_StatUncert");
    RooArgSet *top_params = top_shape_sf->getParameters(*obs_sf);
    RemoveConstantParameters(top_params);

    TF1* top_tf = top_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*top_params));

    double top_sf_low = top_tf->Integral(10, cut)/10;
    double top_sf_low_err = top_tf->IntegralError(10, cut, 0, res->reducedCovarianceMatrix(*top_params).GetMatrixArray()) /10.;
    double top_sf_high = top_tf->Integral(cut, upper)/10;
    double top_sf_high_err = top_tf->IntegralError(cut, upper, 0, res->reducedCovarianceMatrix(*top_params).GetMatrixArray()) /10.;

    cout << "Top: " << top_sf_low << " +/- " << top_sf_low_err << "\t\t" << top_sf_high << " +/- " << top_sf_high_err << endl;
    val_with_error s_top_sf_low = {top_sf_low, top_sf_low_err};
    val_with_error s_top_sf_high = {top_sf_high, top_sf_high_err};
    background top = {s_top_sf_low, s_top_sf_high};


    // Diboson
    RooProduct *vv_shape_sf = (RooProduct*)wspace->obj("vv_sf_overallSyst_x_StatUncert");
    RooArgSet *vv_params = vv_shape_sf->getParameters(*obs_sf);
    RemoveConstantParameters(vv_params);
    TF1* vv_tf = vv_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*vv_params));

    double vv_sf_low = vv_tf->Integral(10, cut)/10;
    double vv_sf_low_err = vv_tf->IntegralError(10, cut, 0, res->reducedCovarianceMatrix(*vv_params).GetMatrixArray()) /10.;
    double vv_sf_high = vv_tf->Integral(cut, upper)/10;
    double vv_sf_high_err = vv_tf->IntegralError(cut, upper, 0, res->reducedCovarianceMatrix(*vv_params).GetMatrixArray()) /10.;

    cout << "VV: " << vv_sf_low << " +/- " << vv_sf_low_err << "\t\t" << vv_sf_high << " +/- " << vv_sf_high_err << endl;
    val_with_error s_vv_sf_low = {vv_sf_low, vv_sf_low_err};
    val_with_error s_vv_sf_high = {vv_sf_high, vv_sf_high_err};
    background vv = {s_vv_sf_low, s_vv_sf_high};

    // Z
    RooProduct *z_shape_sf = (RooProduct*)wspace->obj("z_sf_overallSyst_x_StatUncert_x_sf_z_syst_ShapeSys");
    RooArgSet *z_params = z_shape_sf->getParameters(*obs_sf);
    RemoveConstantParameters(z_params);
    TF1* z_tf = z_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*z_params));

    double z_sf_low = z_tf->Integral(10, cut)/10;
    double z_sf_low_err = z_tf->IntegralError(10, cut, 0, res->reducedCovarianceMatrix(*z_params).GetMatrixArray()) /10.;
    double z_sf_high = z_tf->Integral(cut, upper)/10;
    double z_sf_high_err = z_tf->IntegralError(cut, upper, 0, res->reducedCovarianceMatrix(*z_params).GetMatrixArray()) /10.;

    cout << "z: " << z_sf_low << " +/- " << z_sf_low_err << "\t\t" << z_sf_high << " +/- " << z_sf_high_err << endl;
    val_with_error s_z_sf_low = {z_sf_low, z_sf_low_err};
    val_with_error s_z_sf_high = {z_sf_high, z_sf_high_err};
    background z = {s_z_sf_low, s_z_sf_high};

    // W+Jets
    RooProduct *fake_shape_sf = (RooProduct*)wspace->obj("wjets_sf_overallSyst_x_StatUncert_x_sf_wjets_syst_ShapeSys");
    RooArgSet *fake_params = fake_shape_sf->getParameters(*obs_sf);
    RemoveConstantParameters(fake_params);
    TF1* fake_tf = fake_shape_sf->asTF(RooArgList(*obs_sf), RooArgList(*fake_params));

    double fake_sf_low = fake_tf->Integral(10, cut)/10;
    double fake_sf_low_err = fake_tf->IntegralError(10, cut, 0, res->reducedCovarianceMatrix(*fake_params).GetMatrixArray()) /10.;
    double fake_sf_high = fake_tf->Integral(cut, upper)/10;
    double fake_sf_high_err = fake_tf->IntegralError(cut, upper, 0, res->reducedCovarianceMatrix(*fake_params).GetMatrixArray()) /10.;

    cout << "fake: " << fake_sf_low << " +/- " << fake_sf_low_err << "\t\t" << fake_sf_high << " +/- " << fake_sf_high_err << endl;
    val_with_error s_fake_sf_low = {fake_sf_low, fake_sf_low_err};
    val_with_error s_fake_sf_high = {fake_sf_high, fake_sf_high_err};
    background fake = {s_fake_sf_low, s_fake_sf_high};

    val_with_error dummy_sf = {0,0};

    channel sf = {dummy_sf, sum, top, vv, z, fake};

    return sf;

}


// void IntegralError() {
//     fit_toy("temp.root", 1000000);
// }
