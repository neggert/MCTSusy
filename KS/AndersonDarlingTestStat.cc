#include "AndersonDarlingTestStat.h"

#include "RooStats/RooStatsUtils.h"
#include "RooMsgService.h"
#include "TMath.h"
#include <utility>

using namespace std;

typedef pair<Double_t, Double_t> double_pair;

bool sort_pairs (double_pair i, double_pair j) {return (i.first < j.first);}

Double_t RooStats::AndersonDarlingTestStat::Evaluate( RooAbsData& data, RooArgSet& params) {
    /*
    Find the Anderson-Darling difference between fPdf and the data. This is just
    the maximum over data points of the distance between the CDF and the emprical
    distribution function of the data.
    */

    if (!&data) {
        cout << "problem with data" << endl;
        return 0;
    }

    RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    // Set POI constant so the fit doesn't adjust them
    // RooRealVar* param;
    // TIterator *param_iter = params.createIterator();
    // while (1) {
    //     param = dynamic_cast<RooRealVar*>(param_iter->Next());
    //     if (!param) break;
    //     param->setConstant(kTRUE);
    // }

    //First, find the best fit values
    RooArgSet* constraints = fFitPdf->getParameters(data);
    RooStats::RemoveConstantParameters(constraints);
    fFitPdf->fitTo(data, RooFit::Constrain(*constraints), RooFit::PrintLevel(0), RooFit::Verbose(kFALSE));

    // unet POI constant so the fit doesn't adjust them
    // param_iter = params.createIterator();
    // while (1) {
    //     param = dynamic_cast<RooRealVar*>(param_iter->Next());
    //     if (!param) break;
    //     param->setConstant(kFALSE);
    // }

    Double_t test_stat = 0;

    RooSimultaneous* pdfSim = dynamic_cast<RooSimultaneous*>(fPdf);

    // now check if we have categories to deal with
    if (pdfSim != 0) {
        const RooCategory* cat = dynamic_cast<const RooCategory*>(&(pdfSim->indexCat()));

        // split the data
        TList* datasets = data.split(*cat);
        Int_t num_cats = datasets->GetEntries();

        for(Int_t i= 0; i < num_cats; i++) {
            // set to current category

            const Text_t* cat_name = cat->lookupType(i)->GetName();

            RooAbsData *cat_data = dynamic_cast<RooAbsData*>(datasets->At(i));
            RooAbsPdf *cat_pdf = pdfSim->getPdf(cat_name);

            RooArgSet* observables = cat_pdf->getObservables(cat_data);
            RooRealVar* x = dynamic_cast<RooRealVar*>(observables->first());

            test_stat += EvaluateADDistance(*cat_pdf, *cat_data, *x);

        }
    }
    else {
        RooRealVar* obs = dynamic_cast<RooRealVar*>(fPdf->getObservables(data)->first());
        test_stat = EvaluateADDistance(*fPdf, data, *obs);
    }

    RooMsgService::instance().setGlobalKillBelow(msglevel);

    return test_stat;
}

Double_t RooStats::AndersonDarlingTestStat::EvaluateADDistance( RooAbsPdf& pdf, RooAbsData& data, RooRealVar& observable) {

    vector<double_pair > data_points;
    Int_t n_data = data.numEntries();
    Double_t s_data = data.sumEntries();

    const RooArgSet* datavals;
    RooRealVar* observable_val;

    for (int i = 0; i < n_data; i++) {
        datavals = data.get(i);
        observable_val = (RooRealVar*)(datavals->find(observable.GetName()));
        data_points.push_back(double_pair(observable_val->getVal(), data.weight()));
    }

    sort(data_points.begin(), data_points.end(), sort_pairs);

    Double_t test_stat = 0.;
    Double_t current_cdf_val = 0.;
    Double_t last_cdf_val = 0.;
    Double_t empirical_df = 0.;
    Double_t observableval, bin_prob;

    //CDF of the PDF
    RooAbsReal* cdf;

    for(vector<double_pair>::const_iterator d = data_points.begin();
        d != data_points.end()-1; d++) {
        observableval = ((d+1)->first + d->first)/2.; // d->first is middle of bin, want upper edge.
        observable.setVal(observableval);
        cdf = pdf.createCdf(observable, RooFit::ScanAllCdf());
        current_cdf_val = cdf->getVal();
        empirical_df += d->second/s_data;

        if (current_cdf_val >= 1.0)
            break;

        bin_prob = current_cdf_val-last_cdf_val;

        // from L. Demortier, CDF/ANAL/JET/CDFR/3419
        test_stat += s_data*pow((empirical_df-current_cdf_val), 2)/current_cdf_val/(1.-current_cdf_val)*bin_prob;
        last_cdf_val = current_cdf_val;
    }

    return test_stat;
}
