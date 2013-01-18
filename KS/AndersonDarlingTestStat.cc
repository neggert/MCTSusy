#include "AndersonDarlingTestStat.h"

#include "RooStats/RooStatsUtils.h"
#include "RooMsgService.h"
#include "TMath.h"
#include <utility>

using namespace std;

typedef pair<Double_t, Double_t> double_pair;

bool sort_pairs (double_pair i, double_pair j) {return (i.first < j.first);}

Double_t RooStats::AndersonDarlingTestStat::EvaluateADDistance( RooAbsData& data, RooArgSet& params) {
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

    // get constraint parameters
    RooArgSet* constraints = fPdf->getParameters(data);
    RooStats::RemoveConstantParameters(constraints);

    //First, find the best fit values
    //TODO Can probably save this if the data doesn't change
    fPdf->fitTo(data, RooFit::Constrain(*constraints), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE));


    RooArgSet* observables = fPdf->getObservables(data);
    RooRealVar* x = (RooRealVar*)observables->first(); // The actualy observable variable


    vector<double_pair > data_points;
    Int_t n_data = data.numEntries();
    Double_t s_data = data.sumEntries();

    const RooArgSet* datavals;
    RooRealVar* x_val;

    for (int i = 0; i < n_data; i++) {
        datavals = data.get(i);
        x_val = (RooRealVar*)(datavals->find(x->GetName()));
        data_points.push_back(double_pair(x_val->getVal(), data.weight()));
    }

    sort(data_points.begin(), data_points.end(), sort_pairs);

    Double_t test_stat = 0.;
    Double_t current_cdf_val = 0.;
    Double_t last_cdf_val = 0.;
    Double_t empirical_df = 0.;
    Double_t xval, bin_prob;

    //CDF of the PDF
    RooAbsReal* cdf;

    for(vector<double_pair>::const_iterator d = data_points.begin();
        d != data_points.end()-1; d++) {
        xval = ((d+1)->first + d->first)/2.; // d->first is middle of bin, want upper edge.
        x->setVal(xval);
        cdf = fPdf->createCdf(*x, RooFit::ScanAllCdf());
        current_cdf_val = cdf->getVal();
        empirical_df += d->second/s_data;

        bin_prob = current_cdf_val-last_cdf_val;

        // from L. Demortier, CDF/ANAL/JET/CDFR/3419
        test_stat += s_data*pow((empirical_df-current_cdf_val), 2)/current_cdf_val/(1.-current_cdf_val)*bin_prob;
        last_cdf_val = current_cdf_val;
    }

    RooMsgService::instance().setGlobalKillBelow(msglevel);

    return test_stat;
}
