#include "KolmogorovSmirnovTestStat.h"

#include "RooStats/RooStatsUtils.h"
#include "TMath.h"
#include <utility>

using namespace std;

typedef pair<Double_t, Double_t> double_pair;

bool sort_pairs (double_pair i, double_pair j) {return (i.first < j.first);}

Double_t RooStats::KolmogorovSmirnovTestStat::EvaluateKSDistance( RooAbsData& data, RooArgSet& params) {
    /*
    Find the Kolmogorov-Smirnov difference between fPdf and the data. This is just
    the maximum over data points of the distance between the CDF and the emprical
    distribution function of the data.
    */

    if (!&data) {
        cout << "problem with data" << endl;
        return 0;
    }

    // get constraint parameters
    RooArgSet* constraints = fPdf->getParameters(data);
    RooStats::RemoveConstantParameters(constraints);

    //First, find the best fit values
    //TODO Can probably save this if the data doesn't change
    RooFitResult* bestFit = fPdf->fitTo(data, RooFit::Constrain(*constraints), RooFit::Save(kTRUE));

    RooRealVar* x = (RooRealVar*)params.first(); // The actualy observable variable


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

    Double_t max_distance = 0.;
    Double_t current_cdf_val = 0.;
    Double_t current_distance = 0.;
    Double_t empirical_df = 0.;

    //CDF of the PDF
    RooAbsReal* cdf;

    for(vector<double_pair>::const_iterator d = data_points.begin();
        d != data_points.end(); d++) {

        x->setVal(d->first);
        cdf = fPdf->createCdf(*x, RooFit::ScanAllCdf());
        current_cdf_val = cdf->getVal();
        empirical_df += d->second/s_data;
        current_distance = TMath::Abs(current_cdf_val-empirical_df);
        if (current_distance >  max_distance)
            max_distance = current_distance;
    }

    return max_distance;
}
