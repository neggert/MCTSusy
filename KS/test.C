#include "KolmogorovSmirnovTestStat.h"
#include "AndersonDarlingTestStat.h"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "TFile.h"
#include "TCanvas.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

int test() {

    TFile f("../limits/sig_chi_250_50_combined_meas_model.root");

    RooWorkspace* w = (RooWorkspace*)f.Get("combined");
    ModelConfig* model = (ModelConfig*)w->obj("ModelConfig");

    ((RooRealVar*)model->GetParametersOfInterest()->first())->setVal(0.);
    ((RooRealVar*)model->GetParametersOfInterest()->first())->setConstant(kTRUE);

    RooAbsPdf* pdf = model->GetPdf();

    RooAbsData* data = w->data("obsData");

    RooArgSet obs = *(model->GetObservables());
    RooArgSet poi = *(model->GetParametersOfInterest());


    //First, find the best fit values
    RooArgSet* constraints = pdf->getParameters(data);
    RooStats::RemoveConstantParameters(constraints);
    pdf->fitTo(*data, RooFit::Constrain(*constraints));
    model->SetSnapshot(*(model->GetParametersOfInterest()));

    //Get the AD test statistic for the fit to the data
    AndersonDarlingTestStat ad(*pdf);
    Double_t ts = ad.Evaluate(*data, poi);

    ToyMCSampler* toy = new ToyMCSampler (ad, 100);
    toy->SetPdf(*pdf);
    toy->SetObservables(obs);
    toy->SetGlobalObservables(*(model->GetGlobalObservables()));
    toy->SetParametersForTestStat(*(model->GetParametersOfInterest()));

    RooArgSet params;
    params.add(*(model->GetNuisanceParameters()));
    params.add(*(model->GetParametersOfInterest()));

    toy->GenerateToyData(params, *pdf)->Print();

    SamplingDistribution* sampDist = toy->GetSamplingDistribution(params);

    SamplingDistPlot plot(10);
    plot.AddSamplingDistribution(sampDist);

    TCanvas c1;
    plot.Draw();
    c1.SaveAs("test.pdf");

    Double_t p_value = 1.-sampDist->CDF(ts);

    cout << "Test statistic on data: " << ts << endl;

    cout << "P-value: " << p_value << endl;

    return 0;

}
