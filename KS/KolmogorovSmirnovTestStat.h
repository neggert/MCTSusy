#define ROOSTATS_KolmogorovSmirnovTestStat

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#include "RooStats/TestStatistic.h"

#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

namespace RooStats {
    class KolmogorovSmirnovTestStat : public TestStatistic {
      public:
        KolmogorovSmirnovTestStat() {};
        KolmogorovSmirnovTestStat(RooAbsPdf& pdf) {
            fPdf = &pdf;
        }

        ~KolmogorovSmirnovTestStat() {}

        Double_t Evaluate(RooAbsData& data, RooArgSet& params) {
            return EvaluateKSDistance(data, params);
        }

        const TString GetVarName() const { return TString("KolmogorovSmirnovTestStat");}

        bool PValueIsRightTail(void) const { return true; }

        Double_t EvaluateKSDistance(RooAbsData& data, RooArgSet& params);

      private:
        RooAbsPdf* fPdf;
    };
}
