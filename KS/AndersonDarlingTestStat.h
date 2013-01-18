#define ROOSTATS_AndersonDarlingTestStat

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
    class AndersonDarlingTestStat : public TestStatistic {
      public:
        AndersonDarlingTestStat() {};
        AndersonDarlingTestStat(RooAbsPdf& pdf) {
            fPdf = &pdf;
        }

        ~AndersonDarlingTestStat() {}

        Double_t Evaluate(RooAbsData& data, RooArgSet& params) {
            return EvaluateADDistance(data, params);
        }

        const TString GetVarName() const { return TString("Anderson-Darling Statistic");}

        bool PValueIsRightTail(void) const { return true; }

        Double_t EvaluateADDistance(RooAbsData& data, RooArgSet& params);

      private:
        RooAbsPdf* fPdf;
    };
}
