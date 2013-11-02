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
#include "RooCategory.h"
#include "RooSimultaneous.h"

namespace RooStats {
    class AndersonDarlingTestStat : public TestStatistic {
      public:
        AndersonDarlingTestStat(): fPdf(NULL) { fVarName = "Anderson-Darling Test Statistic"; }
        AndersonDarlingTestStat(RooAbsPdf& pdf) {
            fPdf = &pdf;
            fVarName = "Anderson-Darling Test Statistic";
        }

        ~AndersonDarlingTestStat() {}

        Double_t Evaluate(RooAbsData& data, RooArgSet& paramsOfInterest);

        Double_t Evaluate(RooAbsData& data);

        
        virtual void SetVarName(const char* name) { fVarName = name; }
        virtual const TString GetVarName() const {return fVarName;}


        bool PValueIsRightTail(void) const { return true; }

        Double_t EvaluateADDistance(RooAbsPdf& pdf, RooAbsData& data, RooRealVar& observable);

      private:
        RooAbsPdf* fPdf;
        TString fVarName;

      protected:
        ClassDef(AndersonDarlingTestStat,1)
    };
}
