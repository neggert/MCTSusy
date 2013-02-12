#include "TH2D.h"

void repareHoles (TH2D* hist) {
  for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
      if (hist->GetBinContent (ix,iy) <=0) {
        if (ix > 1 && ix < hist->GetNbinsX() &&
            hist->GetBinContent (ix-1,iy) >1e-40 &&
            hist->GetBinContent (ix+1,iy) >1e-40) {
          hist->SetBinContent (ix,iy, 0.5*(hist->GetBinContent (ix-1,iy)+hist->GetBinContent (ix+1,iy)));
        }
      }
    }
  }

  for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
      if (hist->GetBinContent (ix,iy) <=0) {
        if (iy > 1 && iy < hist->GetNbinsY() &&
            hist->GetBinContent (ix,iy-1) >1e-40 &&
            hist->GetBinContent (ix,iy+1) >1e-40) {
          hist->SetBinContent (ix,iy, 0.5*(hist->GetBinContent (ix,iy-1)+hist->GetBinContent (ix,iy+1)));
        }
      }
    }
  }
  bool thisDone = true;
  while (thisDone) {
    thisDone = false;
    for (int iy = 1; !thisDone && iy <= hist->GetNbinsY(); ++iy) {
      for (int ix = 1; !thisDone && ix <= hist->GetNbinsX(); ++ix) {
        if (hist->GetBinContent (ix,iy) <=0) {
          if (ix > 1 && iy > 1) {
            if (hist->GetBinContent (ix-1,iy-1) > 1e-40 &&
                hist->GetBinContent (ix,iy-1) > 1e-40 &&
                hist->GetBinContent (ix-1,iy) > 1e-40) {
              hist->SetBinContent (ix,iy,
                                   hist->GetBinContent (ix,iy-1)*
                                   hist->GetBinContent (ix-1,iy)/
                                   hist->GetBinContent (ix-1,iy-1));
              thisDone = true;
              break;


            }
          }
        }
      }
    }
  }
}
