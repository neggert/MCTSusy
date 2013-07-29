%module mt2_bisect

%{
	#include "mt2_bisect.h"
%}

%typemap(in) double*(double temp[3]) {   // temp[3] becomes a local variable
  int i;
  if (PyTuple_Check($input)) {
    if (!PyArg_ParseTuple($input,"ddd",temp,temp+1,temp+2)) {
      PyErr_SetString(PyExc_TypeError,"tuple must have 3 elements");
      return NULL;
    }
    $1 = &temp[0];
  } else {
    PyErr_SetString(PyExc_TypeError,"expected a tuple.");
    return NULL;
  }
}

namespace mt2_bisect
{
class mt2
{  
   public:

      mt2();
      void   mt2_bisect();
      void   mt2_massless();
      void   set_momenta(double *pa0, double *pb0, double* pmiss0);
      void   set_mn(double mn);
      double get_mt2();
      void   print();
      int    nevt;
   private:  

      bool   solved;
      bool   momenta_set;
      double mt2_b;

      int    nsols(double Dsq);
      int    nsols_massless(double Dsq);
      inline int    signchange_n( long double t1, long double t2, long double t3, long double t4, long double t5);
      inline int    signchange_p( long double t1, long double t2, long double t3, long double t4, long double t5);
      int scan_high(double &Deltasq_high);
      int find_high(double &Deltasq_high);
      //data members
      double pax, pay, ma, Ea;
      double pmissx, pmissy;
      double pbx, pby, mb, Eb;
      double mn, mn_unscale;
     
      //auxiliary definitions
      double masq, Easq;
      double mbsq, Ebsq;
      double pmissxsq, pmissysq;
      double mnsq;

      //auxiliary coefficients
      double a1, b1, c1, a2, b2, c2, d1, e1, f1, d2, e2, f2;
      double d11, e11, f12, f10, d21, d20, e21, e20, f22, f21, f20;

      double scale;
      double precision;
};

}