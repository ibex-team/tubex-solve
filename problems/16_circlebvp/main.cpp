/* Example circlebvp */

#include "math.h"
#include "tubex.h"
#include "tubex-solve.h"
#include <iomanip>

using namespace std;
 using namespace ibex;
using namespace tubex;

void contract(TubeVector& x)
{
  tubex::Function f("x1", "x2" ,"(-x2;x1)");

  
  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 1.e100)
    ctc_picard.contract(f, x, FORWARD | BACKWARD);
  TubeVector v = f.eval_vector(x);
  
  CtcDeriv ctc_deriv;
 
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.contract(x, v, FORWARD | BACKWARD);
  v = f.eval_vector(x);
  
  //  CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f1);     
  /*
  CtcDynCid* ctc_dyncid = new CtcDynCid(f1);     
  ctc_dyncid->set_fast_mode(true);
  CtcIntegration ctc_integration(f1,ctc_dyncid);

  ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;

  ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;

  delete ctc_dyncid;
  */
  
}

int main()
{
  tubex::Function f("x1", "x2" ,"(-x2;x1)");
  /* =========== PARAMETERS =========== */
  double pi=M_PI;
    Tube::enable_syntheses(false);

    Interval domain(0.,pi);
    //    TubeVector x(domain, 0.005, 2);
    TubeVector x(domain, pi, 2);
    IntervalVector v(2);
    v[0]=Interval(0.,0.);
    v[1]=Interval(-1.e8,1.e8);
    x.set(v, 0.); // ini
    v[0]=Interval(-1.e8,1.e8);
    v[1]=Interval(-1,-1);
    x.set(v,pi);

    double eps=0.001;
  /* =========== SOLVER =========== */
    Vector epsilon(2, eps);
    tubex::Solver solver(epsilon);

    //solver.set_refining_fxpt_ratio(0.999);
    solver.set_refining_fxpt_ratio(2.0);

    solver.set_propa_fxpt_ratio(0.999);
    // solver.set_propa_fxpt_ratio(0.);

    solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.999);

    solver.set_var3b_timept(0);
    solver.set_bisection_timept(3);

    solver.set_trace(1);
    solver.set_max_slices(40000);
    solver.set_refining_mode(0);
    solver.set_contraction_mode(4);
    list<TubeVector> l_solutions = solver.solve(x, f);
    cout << "nb sol " << l_solutions.size() << endl;
    return 0;
}





