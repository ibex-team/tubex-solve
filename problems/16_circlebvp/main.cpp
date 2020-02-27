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
  if (x.volume() > 50000.)
    ctc_picard.contract(f, x, FORWARD | BACKWARD);

  CtcDeriv ctc_deriv;
 
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.contract(x, f.eval_vector(x), FORWARD | BACKWARD);
}

int main()
{
  /* =========== PARAMETERS =========== */
  double pi=M_PI;
    Tube::enable_syntheses(false);

    Interval domain(0.,pi);
    TubeVector x(domain, 0.005, 2);
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

    solver.set_refining_fxpt_ratio(0.9999);

    solver.set_propa_fxpt_ratio(0.999);

    solver.set_cid_fxpt_ratio(0.9999);

    solver.set_cid_propa_fxpt_ratio(0.999);
    solver.set_cid_timept(2);
    solver.set_bisection_timept(2);
    solver.set_trace(1);
    solver.set_max_slices(20000);
    solver.set_refining_mode(2);
    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << "nb sol " << l_solutions.size() << endl;
    return 0;
}





