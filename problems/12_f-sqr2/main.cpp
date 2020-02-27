#include "tubex.h"
#include "tubex-solve.h"
#include <iomanip>

using namespace std;
 using namespace ibex;
using namespace tubex;

/* ode : y=-x^2 */
/* problem defined in Deville,Janssen and Van Hentenryck paper */

void contract(TubeVector& x)
{
  tubex::Function f("x", "-x^2");

  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 50000.0)
    ctc_picard.contract(f, x, FORWARD | BACKWARD );

  CtcDeriv ctc_deriv;
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.contract(x, f.eval_vector(x), FORWARD | BACKWARD);
}

int main()
{
  /* =========== PARAMETERS =========== */

    Tube::enable_syntheses(false);
  
    //    Interval domain(0.,5.0);
    //    TubeVector x(domain, 0.1, 1);
    IntervalVector v(1);
    vector<IntervalVector*> gates; 
    v[0]=Interval(0.1,0.4);
    double volume=0.0;
    double totaltime=0.;
    double step=0.5;
    int nbsteps=10;
    for (int i=0; i< nbsteps; i++){
  /* =========== SOLVER =========== */
      //      Vector epsilon(1,1.0);
      Vector epsilon(1,100.0);

      double t0=i*step;
      double t1=(i+1)*step;
      Interval domain(t0,t1);
      TubeVector x(domain, step/10, 1);
      x.set(v,t0 ); // initial condition
      tubex::Solver solver(epsilon);
      //      solver.set_refining_fxpt_ratio(0.999999);
      solver.set_refining_fxpt_ratio(0.99999);
      //    solver.set_refining_fxpt_ratio(0.98);
      //    solver.set_propa_fxpt_ratio(1.);
      solver.set_propa_fxpt_ratio(0.999);
      //      solver.set_cid_fxpt_ratio(0.9);
      solver.set_cid_fxpt_ratio(0.1);
      solver.set_cid_propa_fxpt_ratio(0.999);
      solver.set_max_slices(20000);
      //solver.set_max_slices(500);
      solver.set_cid_timept(1);
      solver.set_trace(1);
      solver.set_refining_mode(0);
      solver.set_bisection_timept(-1);

      list<TubeVector> l_solutions = solver.solve(x, &contract);
      cout << "time " << (i+1)*step <<  "nb sol " << l_solutions.size() << endl;
      if (l_solutions.size()==1) { cout << " volume " << l_solutions.front().volume() << endl;
	volume+=l_solutions.front().volume();
	totaltime+=solver.solving_time;
       v[0]=l_solutions.front()[0].last_slice()->output_gate();
       IntervalVector* v1 = new IntervalVector(v);
       gates.push_back(v1);}
      else break;
    }

    for (int k=0; k< gates.size(); k++)
      cout << k << "  " << *(gates[k]) << endl;
    cout << " total volume " << volume << " total time " << totaltime << endl;
    return 0;
}
