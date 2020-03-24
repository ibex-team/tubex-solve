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
  ibex::Function f1("x", "-x^2");

  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(true);

  if (x.volume() > 1.e300)
    ctc_picard.contract(f, x, FORWARD );

  if (x.volume() < 1.e300){
    /*
    CtcDeriv ctc_deriv;
    ctc_deriv.set_fast_mode(true);
    ctc_deriv.contract(x, f.eval_vector(x), FORWARD | BACKWARD);
    */
    
    TubeVector v = f.eval_vector(x);
    //    CtcDynCid* ctc_dyncid = new CtcDynCid(f1);     //f2 the function
    CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f1);     //f2 the function
    ctc_dyncid->set_fast_mode(true);
    CtcIntegration ctc_integration(f1,ctc_dyncid);
    ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;
    ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;
    delete ctc_dyncid;
    

  }
   
}

int main()
{
  /* =========== PARAMETERS =========== */

    Tube::enable_syntheses(false);
  
    Interval domain(0.,5.0);
    TubeVector x(domain, 0.2, 1);
    IntervalVector v(1);
    vector<IntervalVector*> gates; 
    v[0]=Interval(0.1,0.4);
    
    Vector epsilon(1,0.2);


    
    x.set(v,0.0 ); // initial condition
    tubex::Solver solver(epsilon);
      //      solver.set_refining_fxpt_ratio(0.999999);
    solver.set_refining_fxpt_ratio(0.99999);
      //    solver.set_refining_fxpt_ratio(0.98);
      //    solver.set_propa_fxpt_ratio(1.);
    solver.set_propa_fxpt_ratio(0.999);
    //    solver.set_var3b_fxpt_ratio(0.9);
    solver.set_var3b_fxpt_ratio(0.);
    solver.set_var3b_propa_fxpt_ratio(0.999);
    solver.set_max_slices(100000);
      //solver.set_max_slices(500);
    solver.set_var3b_timept(1);
    solver.set_trace(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(-1);

    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout <<  "nb sol " << l_solutions.size() << endl;
    if (l_solutions.size()==1) { cout << " volume " << l_solutions.front().volume() << endl;
    }
      
    return 0;
}

