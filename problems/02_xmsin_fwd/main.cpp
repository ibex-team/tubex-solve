/** 
 *  tubex-lib - Examples
 *  Solver testcase
 * ----------------------------------------------------------------------------
 *
 *  \date       2018
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "tubex.h"
#include "tubex-solve.h"
#include <iomanip>

using namespace std;
using namespace ibex;
using namespace tubex;

/*
void contract(TubeVector& x)
{
  tubex::Function f("x", "-sin(x)");
  ibex::Function f1("x", "-sin(x)");
  CtcPicard ctc_picard;
  

  ctc_picard.preserve_slicing(false);
  

  if (x.volume() > 5000)
    ctc_picard.contract(f, x, FORWARD );

  /*
  TubeVector v = f.eval_vector(x);
  CtcDynCid* ctc_dyncid = new CtcDynCid(f1);     
  ctc_dyncid->set_fast_mode(true);
  CtcIntegration ctc_integration(f1,ctc_dyncid);
  
  ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;
  
  ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;
  
  delete ctc_dyncid;

  

  

    
  CtcDeriv ctc_deriv;
  ctc_deriv.preserve_slicing(false);
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.contract(x, f.eval_vector(x), FORWARD |  BACKWARD );
  
}
*/
int main()
{
  /* =========== PARAMETERS =========== */
  tubex::Function f("x", "-sin(x)");
    Tube::enable_syntheses(false);
    Vector epsilon(1,0.002); 
    double tf=10.;
    IntervalVector v(1);
    
    vector<IntervalVector*> gates; 
    Interval domain(0.,10.);
    double volume=0.0;
    double totaltime=0.0;
    double step=10.;
    int nbsteps=1;
    TrajectoryVector truth(domain, tubex::Function("2.*atan(exp(-t)*tan(0.5))"));
    v[0]= truth(Interval(0.))[0];
    
    for (int i=0; i< nbsteps; i++){
      double t0=i*step;
      double t1=(i+1)*step;
      Interval domain(t0,t1);
      TubeVector x(domain,step, 1);
      x.set(v,t0 ); // initial condition
    // Note: use truth(Interval(0.)) instead of truth(0.) for a reliable evaluation
 
  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);
    //    solver.set_refining_fxpt_ratio(0.99999);
    solver.set_refining_fxpt_ratio(2.0);
    //    solver.set_propa_fxpt_ratio(0.99);
    solver.set_var3b_fxpt_ratio(0.);
    //solver.set_var3b_fxpt_ratio(0.99);
    solver.set_var3b_propa_fxpt_ratio(0.99);
    solver.set_var3b_fxpt_ratio(-1);

    solver.set_var3b_timept(1);
    solver.set_max_slices(40000);
    solver.set_refining_mode(0);
    solver.set_contraction_mode(2);
    solver.set_trace(1);
    //    solver.figure()->add_trajectoryvector(&truth, "truth");
    list<TubeVector> l_solutions = solver.solve(x, f);
    if (l_solutions.size()==1) { cout << " volume " << l_solutions.front().volume() << endl;
       volume+=l_solutions.front().volume();
       totaltime+=solver.solving_time;
       v[0]=l_solutions.front()[0].last_slice()->output_gate();
       IntervalVector* v1 = new IntervalVector(v);
       gates.push_back(v1);}
     else break;
    }
    cout << " total volume " << volume << " total time " << totaltime << endl;
    for (int k=0; k< gates.size(); k++)
      cout << k << "  " << *(gates[k]) << endl;
    cout << " truth " << truth(Interval(10.0)) << endl;
  // Checking if this example still works:
  //return (solver.solutions_contain(l_solutions, truth) != NO) ? EXIT_SUCCESS : EXIT_FAILURE;
    return 0;
}
