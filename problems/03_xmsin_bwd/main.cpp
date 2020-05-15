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

using namespace std;
using namespace ibex;
using namespace tubex;


void contract(TubeVector& x)
{
  tubex::Function f("x", "-sin(x)");
  /*
  CtcPicard ctc_picard;
  
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 5000.0)
    ctc_picard.contract(f, x,  BACKWARD);
  
  
  CtcDeriv ctc_deriv;
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.preserve_slicing(false);
  ctc_deriv.contract(x, f.eval_vector(x),FORWARD | BACKWARD);
  */

  
  TubeVector v = f.eval_vector(x);

  CtcDynCid* ctc_dyncid = new CtcDynCid(f);     
  //CtcDynBasic* ctc_dyncid = new CtcDynBasic(f);     
  //CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f);     
  ctc_dyncid->set_fast_mode(true);
  CtcIntegration ctc_integration(f,ctc_dyncid);
  //  cout << " x before call to contract " << x << endl;
  //

  //  ctc_integration.set_picard_mode(true);
  
  ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;
  //  cout << " x after call to fwd contract " << x << endl;
  ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;

  //  cout << " x after end backward " << x << endl;
  delete ctc_dyncid;

}

int main()
{
  /* =========== PARAMETERS =========== */
  tubex::Function f("x", "-sin(x)");
  Tube::enable_syntheses(false);
  Vector epsilon(1,0.005);
  double tf=10.;

  vector<IntervalVector*> gates; 
  Interval domain0(0.,tf);
    TrajectoryVector truth(domain0, tubex::Function("2.*atan(exp(-t)*tan(0.5))"));
     //    x.set(IntervalVector(truth(Interval(tf))), tf); // final condition

    //    v[0]=Interval(exp(Interval(-tf)));  /* =========== SOLVER =========== */
    IntervalVector v(truth(Interval(tf)));
    cout << " v " << v << endl;
    double volume=0.0;
    double totaltime=0.0;
    double step=10.;
    int nbsteps=1;
    TubeVector y(domain0,0.05,1);
    for (int i=0; i< nbsteps; i++){

      double t0=tf-(i+1)*step;
      double t1=tf-i*step;
      Interval domain(t0,t1);
      TubeVector x(domain, step, 1);
      // TubeVector x(domain, step, 1);
      x.set(v,t1 ); // final condition at t1

      tubex::Solver solver(epsilon);
      //      solver.set_refining_fxpt_ratio(0.999);
      solver.set_refining_fxpt_ratio(2.0);
      //      solver.set_propa_fxpt_ratio(0.99);
      solver.set_propa_fxpt_ratio(0);
      solver.set_var3b_propa_fxpt_ratio(0.99);
      //solver.set_var3b_fxpt_ratio(0.99);
      solver.set_var3b_fxpt_ratio(-1.);
      solver.set_bisection_timept(-2);
      solver.set_var3b_timept(0);
      solver.set_refining_mode(0);
      solver.set_contraction_mode(2);
      solver.set_max_slices(40000);
      solver.set_trace(1);
      //    solver.figure()->add_trajectoryvector(&truth, "truth");
      list<TubeVector> l_solutions = solver.solve(x, f);
      //      list<TubeVector> l_solutions = solver.solve(x, &contract);
      cout << "time " << (i+1)*step << endl ;
      cout << l_solutions.size() << endl;
      if (l_solutions.size()==1) 
	{ cout << " volume " << l_solutions.front().volume() << endl;
	volume+=l_solutions.front().volume();
	totaltime+=solver.solving_time;
	v[0]=l_solutions.front()[0].first_slice()->input_gate();
	IntervalVector* v1 = new IntervalVector(v);
	gates.push_back(v1);}
      else {cout << " error : more than one solution " << endl; break;}
    }

    for (int k=0; k< gates.size() ; k++)
      cout << k << "  " << *(gates[k]) << endl;
    cout << " total volume " << volume << " total solving time " << totaltime << endl;
    for (int k=gates.size()-1; k >=0; k--)
      delete gates[k];
    cout << " truth(0) " << truth(Interval(0.0)) << endl;
    cout << " truth(10) " << truth(Interval(10.0)) << endl;

    return 0;

  // Checking if this example sill works:
  //return (l_solutions.size() == 1
    //       && solver.solutions_contain(l_solutions, truth) == YES) ? EXIT_SUCCESS : EXIT_FAILURE;
}
