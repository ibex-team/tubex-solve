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
  TFunction f("x", "-x");

  CtcPicard ctc_picard;
  
  ctc_picard.preserve_slicing(false);
   if (x.volume()> 500.0) 
     ctc_picard.contract(f, x, TimePropag::BACKWARD);
  
   //  cout << " after picard " << x << endl;
   //   cout << " volume " << x.volume()  << endl;
    
   
  CtcDeriv ctc_deriv;
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.preserve_slicing(false);
  ctc_deriv.contract(x, f.eval_vector(x),TimePropag::FORWARD | TimePropag::BACKWARD);

  /*
   TubeVector v = f.eval_vector(x);
   CtcDynCid* ctc_dyncid = new CtcDynCid(f);     
   ctc_dyncid->set_fast_mode(true);
   CtcIntegration ctc_integration(f,ctc_dyncid);
  
   ctc_integration.contract(x,v,x[0].domain().lb(),TimePropag::FORWARD) ;
  
   ctc_integration.contract(x,v,x[0].domain().ub(),TimePropag::BACKWARD) ;
  
   delete ctc_dyncid;
   */

}

int main()
{
  /* =========== PARAMETERS =========== */
  cout << " avant appel tubex function " << endl;

  TFunction f("x", "-x");
cout << " apres appel tubex function " << endl;
    Tube::enable_syntheses(false);
    Vector epsilon(1, 0.005);
    double tf=10.;
    IntervalVector v(1);
    vector<IntervalVector*> gates; 
    Interval domain0(0.,tf);
    TrajectoryVector truth(domain0, TFunction("exp(-t)"));
    //    x.set(IntervalVector(truth(Interval(tf))), tf); // final condition
    v[0]=Interval(exp(Interval(-tf)));  
    /* =========== SOLVER =========== */
    double volume=0.0;
    double totaltime=0.0;
    double step=10.;
    int nbsteps=1;


    for (int i=0; i< nbsteps; i++){

      double t0=tf-(i+1)*step;
      double t1=tf-i*step;
      Interval domain(t0,t1);
      //      TubeVector x(domain, step/10000, 1);
      TubeVector x(domain,  1);
      //      TubeVector x(domain, 10., 1);
      x.set(v,t1 ); // final condition at t1
      cout << " avant solver " << endl;
      tubex::Solver solver(epsilon);
      cout << " apres solver " << endl;
      //      solver.set_refining_fxpt_ratio(0.9999);
      solver.set_refining_fxpt_ratio(2.0);
      //      solver.set_propa_fxpt_ratio(0.9999);
      //solver.set_propa_fxpt_ratio(0.99);
      solver.set_propa_fxpt_ratio(0.);

      solver.set_var3b_fxpt_ratio(-1);
      //      solver.set_var3b_fxpt_ratio(0.99);
      solver.set_var3b_propa_fxpt_ratio(0.99);
      //
      solver.set_var3b_timept(0);
      solver.set_trace(1);
      solver.set_max_slices(40000);
      solver.set_refining_mode(0);
      solver.set_bisection_timept(-2);
      solver.set_contraction_mode(0);
      //    solver.figure()->add_trajectoryvector(&truth, "truth");
      cout << " avant solver " << endl;
      list<TubeVector> l_solutions = solver.solve(x, f);
      cout << " apres solve " << endl;
      //      cout << "time " << (i+1)*step <<  "nb sol " << l_solutions.size() << endl;
      if (l_solutions.size()==1) { cout << " volume " << l_solutions.front().volume() << endl;
	volume+=l_solutions.front().volume();
	totaltime+=solver.solving_time;
	v[0]=l_solutions.front()[0].first_slice()->input_gate();
	IntervalVector* v1 = new IntervalVector(v);
	gates.push_back(v1);}
      else {cout << " error : more than one solution " << endl; break;}
    }
    for (int k=gates.size()-1; k >=0; k--)
      cout << k << "  " << *(gates[k]) << endl;
    cout << " total volume " << volume << " total solving time " << totaltime << endl;
    cout << " gate at t=0.0 " << *(gates[gates.size()-1]) << endl;
    for (int k=gates.size()-1; k >=0; k--)
      delete gates[k];

    return 0;

  // Checking if this example still works:
  //return (l_solutions.size() == 1
    //       && solver.solutions_contain(l_solutions, truth) == YES) ? EXIT_SUCCESS : EXIT_FAILURE;
}
