/** 
 *  tubex-lib - Examples
 *  Solver testcase 17
 * ----------------------------------------------------------------------------
 *
 *  \date       2019
 *  \author     Bertrand Neveu
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


#include "math.h"
#include "tubex.h"
#include "tubex-solve.h"

using namespace std;
using namespace ibex;
using namespace tubex;

void contract(TubeVector& x)
{
  tubex::Function f("x1", "x2" ,"(x2;x2/0.2)");


  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 1.e300)
    ctc_picard.contract(f, x, FORWARD | BACKWARD);
   
  TubeVector v = f.eval_vector(x);

  CtcDeriv ctc_deriv;
  ctc_deriv.set_fast_mode(true);

  ctc_deriv.contract(x, v, FORWARD | BACKWARD);

 
  
  v=f.eval_vector(x);

  
  CtcDynCid* ctc_dyncid = new CtcDynCid(f);     
  //  CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f);     
  ctc_dyncid->set_fast_mode(true);
  CtcIntegration ctc_integration(f,ctc_dyncid);
  

  
  // cout << x << " volume before " << x.volume() << endl;;
 
  
  ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;
  ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;
  //cout << x << " volume after " << x.volume() << endl;;
  
  delete ctc_dyncid;
 

  
  
}

int main()
{
  tubex::Function f("x1", "x2" ,"(x2;x2/0.2)");
  /* =========== PARAMETERS =========== */

    Tube::enable_syntheses(false);

    Interval domain(0.,1.);
    TubeVector x(domain, 0.1, 2);
    IntervalVector v(2);
    v[0]=Interval(1.,1.);
    //    v[1]=Interval(-18.,18.);
    v[1]=Interval(-10.,10.);
    x.set(v, 0.); // ini
    v[0]=Interval(0.,0.);
    //v[1]=Interval(-18.,18.);
    v[1]=Interval(-10.,10.);
    x.set(v,1.);

    double eps=0.02;

  /* =========== SOLVER =========== */
    Vector epsilon(2, eps);


    tubex::Solver solver(epsilon);

    //    solver.set_refining_fxpt_ratio(0.9999);
    solver.set_refining_fxpt_ratio(2.0);
    solver.set_propa_fxpt_ratio(0.9999);
    
    //    solver.set_propa_fxpt_ratio(0.1);
    solver.set_var3b_fxpt_ratio(0.999);
    //    solver.set_var3b_fxpt_ratio(0.);

      
    solver.set_var3b_timept(0);
    solver.set_bisection_timept(3);
    solver.set_max_slices(20000);
    solver.set_refining_mode(0);
    solver.set_contraction_mode(2);
    solver.set_trace(1);
    list<TubeVector> l_solutions = solver.solve(x, f);
    cout << "nb sol " << l_solutions.size() << endl;
    return 0;
}
