/** 
 *  tubex-lib - Examples
 *  Solver testcase 15
 * ----------------------------------------------------------------------------
 *
 *  \date       2019
 *  \author     Bertrand Neveu
 *  \copyright  Copyright 2019 Bertrand Neveu
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
  
  // Differential equation

  tubex::Function f("x1", "x2", "(x2;-x1)");
  ibex::Function f1("x1", "x2", "(x2;-x1)");

    CtcPicard ctc_picard;
    ctc_picard.preserve_slicing(true);
    if (x.volume() > 1.e100)
      ctc_picard.contract(f, x);

    if (x.volume() < 1.e100){
    TubeVector v = f.eval_vector(x);
    
    
      CtcDeriv ctc_deriv;
      ctc_deriv.set_fast_mode(true);
      ctc_deriv.contract(x, v);
      v=f.eval_vector(x);
      
      //      CtcDynCid* ctc_dyncid = new CtcDynCid(f1);     
      CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f1);     
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
  double pi=M_PI;
    Tube::enable_syntheses(false);
    int n = 2;

    Vector epsilon(n, 0.0005);
    Interval domain(0.,pi/2);
    TubeVector x(domain, n);
    IntervalVector x0(2);
    IntervalVector x1(2);
    x0[0]=Interval(0.,0.);
    x0[1]=Interval(-1.e8,1.e8);
    x.set(x0,0.);
    x1[0]=Interval(2.,2.);
    x1[1]=Interval(-1.e8,1.e8);
    x.set(x1,pi/2);
    cout << " x " << x << endl;
    
  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);
    //    solver.set_refining_fxpt_ratio(0.9999);
    //    solver.set_refining_fxpt_ratio(0.9995);
    solver.set_refining_fxpt_ratio(2.0);
    solver.set_propa_fxpt_ratio(0.9999);
    solver.set_var3b_fxpt_ratio(0.9999);
    //    solver.set_var3b_fxpt_ratio(0.);
    solver.set_var3b_propa_fxpt_ratio(0.9999);
    solver.set_trace(1);
    solver.set_var3b_timept(0);
    solver.set_bisection_timept(3);
    solver.set_max_slices(20000);

    solver.set_refining_mode(3);
    list<TubeVector> l_solutions = solver.solve(x, &contract);

    return (0);
}
