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
  tubex::Function f("y1", "y2", "(-0.7*y1 ; 0.7*y1 - (ln(2)/5.)*y2)");

  
  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 50000.)
    ctc_picard.contract(f, x, FORWARD | BACKWARD);
  TubeVector v = f.eval_vector(x);
  /*
  CtcDeriv ctc_deriv;
  ctc_deriv.preserve_slicing(false);
  ctc_deriv.contract(x, v, FORWARD | BACKWARD);
  */
      
  CtcDynCid* ctc_dyncid = new CtcDynCid(f);     
  //CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f);     
  ctc_dyncid->set_fast_mode(true);
  CtcIntegration ctc_integration(f,ctc_dyncid);

  ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;

  ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;

  delete ctc_dyncid;
  
  
  // Check if the following is useful:   seems useless (BN)
  /*
  CtcEval ctc_eval;
  Interval t(1.,3.);
  ctc_eval.contract(Interval(1.,3.), Interval(1.1,1.3), x[1], v[1]);
  */
}

int main()
{
  /* =========== PARAMETERS =========== */
  tubex::Function f("y1", "y2", "(-0.7*y1 ; 0.7*y1 - (ln(2)/5.)*y2)");
    Tube::enable_syntheses(false);
    int n = 2;
    Interval domain(0.,6.);
    TubeVector x(domain, 0.003, IntervalVector(n, Interval(-9999.,9999.))); // todo: remove bounds
    //    TubeVector x(domain, 0.1, IntervalVector(n, Interval(-9999.,9999.))); // todo: remove bounds
    
    //Vector epsilon(n); epsilon[0] = 0.15; epsilon[1] = 0.15;
    //    Vector epsilon(n); epsilon[0] = 0.03; epsilon[1] = 0.03;
    Vector epsilon(n); epsilon[0] = 0.04 ; epsilon[1] = 0.04;

    // Boundary condition:
    IntervalVector init = x(x.domain().lb());
    init[0] = 1.25;
    x.set(init, x.domain().lb());

    // Additional restriction (maximum value):
    Interval domain_restriction(1.,3.);
    IntervalVector max_restriction(2);
    max_restriction[1] = Interval(1.1,1.3);
    x.set(max_restriction, domain_restriction);

  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);
    //    solver.set_refining_fxpt_ratio(0.99999);
    //    solver.set_refining_fxpt_ratio(0.999);
    solver.set_refining_fxpt_ratio(2.0);
    solver.set_propa_fxpt_ratio(0.99999);

    solver.set_var3b_fxpt_ratio(0.1);
    //    solver.set_var3b_fxpt_ratio(0.);
    

    solver.set_var3b_timept(2);
 
    solver.set_trace(1);
    solver.set_max_slices(40000);
    solver.set_refining_mode(3);
    solver.set_contraction_mode(2);
    solver.set_bisection_timept(-2);
    // Displaying the additional restriction:
    //    solver.figure()->draw_box(domain_restriction, max_restriction, "blue");

    list<TubeVector> l_solutions = solver.solve(x, f);

    
  return EXIT_SUCCESS;
}
