/** 
 *  tubex-lib - Examples
 *  Solver testcas 04
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
#include "ibex_CtcHC4.h"
#include "ibex_SystemFactory.h"

using namespace std;
using namespace ibex;
using namespace tubex;

void contract(TubeVector& x)
{
  // Boundary constraints

    Variable vx0, vx1;
    SystemFactory fac;
    fac.add_var(vx0);
    fac.add_var(vx1);
    fac.add_ctr(sqr(vx0) + sqr(vx1) = 1);
    System sys(fac);
    ibex::CtcHC4 hc4(sys);
    IntervalVector bounds(2);
    bounds[0] = x[0](0.);
    bounds[1] = x[0](1.);
    //    cout << " bounds " <<  bounds[0] << " " << bounds[1] << endl;
    hc4.contract(bounds);
    x.set(IntervalVector(bounds[0]), 0.);
    x.set(IntervalVector(bounds[1]), 1. );
  

  // Differential equation

    tubex::Function f("x", "x");

    //    cout << " x before Picard " << x << x.volume() << endl;
    //    cout << " first slice " << *(x[0].first_slice()) << endl;
    CtcPicard ctc_picard;
    
    /*    ctc_picard.preserve_slicing(true);
    if (x.volume() > 1.e100)
      ctc_picard.contract(f, x);
    */
    //    cout << " x after Picard " << x << x.volume() << endl;
    //    cout << " first slice " << *(x[0].first_slice()) << endl;

    /*
    CtcDeriv ctc_deriv;
    ctc_deriv.preserve_slicing(false);
    ctc_deriv.set_fast_mode(true);
    ctc_deriv.contract(x, f.eval_vector(x));
    */
    
    /*
    TubeVector v = f.eval_vector(x);
    CtcDynCid* ctc_dyncid = new CtcDynCid(f);     
    //CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f);     
    ctc_dyncid->set_fast_mode(true);
    CtcIntegration ctc_integration(f,ctc_dyncid);

    ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;

    ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;

    delete ctc_dyncid;
    */
     
}

int main()
{
    tubex::Function f("x", "x");
  /* =========== PARAMETERS =========== */

    Tube::enable_syntheses(false);
    int n = 1;
    
    Vector epsilon(n, 0.0005);
    Interval domain(0.,1.);
    //    TubeVector x(domain, n, Interval (-1.e100,1.e100));
    TubeVector x(domain, n);
    TrajectoryVector truth1(domain, tubex::Function("exp(t)/sqrt(1+exp(2))"));
    TrajectoryVector truth2(domain, tubex::Function("-exp(t)/sqrt(1+exp(2))"));

  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);
    //    solver.set_refining_fxpt_ratio(0.9995);
    //    solver.set_refining_fxpt_ratio(0.9999);
    solver.set_refining_fxpt_ratio(2.0);
    //    solver.set_refining_fxpt_ratio(0.9999);
    //    solver.set_propa_fxpt_ratio(0.999);
    solver.set_propa_fxpt_ratio(0.);
    solver.set_var3b_propa_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(0.999);
    // solver.set_var3b_fxpt_ratio(0.);

    solver.set_trace(1);
    solver.set_var3b_timept(2);
    solver.set_bisection_timept(3);
    solver.set_max_slices(10000);
    solver.set_refining_mode(0);
    solver.set_contraction_mode(2);
    //    solver.figure()->add_trajectoryvector(&truth1, "truth1");
    //    solver.figure()->add_trajectoryvector(&truth2, "truth2");
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);


  // Checking if this example still works:
  return (solver.solutions_contain(l_solutions, truth1) == YES
       && solver.solutions_contain(l_solutions, truth2) == YES) ? EXIT_SUCCESS : EXIT_FAILURE;
}
