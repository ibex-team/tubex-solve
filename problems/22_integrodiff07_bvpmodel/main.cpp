/** 
 *  tubex-lib - Examples
 *  Solver testcase
 * ----------------------------------------------------------------------------
 *
 *  \date       2020
 *  \author     Bertrand Neveu
 *  \copyright  Copyright 2020 Simon Rohou
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
    hc4.contract(bounds);
    IntervalVector x0(2);
    IntervalVector x1(2);
    x0[0]=bounds[0];
    x0[1]=0.0;
    x1[0]=bounds[1];
    x1[1]=x[1](1.);

    x.set(x0, 0.);
    x.set(x1, 1.);
}
/*
  // Differential equation
    tubex::Function f("x1", "x2", "(1-2*x1-5*x2;x1)");

    CtcPicard ctc_picard;
    ctc_picard.preserve_slicing(true);
    ctc_picard.contract(f, x, FORWARD | BACKWARD);

    CtcDeriv ctc_deriv;
    ctc_deriv.preserve_slicing(true);
    ctc_deriv.set_fast_mode(true);

    //    cout << f.eval_vector(x) << endl;
    ctc_deriv.contract(x, f.eval_vector(x), FORWARD | BACKWARD);
    //    cout << " after ctc deriv  " << x << endl;
}

*/


int main()
{
  /* =========== PARAMETERS =========== */
    TFunction f("x1", "x2", "(1-2*x1-5*x2;x1)");

    Tube::enable_syntheses(false);
    int n = 2;
    Vector epsilon(n, 0.02);
    Interval domain(0.,1.);
    TubeVector x(domain, n); 
    IntervalVector x0(2);
    IntervalVector x1(2);
    
    // x0[0]=Interval(-100,100);
    x0[1]=Interval(0,0);

    x.set(x0,0.);
    //    x1[0]=Interval(-100,100);
    // x1[1]=Interval(-10,10);
    //    x.set(x1,1.0);
    

    TrajectoryVector truth1(domain, TFunction("(exp(-t)*(-(cos(2*t)*(-1 + cos(4) + 2*sin(4) + 4*exp(1)*sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4))))) + sin(2*t)*(2 + 2*cos(4) - sin(4) + 2*exp(1)*(2*exp(1) + sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4)))))))/(5 + 3*cos(4) + 8*exp(2) - 4*sin(4))"));
    TrajectoryVector truth2(domain, TFunction("(exp(-t)*(-(cos(2*t)*(-1 + cos(4) + 2*sin(4) - 4*exp(1)*sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4))))) + sin(2*t)*(2 + 2*cos(4) + 4*exp(2) - sin(4) - 2*exp(1)*sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4))))))/(5 + 3*cos(4) + 8*exp(2) - 4*sin(4))"));

  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);
    //    solver.set_refining_fxpt_ratio(0.99999);
    solver.set_refining_fxpt_ratio(0.9);
    solver.set_propa_fxpt_ratio(0.);
    //solver.set_propa_fxpt_ratio(0.999);
    //    solver.set_var3b_fxpt_ratio(0.99);
    //    solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.999);
    solver.set_trace(1);
    solver.set_max_slices(400);
    solver.set_refining_mode(0);
    solver.set_var3b_timept(0);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(2);
    //    solver.figure()->add_trajectoryvector(&truth1, "truth1");
    //    solver.figure()->add_trajectoryvector(&truth2, "truth2");
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);

  // Checking if this example still works:
    return 0;
}

