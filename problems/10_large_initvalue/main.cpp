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
  tubex::Function f("x", "-x");
  ibex::Function f1("x", "-x");

  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 50000.0)
    ctc_picard.contract(f, x);


  TubeVector v = f.eval_vector(x);
  /*
   CtcCidSlicing ctc_cidslicing (f1);
   

   ctc_cidslicing.preserve_slicing(false);
   ctc_cidslicing.contract(x,v,FORWARD,false);
  */
   //   ctc_cidslicing.contract(x,v,BACKWARD,false);



  CtcDeriv ctc_deriv;
  ctc_deriv.preserve_slicing(false);
  ctc_deriv.contract(x, f.eval_vector(x));

  
}

int main()
{
  /* =========== PARAMETERS =========== */

    Tube::enable_syntheses(false);
    int n = 1;
    //    Vector epsilon(n, 0.1);
    Vector epsilon(n, 10.);
    Interval domain(0.,1.);
    TubeVector x(domain,0.1, n);
    x.set(IntervalVector(n, Interval(0.5,1.)*exp(Interval(-0.))), 0.); // initial condition
    TrajectoryVector truth1(domain, tubex::Function("1.0*exp(-t)"));
    TrajectoryVector truth2(domain, tubex::Function("0.5*exp(-t)"));

  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);
    solver.set_refining_fxpt_ratio(0.9999);
    solver.set_propa_fxpt_ratio(0.9999);
    //solver.set_cid_fxpt_ratio(0.9);
    solver.set_cid_fxpt_ratio(0.);

    solver.set_cid_propa_fxpt_ratio(0.999);
    solver.set_cid_timept(1);
    solver.set_trace(1);
    solver.set_max_slices(20000);
    solver.set_refining_mode(3);
    //    solver.figure()->add_trajectoryvector(&truth1, "truth1");
    //    solver.figure()->add_trajectoryvector(&truth2, "truth2");
    list<TubeVector> l_solutions = solver.solve(x, &contract);


  // Checking if this example still works:
  Tube hull = TubeVector::hull(l_solutions)[0];
  Tube f_hull = Tube(domain, 0.0001, tubex::Function("[0.5,1.0]*exp(-t)"));
  return f_hull.is_subset(hull) ? EXIT_SUCCESS : EXIT_FAILURE;
}
