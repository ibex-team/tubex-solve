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

class FncDelayCustom : public tubex::Fnc
{
  public: 

    FncDelayCustom(double delay) : Fnc(1, 1, true), m_delay(delay) { }
    const Interval eval(int slice_id, const TubeVector& x) const { cout << "not defined 1" << endl; }
    const Interval eval(const Interval& t, const TubeVector& x) const { cout << "not defined 2" << endl; }
    const Interval eval(const IntervalVector& x) const { cout << "not defined 3" << endl; }
    const IntervalVector eval_vector(const IntervalVector& x) const { cout << "not defined 4" << endl; }

    const IntervalVector eval_vector(int slice_id, const TubeVector& x) const
    {
      Interval t = x[0].slice(slice_id)->domain();
      return eval_vector(t, x);
    }

    const IntervalVector eval_vector(const Interval& t, const TubeVector& x) const
    {
      IntervalVector eval_result(x.size(), Interval::EMPTY_SET);

      if((t - m_delay).lb() <= x.domain().lb())
        eval_result |= x(t);

      if((t - m_delay).ub() >= x.domain().lb())
        eval_result |= exp(m_delay) * x((t - m_delay) & x.domain());

      return eval_result;
    }

  protected:

    double m_delay = 0.;
};

void contract(TubeVector& x)
{
  double delay = 0.5;
  FncDelayCustom f(delay);

  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 50000.)
    ctc_picard.contract(f, x, FORWARD | BACKWARD);

  // Computing the derivative v from a delay of x
  CtcDelay ctc_delay;
  TubeVector v(x, IntervalVector(x.size()));
  ctc_delay.contract(delay, x, v);
  v *= exp(delay);

  CtcDeriv ctc_deriv;
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.preserve_slicing(false);
  ctc_deriv.contract(x, v, FORWARD | BACKWARD);
}

int main()
{
  /* =========== PARAMETERS =========== */

    Tube::enable_syntheses(false);
    int n = 1;
    Vector epsilon(n, 100.);
    Interval domain(0.,5.);
    TubeVector x(domain, n);
    TrajectoryVector truth(domain, tubex::Function("exp(t)"));
    //delete?
    double t_value = domain.lb();
    Interval init_value(exp(t_value));
    IntervalVector v(1);
    vector<IntervalVector*> gates;
    v[0]= init_value;
    double volume=0.0;
    double totaltime=0.0;
    //    double step=0.0005;
    //    int nbsteps=10000;
    double step=5;
    int nbsteps=1;
  /* =========== SOLVER =========== */
    for (int i=0; i< nbsteps; i++){
      tubex::Solver solver(epsilon);
      solver.set_trace(0);
      double t0=i*step;
      double t1=(i+1)*step;
      Interval domain(t0,t1);
      TubeVector x(domain, step/1000, 1);
      x.set(v[0],t0);
      solver.set_refining_fxpt_ratio(0.99999);
      solver.set_propa_fxpt_ratio(1.);
      solver.set_cid_fxpt_ratio(0.);
      solver.set_cid_propa_fxpt_ratio(0.99);
      solver.set_max_slices(10000);
      solver.set_refining_mode(0);
      solver.set_trace(1);
      //      solver.set_trace(1);
      //    solver.figure()->add_trajectoryvector(&truth, "truth");
      list<TubeVector> l_solutions = solver.solve(x, &contract);
      //  cout << "time " << (i+1)*step <<  "nb sol " << l_solutions.size() << endl;
      if (l_solutions.size()==1) { 
	// cout << " volume " << l_solutions.front().volume() << endl;
	volume+=l_solutions.front().volume();
	totaltime+=solver.solving_time;
	v[0]=l_solutions.front()[0].last_slice()->output_gate();
	IntervalVector* v1 = new IntervalVector(v);
	gates.push_back(v1);}
      else break;
    }
    /*
    for (int k=0; k < gates.size() ; k++)
      cout << k << "  " << *(gates[k]) << endl;
    */
    cout << " last gate " << *gates[gates.size()-1] << endl;
    cout << " total volume " << volume << " total time " << totaltime << endl;
    return 0;

  // Checking if this example still works:
    //  return (solver.solutions_contain(l_solutions, truth) == YES) ? EXIT_SUCCESS : EXIT_FAILURE;
}
