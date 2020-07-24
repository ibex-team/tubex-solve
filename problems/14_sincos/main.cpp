#include "tubex.h"
#include "tubex-solve.h"
#include <iomanip>

using namespace std;
 using namespace ibex;
using namespace tubex;

void contract(TubeVector& x)
{
   TFunction f("x1", "x2" ,"(-10*(x1-sin(x2))+cos(x2);1)");

  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 50000.)
    ctc_picard.contract(f, x, TimePropag::FORWARD);

  CtcDeriv ctc_deriv;
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.contract(x, f.eval_vector(x), TimePropag::FORWARD | TimePropag::BACKWARD);
}

int main()
{
  TFunction f("x1", "x2" ,"(10*(x1-sin(x2))+cos(x2);1)");
  /* =========== PARAMETERS =========== */
    Tube::enable_syntheses(false);
    IntervalVector v(2);
    vector<IntervalVector*> gates;    
    
    v[0]=Interval(0.0,0.0);
    v[1]=Interval(0.0,0.0);

    double volume=0.0;
    double totaltime=0.0;

    double step=1.;
    int nbsteps=1;
    for (int i=0; i< nbsteps; i++){

      Vector epsilon(2, 0.);

      double t0=i*step;
      if (t0>5) break;
      double t1=(i+1)*step;
      if (t1>5) t1=5;
      Interval domain(t0,t1);
      //      cout << " domain " << domain << endl;
      TubeVector x(domain, step, 2);
      x.set(v,t0 ); // initial condition
      //      v[0]=Interval();
      //      v[1]=Interval(t1,t1);
      //      x.set(v,t1 ); 
      /*
      for (int i=1; i<10; i++){
      v[0]=Interval();
      v[1]=Interval(i*t1/10.,i*t1/10.);
      x.set(v,i*t1/10);
      }
      */

  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);

    //solver.set_refining_fxpt_ratio(0.999);
    solver.set_refining_fxpt_ratio(2);

    //solver.set_propa_fxpt_ratio(0.9999);
    solver.set_propa_fxpt_ratio(0.);


    //solver.set_var3b_fxpt_ratio(0.);
    //    solver.set_var3b_fxpt_ratio(0.9999);
    solver.set_var3b_fxpt_ratio(-1);
    //solver.set_var3b_propa_fxpt_ratio(0.9999);
    solver.set_var3b_timept(1);
    solver.set_max_slices(400000);
    solver.set_refining_mode(0);
    solver.set_contraction_mode(2);
    solver.set_trace(1);
    solver.set_bisection_timept(-2);

    list<TubeVector> l_solutions = solver.solve(x, f);
    //    cout << "time " << (i+1)*step << " nb sol " << l_solutions.size() << endl;
    if (l_solutions.size()==1) {// cout << " volume " << l_solutions.front().volume() << endl;
      volume+=l_solutions.front().volume();
      totaltime+=solver.solving_time;
      v[0]=l_solutions.front()[0].last_slice()->output_gate();
      v[1]=l_solutions.front()[1].last_slice()->output_gate();
      IntervalVector* v1 = new IntervalVector(v);
      gates.push_back(v1);}
    else break;
    }
    /*
    for (int k=0; k< gates.size(); k++)
      cout << k << "  " << *(gates[k]) << endl;
    */
    cout << " last gate " << *(gates[gates.size()-1]) << endl;
    cout << " total volume " << volume << " total time " << totaltime << endl;
    return 0;
}
