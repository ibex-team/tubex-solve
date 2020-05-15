/** 
 *  tubex-lib - Examples
 *  Solver testcase 20
 *  \date       2020
 *  \author     Bertrand Neveu
 */

#include "tubex.h"
#include "tubex-solve.h"
#include <iomanip>

using namespace std;
using namespace ibex;
using namespace tubex;

void contract(TubeVector& x)
{
  tubex::Function f("x1", "x2" ,"(-x2+0.1*x1*(1-x1^2-x2^2);x1+0.1*x1*(1-x1^2-x2^2))");

  CtcPicard ctc_picard;
  ctc_picard.preserve_slicing(false);
  if (x.volume() > 1.e100) 
    ctc_picard.contract(f, x, FORWARD);

   if (x.volume() < 1.e100){
     TubeVector v = f.eval_vector(x);

     CtcDeriv ctc_deriv;
     ctc_deriv.set_fast_mode(true);
     ctc_deriv.contract(x, v, FORWARD | BACKWARD);
     v=f.eval_vector(x);

  
     CtcDynCid* ctc_dyncid = new CtcDynCid(f);     
     ctc_dyncid->set_fast_mode(true);
     CtcIntegration ctc_integration(f,ctc_dyncid);
     ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;
     ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;
     delete ctc_dyncid;


   }   

}

int main()
{
  tubex::Function f("x1", "x2" ,"(-x2+0.1*x1*(1-x1^2-x2^2);x1+0.1*x1*(1-x1^2-x2^2))");
  /* =========== PARAMETERS =========== */
    Tube::enable_syntheses(false);
    IntervalVector v(2);
    vector<IntervalVector*> gates;    
    
    v[0]=Interval(1.0,1.0).inflate(0.1);
    v[1]=Interval(0.0,0.0).inflate(0.1);

    double volume=0.0;
    double totaltime=0.0;

    double step=4.;  
    int nbsteps=1;
    for (int i=0; i< nbsteps; i++){

      Vector epsilon(2, 1.);

      double t0=i*step;
      if (t0>5) break;
      double t1=(i+1)*step;
      if (t1>5) t1=5;
      Interval domain(t0,t1);
      //      cout << " domain " << domain << endl;
      TubeVector x(domain, step, 2);
      x.set(v,t0 ); // initial condition
      
      /* =========== SOLVER =========== */

      tubex::Solver solver(epsilon);

      //    solver.set_refining_fxpt_ratio(0.99999);
      solver.set_refining_fxpt_ratio(2.0);

      //      solver.set_propa_fxpt_ratio(0.9999);
      solver.set_propa_fxpt_ratio(0.99);
      //      solver.set_var3b_propa_fxpt_ratio(0.999);

      solver.set_var3b_fxpt_ratio(-1);
      //solver.set_var3b_fxpt_ratio(0.9);
      //      solver.set_var3b_timept(1);
      solver.set_max_slices(2000);
      solver.set_refining_mode(0);
      solver.set_trace(1);
      solver.set_bisection_timept(-1);
      solver.set_contraction_mode(2);
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
