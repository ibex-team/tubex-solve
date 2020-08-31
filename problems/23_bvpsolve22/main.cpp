/** 
 *  tubex-lib - Examples
 *  Solver testcase  23  bvpsolve22 (case ksi=0.1)
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


using namespace std;
using namespace ibex;
using namespace tubex;



int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
     TFunction f("x1", "x2" ,"(x2;-10*(x2+x1^2))");
     //tubex::Function f("x1", "x2" ,"(x2;-100*(x2+x1^2))");
    Interval domain(0.,1.);
    TubeVector x(domain,2);
    IntervalVector v(2);
    v[0]=Interval(0.,0.);
    v[1]=Interval(-20,20);
    x.set(v, 0.); // ini
    v[0]=Interval(0.5,0.5);
    v[1]=Interval(-20,20);

    x.set(v,1.);

    double eps0=0.05;
    double eps1=0.05;

    /* =========== SOLVER =========== */
    Vector epsilon(2);
	epsilon[0]=eps0;
	epsilon[1]=eps1;


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.);
    solver.set_propa_fxpt_ratio(0.999);
    //solver.set_propa_fxpt_ratio(0.);

    solver.set_var3b_fxpt_ratio(0.999);
    //solver.set_var3b_fxpt_ratio(-1);

    solver.set_var3b_external_contraction(false);
    solver.set_var3b_propa_fxpt_ratio(0.999);
    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(5000);
    //solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(2);
    solver.set_stopping_mode(0);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x,f, &contract);
    list<TubeVector> l_solutions = solver.solve(x,f);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;
    
}

