/* ============================================================================
 *  tubex-lib - Solver class
 * ============================================================================
 *  Copyright : Copyright 2017 Simon Rohou
 *  License   : This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 *
 *  Author(s) : Simon Rohou, Bertrand Neveu
 *  Bug fixes : -
 *  Created   : 2018
 * ---------------------------------------------------------------------------- */

#ifndef __TUBEX_SOLVER_H__
#define __TUBEX_SOLVER_H__

#include <list>
#include <vector>

#include "tubex_TubeVector.h"
#include "tubex_TrajectoryVector.h"
#include "tubex_VIBesFigTubeVector.h"
#include "tubex_CtcPicard.h"

#include "tubex_CtcDeriv.h"
#include "tubex_CtcIntegration.h"
#include "tubex_CtcDynCid.h"
#include "tubex_CtcDynCidGuess.h"
#include "tubex_CtcDynBasic.h"

using namespace std;
namespace tubex
{
  class Solver
  {
  public:

      Solver(const ibex::Vector& max_thickness);
      ~Solver();
      /* Ratios used for stopping fixed point algorithms : all ratios are about the tube volume.
	 For all these ratios, their possible values are :
	 -1 : no call of the fixed point algorithm
	 0  : one only call of the fixed point algorithm
	 a value x between 0 and 1: the algorithm is called while the volume is reduced by a factor more than 1-x : for example, x=0.9 : the algorithm is called while the volume is reduced by more than 10%.
	 a value > 1 : only used for refining_fxpt_ratio to ensure that the slicing will be done until max_slices is reached.
      */

      /* slicing fixed point ratio : used for stopping the slicing*/
      void set_refining_fxpt_ratio(float refining_fxpt_ratio); 

      /* contraction fixed point ratio : used by the algorithm managing the whole sequence of contractors (external contract fonction and ODE contractor) called after a slicing and after a bisection*/
      void set_propa_fxpt_ratio(float propa_fxpt_ratio);
      /* var3b fix point ratio (for managing the sequence of var3b calls) */
      void set_var3b_fxpt_ratio(float var3b_fxpt_ratio);
      /* for managing the sequence of contractors used by the var3b subcontraction. (ctcDeriv and possibly CtcVnode).*/
      void set_var3b_propa_fxpt_ratio(float propa_fxpt_ratio);

      /* Slice choice for var3b
      where do the var3b contraction : 
      1 domain.ub(); 
      -1 domain.lb(); 
      0 max_diam_gate(); 
      2 randomly domain.ub() or domain.lb() */
      void set_var3b_timept(int var3b_timept); 

      /* Time choice for bisection 
      1 domain.ub(); 
      -1 domain.lb(); 
      0 max_diam_gate(); 
      2 randomly domain.ub() or domain.lb(); 
      3 round robin  domain.ub() and domain.lb() ; 
      -2 no bisection */
      void set_bisection_timept (int bisection_timept); 

      /* slicing limit :  no more refining when  max_slices is  reached */
      void set_max_slices(int max_slices); 

      /* refining mode : which slices to refine
         0 : all slices ; 
         1 : one slice (the steepest) ; 
         2 : the slices with a difference between input and output gates greater  than average ; 
         3 : the slices with a difference between input and output gates greater  than the median.
      */
      void set_refining_mode(int refining_mode); 
 
      /* contraction mode : the ODE contractor called   
      0 for CtcDynBasic ;
      1 for CtcDynCid ;
      2 for CtcDynCidGuess ;
      4 for CtcDeriv ;*/
      void set_contraction_mode(int contraction_mode) ;

      /* stopping mode : the stopping criterion in one branch of the search tree */
      void set_stopping_mode(int stopping_mode) ; // 0 for max tube diam , 1 for max gate diam

      /* calling external contraction (function contract) during var3b subcontractions 
      0 for not calling external contraction during var3b; 
      1 for calling external contractor during var3b used for calling  ctcvnode inside var3b contractions.*/
      void set_var3b_external_contraction (bool external_contraction) ; 

      /* trace during search : 
       0 for no trace
       1 for trace message at each refining, bisection and solution 
      */
      void set_trace(int trace);  

     
      /* the solve method, it has for parameters a tube vector x0 , and 3 possibilities
         -  a tube vector contractor ctc_func (for general problems as Integrodifferential problems and/or for using ctcvnode )
         - a differential function (TFnc computing the derivative of a tube vector) (for pure ODEs)
         - a differential function and a tube vector contractor (for ODE problems with side constraints and/or calls to ctcVnode)
         The returned results are tubes containing the solutions. Two tubes have at least a disjoint gate.
      */

      const std::list<TubeVector> solve(const TubeVector& x0, TFnc & f,void (*ctc_func)(TubeVector&, double t0, bool incremental)=NULL );
      const std::list<TubeVector> solve(const TubeVector& x0, TFnc* f,void (*ctc_func)(TubeVector&, double t0, bool incremental));
      const std::list<TubeVector> solve(const TubeVector& x0, void (*ctc_func)(TubeVector&, double t0, bool incremental));

      VIBesFigTubeVector* figure();
      static const ibex::BoolInterval solutions_contain(const std::list<TubeVector>& l_solutions, const TrajectoryVector& truth);
      /* the solving time of a solve call */
      double solving_time;

  protected:
      double one_finite_gate(const TubeVector &x);
      bool empty_intersection(TubeVector& t1, TubeVector& t2);
      void clustering(std::list<std::pair<int,TubeVector> >& l_tubes);
      void clustering(std::list<TubeVector>& l_tubes);
      bool stopping_condition_met(const TubeVector& x);
      bool gate_stopping_condition(const TubeVector& x);
      bool diam_stopping_condition(const TubeVector& x);
      bool fixed_point_reached(double volume_before, double volume_after, float fxpt_ratio);

      void bisection (const TubeVector &x, list<pair<pair<int,double>,TubeVector> > &s, int level);
      void fixed_point_contraction (TubeVector &x, TFnc* f, void (*ctc_func)(TubeVector&, double t0, bool incremental), float propa_fxpt_ratio, bool incremental, double t0, bool v3b=false);
      void contraction (TubeVector &x, TFnc * f,
			void (*ctc_func) (TubeVector&,double t0,bool incremental),
			bool incremental, double t0 , bool v3b);
      void deriv_contraction (TubeVector &x, const TFnc& f);
      void integration_contraction(TubeVector &x, const TFnc& f, double t0, bool incremental);
      void picard_contraction (TubeVector &x, const TFnc& f);
      void var3b(TubeVector &x,TFnc* f, void (*ctc_func)(TubeVector&, double t0, bool incremental));

      bool refining (TubeVector &x);
      double average_refining_threshold(const TubeVector &x, 
				vector<double>& slice_step, vector<double>& t_refining);
      double median_refining_threshold(const TubeVector &x, 
				vector<double>& slice_step, vector<double>& t_refining);
      bool refining_with_threshold(TubeVector & x, int nb_slices);
      void refining_all_slices(TubeVector & x, int nb_slices);

      void bisection_guess (TubeVector & x, TFnc& f);
      std::pair<int,std::pair<double,double>> bisection_guess(TubeVector x, TubeVector v, DynCtc* slice_ctr, TFnc& fnc, int variant);
      std::pair<int,std::pair<double,double>> bisection_guess(TubeVector x, TubeVector v, DynCtc* slice_ctr, TFnc& fnc);

      ibex::Vector m_max_thickness = ibex::Vector(1);
      float m_refining_fxpt_ratio = 0.005;
      float m_propa_fxpt_ratio = 0.0 ; // one call no fix point
      float m_var3b_fxpt_ratio = -1;   // no var3b
      float m_var3b_propa_fxpt_ratio = 0.0;
      /* Internal parameters for var3b algorithm */
      float m_var3b_bisection_minrate = 0.0001;
      float m_var3b_bisection_maxrate = 0.4;
      int m_var3b_bisection_ratefactor=2;
      void print_solutions(const list<TubeVector> & l_solutions);
      void display_solutions(const list<TubeVector> & l_solutions);
      int m_var3b_timept=0;
      int m_bisection_timept=0;
      int m_trace=0;
      int m_max_slices=5000;
      int m_refining_mode=0;
      int m_contraction_mode=0; 
      int m_stopping_mode=0;
      bool m_var3b_external_contraction=true;
 
     
      /* number of bisections */
      int bisections=0; 

      // Embedded graphics
      VIBesFigTubeVector *m_fig = NULL;
  };
}

#endif
