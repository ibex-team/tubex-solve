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
//#include "ibex_BoolInterval.h"
using namespace std;
namespace tubex
{
  class Solver
  {
  public:

      Solver(const ibex::Vector& max_thickness);
      ~Solver();

      // Ratios used for stopping fixed point algorithms 
      void set_refining_fxpt_ratio(float refining_fxpt_ratio);
      void set_propa_fxpt_ratio(float propa_fxpt_ratio);
      void set_cid_fxpt_ratio(float cid_fxpt_ratio);
      void set_cid_propa_fxpt_ratio(float cid_propa_fxpt_ratio);

      // Slice choice for cid
      void set_cid_timept(int cid_timept); // where do the cid contraction : 1 domain.ub(); -1 domain.lb(); 0 max_diam_gate(); 2 randomly domain.ub() or domain.lb()

      // Slice choice for bisection
      void set_bisection_timept (int bisection_timept); // where to bisect : 1 domain.ub(); -1 domain.lb(); 0 max_diam_gate; 2 randomly domain.ub() or domain.lb()

      // no more refining when a number of slices greater than max_slices is already reached
      void set_max_slices(int max_slices); 

      // refining mode : which slices to refine
      void set_refining_mode(int refining_mode); // 0 : all slices ; 1 : one slice ; 2 : the slices with a difference more than average

      void set_trace(int trace);
      double solving_time;

      const std::list<TubeVector> solve(const TubeVector& x0, void (*ctc_func)(TubeVector&));
      VIBesFigTubeVector* figure();
      static const ibex::BoolInterval solutions_contain(const std::list<TubeVector>& l_solutions, const TrajectoryVector& truth);


  protected:
      bool empty_intersection(TubeVector& t1, TubeVector& t2);
      void clustering(std::list<std::pair<int,TubeVector> >& l_tubes);
      void clustering(std::list<TubeVector>& l_tubes);
      bool stopping_condition_met(const TubeVector& x);
      bool fixed_point_reached(double volume_before, double volume_after, float fxpt_ratio);
      void propagation(TubeVector &x, void (*ctc_func)(TubeVector&), float propa_fxpt_ratio);
      void cid(TubeVector &x, void (*ctc_func)(TubeVector&));
      bool refining (TubeVector &x);
      double refining_threshold(const TubeVector &x, 
				vector<double>& slice_step, vector<double>& t_refining);
      ibex::Vector m_max_thickness = ibex::Vector(1);
      float m_refining_fxpt_ratio = 0.005;
      float m_propa_fxpt_ratio = 0.005;
      float m_cid_fxpt_ratio = 0.005;
      float m_cid_propa_fxpt_ratio = 0.005;
      /* Internal parameters for cid algorithm */
      float m_cid_bisection_minrate = 0.0001;
      float m_cid_bisection_maxrate = 0.4;
      int m_cid_bisection_ratefactor=2;

      int m_cid_timept=0;
      int m_bisection_timept=0;
      int m_trace=0;
      int m_max_slices=5000;
      int m_refining_mode=0; 
      // Embedded graphics
      VIBesFigTubeVector *m_fig = NULL;
  };
}

#endif