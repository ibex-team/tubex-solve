/* ============================================================================
 *  tubex-lib - CtcPicard class
 * ============================================================================
 *  Copyright : Copyright 2017 Simon Rohou
 *  License   : This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 *
 *  Author(s) : Simon Rohou
 *  Bug fixes : -
 *  Created   : 2015
 * ---------------------------------------------------------------------------- */

#ifndef CtcPicard_HEADER
#define CtcPicard_HEADER

#include "tubex_Ctc.h"
#include "tubex_Tube.h"
#include "tubex_TubeSlice.h"

namespace tubex
{
  /**
   * \brief CtcPicard class.
   */
  class CtcPicard : Ctc
  {
    public:

      CtcPicard(float delta = 1.1);

      // Tube
      bool contract(const ibex::Function& f, Tube& x);
      bool contract(const ibex::Function& f, std::vector<Tube*>& x);

      // TubeSlice
      bool contract(const ibex::Function& f, TubeSlice& x);
      bool contract(const ibex::Function& f, std::vector<TubeSlice*>& x);

    protected:

      bool contract(const ibex::Function& f,
                    ibex::Interval& x, const ibex::Interval& x0,
                    const ibex::Interval& h);

      bool contract(const ibex::Function& f,
                    ibex::IntervalVector& x, const ibex::IntervalVector& x0,
                    const ibex::Interval& h);

      float m_delta;
  };
}

#endif