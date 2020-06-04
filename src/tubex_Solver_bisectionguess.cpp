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


#include <time.h>
#include "tubex_Solver.h"
#include <cmath>
#include<float.h>



using namespace std;
using namespace ibex;

namespace tubex
{
  void Solver::bisection_guess (TubeVector & x, tubex::Fnc& f){
    TubeVector v = f.eval_vector(x);
    CtcDynCid ctc(f);
    std::pair< int,std::pair<double,double> > guess = bisection_guess (x,v,&ctc,f,2);
    cout << " var " << guess.first << " time  " << guess.second.first << " point " << guess.second.second <<  endl;
    if (guess.first >= 0 && x.size() > 1) {
      cout << " var 0  " <<  *(x[0].first_slice()) <<  endl;
      cout << 	*(x[0].last_slice()) << endl;
      cout << " var 1 " <<  *(x[1].first_slice()) <<  endl;
      cout << 	*(x[1].last_slice()) << endl;


    }
   
  }

  

  std::pair<int,std::pair<double,double>> Solver::bisection_guess(TubeVector x, TubeVector v, Ctc* slice_ctr, tubex::Fnc& fnc, int variant){
    //variant 0 -> return immediately as soon as we find a potential gate
		//variant 1 -> return the largest gate in a slice
		//variant 2 -> return the largest gate in the complete tube

		/*variable - time of bisection - bisection point*/
		pair <int, pair <double,double> > bisection;

		double t_bisection;  //time of bisection (in t)
		double x_bisection; // value of bisection (in x)

		/*init pair*/
		bisection = make_pair(-1,make_pair(-1,-1));

		/*check if everything is ok*/
		assert(x.size() == v.size());
		assert(x.domain() == v.domain());
		assert(TubeVector::same_slicing(x, v));

		/*init all the tubes*/
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;
		/*auxiliar tubes*/
		vector<Slice*> aux_x_slice; TubeVector aux_x = x;
		vector<Slice*> aux_v_slice; TubeVector aux_v = v;


		double max_diameter = -1;
		double gate_diam;
		for (int it = 0 ; it < 2 ; it++){

			/*push slices for forward phase*/
			x_slice.clear(); v_slice.clear();
			aux_x_slice.clear(); aux_v_slice.clear();

			//for forward
			TPropagation t_propa;
			if (it == 0){
				t_propa = FORWARD;
				for (int i = 0 ; i < x.size() ; i++){
					x_slice.push_back(x[i].first_slice());
					v_slice.push_back(v[i].first_slice());
					aux_x_slice.push_back(aux_x[i].first_slice());
					aux_v_slice.push_back(aux_v[i].first_slice());
				}
			}
			//for backward
			else{
				t_propa = BACKWARD;
				for (int i = 0 ; i < x.size() ; i++){
					x_slice.push_back(x[i].last_slice());
					v_slice.push_back(v[i].last_slice());
					aux_x_slice.push_back(aux_x[i].last_slice());
					aux_v_slice.push_back(aux_v[i].last_slice());
				}
			}

			while (x_slice[0] != NULL){
				for (int i = 0 ; i < x.size() ; i++){
					if (t_propa & FORWARD){
						x_bisection = aux_x_slice[i]->output_gate().mid();
						t_bisection = aux_x_slice[i]->domain().ub();
						gate_diam = aux_x_slice[i]->output_gate().diam();
						aux_x_slice[i]->set_output_gate(x_bisection);
					}
					else if (t_propa & BACKWARD){
						x_bisection = aux_x_slice[i]->input_gate().mid();
						t_bisection = aux_x_slice[i]->domain().lb();
						gate_diam = aux_x_slice[i]->input_gate().diam();
						aux_x_slice[i]->set_input_gate(x_bisection);
					}
					if(dynamic_cast <CtcDynCid*> (slice_ctr)){
						CtcDynCid * cid = dynamic_cast <CtcDynCid*> (slice_ctr);
						cid->contract(aux_x_slice,aux_v_slice,t_propa);
					}

					else if(dynamic_cast <CtcDynCidGuess*> (slice_ctr)){
						CtcDynCidGuess * cidguess = dynamic_cast <CtcDynCidGuess*> (slice_ctr);
						cidguess->contract(aux_x_slice,aux_v_slice,t_propa);
					}
					else if(dynamic_cast <CtcDynBasic*> (slice_ctr)){
						CtcDynBasic * basic = dynamic_cast <CtcDynBasic*> (slice_ctr);
						basic->contract(aux_x_slice,aux_v_slice,t_propa);
					}

					for (int k = 0 ; k < aux_x_slice.size() ; k++){
						if (aux_x_slice[k]->is_empty()){
							if (variant == 0){
								bisection.first = i;
								bisection.second.first = t_bisection;
								bisection.second.second = x_bisection;
								return bisection;
							}
							else {
								if (gate_diam > max_diameter){
									bisection.first = i;
									bisection.second.first = t_bisection;
									bisection.second.second = x_bisection;
									max_diameter = gate_diam;
									break;
								}
							}
						}
					}
					//restore domains
					for (int k = 0 ; k < aux_x_slice.size() ; k++){
						aux_x_slice[k]->set_envelope(x_slice[k]->codomain());
						aux_x_slice[k]->set_input_gate(x_slice[k]->input_gate());
						aux_x_slice[k]->set_output_gate(x_slice[k]->output_gate());
					}
				}

				if (variant == 1){
					if (bisection.first != -1)
						return bisection;
				}

				if (t_propa & FORWARD){
					for (int i = 0 ; i < x.size() ; i++){
						x_slice[i] = x_slice[i]->next_slice();
						v_slice[i] = v_slice[i]->next_slice();
						aux_x_slice[i] = aux_x_slice[i]->next_slice();
						aux_v_slice[i] = aux_v_slice[i]->next_slice();
					}
				}
				else if (t_propa & BACKWARD){
					for (int i = 0 ; i < x.size() ; i++){
						x_slice[i] = x_slice[i]->prev_slice();
						v_slice[i] = v_slice[i]->prev_slice();
						aux_x_slice[i] = aux_x_slice[i]->prev_slice();
						aux_v_slice[i] = aux_v_slice[i]->prev_slice();
					}
				}
			}
		}
		/*if variant 2 is selected, it will return here*/
		return bisection;
	}





  std::pair<int,std::pair<double,double>> Solver::bisection_guess(TubeVector x, TubeVector v, Ctc* slice_ctr, tubex::Fnc& fnc){

		/*variable - time of bisection - bisection point*/
		pair <int, pair <double,double> > bisection;

		double t_bisection;  //time of bisection (in t)
		double x_bisection; // value of bisection (in x)

		/*init pair*/
		bisection = make_pair(-1,make_pair(-1,-1));

		/*check if everything is ok*/
		assert(x.size() == v.size());
		assert(x.domain() == v.domain());
		assert(TubeVector::same_slicing(x, v));



		/*init all the tubes*/
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;
		/*auxiliar tubes*/
		vector<Slice*> aux_x_slice; TubeVector aux_x = x;
		vector<Slice*> aux_v_slice; TubeVector aux_v = v;


		for (int it = 0 ; it < 2 ; it++){

			/*push slices for forward phase*/
			x_slice.clear(); v_slice.clear();
			aux_x_slice.clear(); aux_v_slice.clear();

			//for forward
			TPropagation t_propa;
			if (it == 0){
				t_propa = FORWARD;
				for (int i = 0 ; i < x.size() ; i++){
					x_slice.push_back(x[i].first_slice());
					v_slice.push_back(v[i].first_slice());
					aux_x_slice.push_back(aux_x[i].first_slice());
					aux_v_slice.push_back(aux_v[i].first_slice());
				}
			}
			//for backward
			else{
				t_propa = BACKWARD;
				for (int i = 0 ; i < x.size() ; i++){
					x_slice.push_back(x[i].last_slice());
					v_slice.push_back(v[i].last_slice());
					aux_x_slice.push_back(aux_x[i].last_slice());
					aux_v_slice.push_back(aux_v[i].last_slice());
				}
			}

			while (x_slice[0] != NULL){
				for (int i = 0 ; i < x.size() ; i++){
					if (t_propa & FORWARD){
						x_bisection = aux_x_slice[i]->output_gate().mid();
						t_bisection = aux_x_slice[i]->domain().ub();
						aux_x_slice[i]->set_output_gate(x_bisection);
					}
					else if (t_propa & BACKWARD){
						x_bisection = aux_x_slice[i]->input_gate().mid();
						t_bisection = aux_x_slice[i]->domain().lb();
						aux_x_slice[i]->set_input_gate(x_bisection);
					}
					if(dynamic_cast <CtcDynCid*> (slice_ctr)){
						CtcDynCid * cid = dynamic_cast <CtcDynCid*> (slice_ctr);
						cid->contract(aux_x_slice,aux_v_slice,t_propa);
					}

					else if(dynamic_cast <CtcDynCidGuess*> (slice_ctr)){
						CtcDynCidGuess * cidguess = dynamic_cast <CtcDynCidGuess*> (slice_ctr);
						cidguess->contract(aux_x_slice,aux_v_slice,t_propa);

					}
					else if(dynamic_cast <CtcDynBasic*> (slice_ctr)){
						CtcDynBasic * basic = dynamic_cast <CtcDynBasic*> (slice_ctr);
						basic->contract(aux_x_slice,aux_v_slice,t_propa);

					}

					for (int i = 0 ; i < aux_x_slice.size() ; i++){
						if (aux_x_slice[i]->is_empty()){
							bisection.first = i;
							bisection.second.first = t_bisection;
							bisection.second.second = x_bisection;
							return bisection;
						}
					}
					//restore domains
					for (int i = 0 ; i < aux_x_slice.size() ; i++){
						aux_x_slice[i]->set_envelope(x_slice[i]->codomain());
						aux_x_slice[i]->set_input_gate(x_slice[i]->input_gate());
						aux_x_slice[i]->set_output_gate(x_slice[i]->output_gate());
					}

				}

				if (t_propa & FORWARD){
					for (int i = 0 ; i < x.size() ; i++){
						x_slice[i] = x_slice[i]->next_slice();
						v_slice[i] = v_slice[i]->next_slice();
						aux_x_slice[i] = aux_x_slice[i]->next_slice();
						aux_v_slice[i] = aux_v_slice[i]->next_slice();
					}
				}
				else if (t_propa & BACKWARD){
					for (int i = 0 ; i < x.size() ; i++){
						x_slice[i] = x_slice[i]->prev_slice();
						v_slice[i] = v_slice[i]->prev_slice();
						aux_x_slice[i] = aux_x_slice[i]->prev_slice();
						aux_v_slice[i] = aux_v_slice[i]->prev_slice();
					}
				}
			}
		}
		/*there is not candidate to bisect (to check, the variable is -1)*/
		return bisection;
	}


}
