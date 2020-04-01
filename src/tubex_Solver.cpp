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
#include "tubex_Exception.h"
#define GRAPHICS 0
#include <cmath>
#include<float.h>



using namespace std;
using namespace ibex;

namespace tubex
{
  Solver::Solver(const Vector& max_thickness)
  {
    m_max_thickness = max_thickness;
    
    #if GRAPHICS // embedded graphics
      vibes::beginDrawing();
      m_fig = new VIBesFigTubeVector("Solver");
      m_fig->set_properties(100, 100, 700, 500);
    #endif
  }

  Solver::~Solver()
  {
    #if GRAPHICS
      delete m_fig;
      vibes::endDrawing();
    #endif
  }

  void Solver::set_refining_fxpt_ratio(float refining_fxpt_ratio)
  {
    //    assert(Interval(0.,1.).contains(refining_fxpt_ratio));
    m_refining_fxpt_ratio = refining_fxpt_ratio;
  }

  void Solver::set_propa_fxpt_ratio(float propa_fxpt_ratio)
  {
    assert(Interval(0.,1.).contains(propa_fxpt_ratio));
    m_propa_fxpt_ratio = propa_fxpt_ratio;
  }

  void Solver::set_var3b_fxpt_ratio(float var3b_fxpt_ratio)
  {
    assert(Interval(0.,1.).contains(var3b_fxpt_ratio));
    m_var3b_fxpt_ratio = var3b_fxpt_ratio;
  }

   void Solver::set_var3b_propa_fxpt_ratio(float var3b_propa_fxpt_ratio)
  {
    assert(Interval(0.,1.).contains(var3b_propa_fxpt_ratio));
    m_var3b_propa_fxpt_ratio = var3b_propa_fxpt_ratio;
  }

  void Solver::set_var3b_timept(int var3b_timept)
  {
    m_var3b_timept = var3b_timept;
  }

  void Solver::set_bisection_timept(int bisection_timept)
  {
    m_bisection_timept = bisection_timept;
  }

  void Solver::set_trace(int trace)
  {
    m_trace = trace;
  }


  void Solver::set_max_slices(int max_slices)
  {
    m_max_slices=max_slices;
  }
  

  void Solver::set_refining_mode(int refining_mode)
  {
    m_refining_mode=refining_mode;
  }
    
  void Solver::set_contraction_mode(int contraction_mode)
  {
    m_contraction_mode=contraction_mode;
  }


  bool Solver::refining(TubeVector& x)
  { // no refining if max_slices is already reached
    int nb_slices=x[0].nb_slices();
    if (nb_slices >= m_max_slices)
      return false;

    //    cout << " volume before refining " << x.volume() << endl;

    if  (m_refining_mode == 1){  // one slice is refined
      
      // double t_refining = x[0].wider_slice()->domain().mid() // the widest slice            
      double t_refining= x.steepest_slice()->domain().mid();    
      x.sample(t_refining);
      // cout << "refining point " << t_refining << endl;
    }
    else if //(m_refining_mode==0 || x.volume() >= DBL_MAX)
      (m_refining_mode==0)
      { // all slices are refined 
	
	vector<double> t_refining;

	for (const Slice*s= x[0].first_slice(); s!=NULL; s=s->next_slice()){
	  t_refining.push_back(s->domain().mid());

	}
	//	cout << " refining " << t_refining.size() << endl;

	for (int k=0; k<t_refining.size(); k++){
	  if (k+nb_slices >= m_max_slices) break;
	  x.sample(t_refining[k]);
	}
  }
	    
    else if (m_refining_mode== 2 || m_refining_mode== 3){
      // the refining is focused on "slices " with a larger than average (or median) max difference (in all dimensions)  between input and output gates
      
      vector<double> t_refining;
      vector<double>  slice_step;

      double step_threshold = refining_threshold(x, slice_step, t_refining);
      //      cout << " step threshold " << step_threshold << endl;
     
      int new_slices=0;
      for (int k=0; k<t_refining.size(); k++){
	//       	cout << "refining " << t_refining[k] << " " << slice_step[k] << endl;
	if( (nb_slices+ new_slices) >= m_max_slices) break;
	if (slice_step[k] >= step_threshold)
	    {//cout << "refining+ " << t_refining[k] << endl;
	      x.sample(t_refining[k]); new_slices++;
	    }
      }
      if (nb_slices < m_max_slices && x[0].nb_slices() == nb_slices)  // patch to avoid infinite loops (the selected slices could not be bisected)
	for (int k=0; k<t_refining.size(); k++){
	  x.sample(t_refining[k]);
          nb_slices++;
	  if (nb_slices >= m_max_slices) break;
	}
    }
    return true;
  }


  double Solver::refining_threshold(const TubeVector &x, vector<double> & slice_step, vector<double>& t_refining) {
        int nbsteps=0;
	double step_threshold;
	vector<double> stepmed;
	const Slice* s[x.size()];
	for(int k=0; k< x.size(); k++)
	  s[k]=x[k].first_slice();

	for (const Slice*slice=s[0]; slice!=NULL; slice=slice->next_slice()){

	  t_refining.push_back(slice->domain().mid());
	  //	       cout << *slice << endl;
	  double step_max= fabs(slice->output_gate().mid() - slice->input_gate().mid());
	  //	  cout << " step  " <<  t_refining.size() << " " << step_max << " input " << slice->input_gate() << " output " << slice->output_gate() << endl;


	  for (int k=1; k< x.size(); k++){
	    //	    cout << " step k " << k << " " << fabs(s[k]->output_gate().mid() - s[k]->input_gate().mid()) << endl;
	    step_max=std::max(step_max,fabs(s[k]->output_gate().mid() - s[k]->input_gate().mid()));
	    //	    cout << " step_max " << k << "  " << step_max << endl;
	    s[k]=s[k]->next_slice();

	  }
	  //	  cout << "step_max " << step_max << endl;
	  slice_step.push_back(step_max);
	  if (m_refining_mode==3) stepmed.push_back(step_max); // storage for computing the median value
	  if (m_refining_mode==2)
	    if (step_max < DBL_MAX)  // to not take into account infinite gates in average computation
	      {
		nbsteps++;
		step_threshold=(step_threshold*(nbsteps-1)+step_max)/nbsteps;
	      }
	  
	  //	  cout << "step_threshold  " << step_threshold << endl;
	}
	if (m_refining_mode==3){
	  sort(stepmed.begin(),stepmed.end());
	  step_threshold =stepmed[stepmed.size()/2];
	}
	//       	cout << " step_threshold " << step_threshold << endl;
	//	cout << " step_max " << stepmed[stepmed.size()-1] << endl;
	return step_threshold;
  }

  const list<TubeVector> Solver::solve(const TubeVector& x0,  tubex::Function& f, void (*ctc_func)(TubeVector&)) { return (solve(x0,&f,ctc_func));}

  const list<TubeVector> Solver::solve(const TubeVector& x0, void (*ctc_func)(TubeVector&)) {return (solve(x0,NULL, ctc_func));}


  const list<TubeVector> Solver::solve(const TubeVector& x0, tubex::Function* f, void (*ctc_func)(TubeVector&))
  {
    solving_time=0.0;
    assert(x0.size() == m_max_thickness.size());

    int i = 0;
    clock_t t_start = clock();

    #if GRAPHICS
    m_fig->show(true);
    #endif
    double t_init=x0[0].domain().lb();
    int prev_level = 0;
    list<pair<pair<int,double>,TubeVector> > s;
    s.push_back(make_pair(make_pair(0,t_init), x0));
    list<TubeVector> l_solutions;

    while(!s.empty())
    {
      int level = s.front().first.first;
      double t_bisect= s.front().first.second;
      /*
      if(level != prev_level && s.size() >= 20)

      {
        cout << "clustering (" << s.size() << " items)" << endl;
	clustering(s);
        cout << "after clustering (" << s.size() << " items)" << endl;
        prev_level = level;
      }
	*/
      TubeVector x = s.front().second;
      s.pop_front();

      bool emptiness;
      double volume_before_refining;
      
      //      cout << " before propagation " << x << endl;
      //      propagation(x, ctc_func, m_propa_fxpt_ratio);
      if (level==0)
	propagation(x, f, ctc_func, m_propa_fxpt_ratio, false, t_bisect);
      else
	propagation(x, f, ctc_func, m_propa_fxpt_ratio, true, t_bisect);
      //      cout << " after propagation " << x << endl;
      
      emptiness = x.is_empty();
      double volume_before_var3b;
      if (!emptiness && m_var3b_fxpt_ratio){
	do
	  { 
	    // cout << " volume before var3b "  << x.volume() << endl;
	    //	    var3b(x, ctc_func);
	    var3b(x, f, ctc_func);
	    // cout << " volume after var3b "  << x.volume() << endl;
	    emptiness = x.is_empty();
	  }
	while(!emptiness  
	      && !fixed_point_reached(volume_before_var3b, x.volume(), m_var3b_fxpt_ratio));
      }

      
      emptiness = x.is_empty();
      if (! emptiness)
	do
      {
        volume_before_refining = x.volume();
	//  cout << " volume before refining " <<  volume_before_refining << endl;
        // 1. Refining

	if(m_refining_fxpt_ratio != 0.)
	  if (! refining(x))
	    {
	      cout << " end refining " <<  volume_before_refining << " after  " << x.volume() << endl; 
	      break;}
	cout << " nb_slices after refining step " << x[0].nb_slices() << endl;
	// 2. Propagations up to the fixed point
	//	propagation(x, ctc_func, m_propa_fxpt_ratio);
	propagation(x, f, ctc_func, m_propa_fxpt_ratio, false, x[0].domain().lb());
	// 3.      
	emptiness = x.is_empty();
	double volume_before_var3b;
	if (!emptiness && m_var3b_fxpt_ratio){
	  do
	    { 
	      // cout << " volume before var3b "  << x.volume() << endl;
	      //	      var3b(x, ctc_func);
	      var3b(x, f, ctc_func);
	      // cout << " volume after var3b "  << x.volume() << endl;
	      emptiness = x.is_empty();
	    }
	  while(!emptiness  
		&& !fixed_point_reached(volume_before_var3b, x.volume(), m_var3b_fxpt_ratio));
	  }
       	cout << " volume after refining " <<  x.volume() << endl;
      }
      
      while(!emptiness
	    && !stopping_condition_met(x)
	    && ( x.volume() >= DBL_MAX  || !fixed_point_reached(volume_before_refining, x.volume(), m_refining_fxpt_ratio)));

      // 4. Bisection
      emptiness=x.is_empty();
      if(!emptiness)
        {
          if(stopping_condition_met(x) || m_bisection_timept==-2 )
          {
            l_solutions.push_back(x);
	    /*
            #if GRAPHICS // displaying solution

	      ostringstream o; o << "solution_" << i;
	      m_fig->add_tubevector(&l_solutions.back(), o.str());
              m_fig->show(true);
            #endif
	    */
              i++;
	      if (m_trace) cout << "solution_" << i <<  " vol  " << x.volume() << " max thickness " << x.max_diam() << endl;
          }

          else
          {
            cout << "Bisection... (level " << level << ")" << endl;
	    double t_bisection;
	    
	    if (m_bisection_timept==0)
	      x.max_gate_diam(t_bisection);
	    else if (m_bisection_timept==1)
	      t_bisection=x[0].domain().ub();
	    else if  (m_bisection_timept==-1)
	      t_bisection=x[0].domain().lb();
	    else if  (m_bisection_timept==2){
	      if (rand()%2)
		t_bisection=x[0].domain().lb();
	      else
		t_bisection=x[0].domain().ub();
	    }
	      else if  (m_bisection_timept==3){
	      if (level%2)
		t_bisection=x[0].domain().lb();
	      else
		t_bisection=x[0].domain().ub();
	    }

	 


	    level++;
            try{
	      pair<TubeVector,TubeVector> p_x = x.bisect(t_bisection);
	  
	      //   s.push_back(make_pair(level, p_x.second));   // breadth first variant
	      //   s.push_back (make_pair(level, p_x.first));
	    
	      s.push_front(make_pair(make_pair(level,t_bisection), p_x.second));       // depth first variant
	      s.push_front(make_pair(make_pair(level,t_bisection), p_x.first));

	    }
	    catch (Exception &)   // when the bisection time was not bisectable, change to largest_slice
	      {	 
		
                t_bisection=x.largest_slice()->domain().mid();
	        pair<TubeVector,TubeVector> p_x = x.bisect(t_bisection);

	    //   s.push_back(make_pair(level, p_x.second));   // breadth first variant
	    //   s.push_back (make_pair(level, p_x.first));
	    
		s.push_front(make_pair(make_pair(level,t_bisection), p_x.second));       // depth first variant
		s.push_front(make_pair(make_pair(level,t_bisection), p_x.first));
	      }


	    cout << " t_bisection " << t_bisection << " x volume " << x.volume() << " nb_slices " << x.nb_slices() << endl;
             
     	  }
    	}
	  
    }
    //       cout << "Solutions: " << l_solutions.size() << "  (" << (int)((double)(clock() - t_start)/CLOCKS_PER_SEC) << "s)   " << flush;
    
    if (m_trace){
      cout << endl;
      printf("Solving time : %.2fs\n", (double)(clock() - t_start)/CLOCKS_PER_SEC);
    }

    int j = 0;
    list<TubeVector>::iterator it;
    for(it = l_solutions.begin(); it != l_solutions.end(); ++it)
    {
      j++;
      
      if (m_trace) 
	cout << "  " << j << ": "
           << *it <<  ", ti↦" << (*it)(it->domain().lb()) <<  ", tf↦" << (*it)(it->domain().ub())
           << " (max diam: " << it->max_diam() << ")"
	   << " volume : " << it->volume() 
           << endl;
      
    }
    while (l_solutions.size()>1){
      int k = l_solutions.size();
      clustering(l_solutions);
      if (k==l_solutions.size())
	{ if (m_trace) cout << " end of clustering " << endl;
	  break;}
      j=0;
      for(it = l_solutions.begin(); it != l_solutions.end(); ++it)
	{
	  j++;
	  if (m_trace)
	    cout << "  " << j << ": "
	       << *it << ", ti↦" << (*it)(it->domain().lb())  <<  ", tf↦" << (*it)(it->domain().ub())
	       << " (max diam: " << it->max_diam() << ")"
	       << " volume : " << it->volume() 
	       << endl;
	}

    }

    for(it = l_solutions.begin(); it != l_solutions.end(); ++it){
            #if GRAPHICS // displaying solution
              i++;
	      ostringstream o; o << "solution_" << i;
	      const TubeVector* tv= &(*it);
	      m_fig->add_tubevector(tv, o.str());
              m_fig->show(true);
            #endif
	    
	      }
    double total_time =  (double)(clock() - t_start)/CLOCKS_PER_SEC;
    solving_time=total_time;
    if (m_trace)    printf("Total time with clustering: %.2fs\n", (double)(clock() - t_start)/CLOCKS_PER_SEC);
    return l_solutions;
    }
  
  // clustering during the search :  not used 
  void Solver::clustering(list<pair<int,TubeVector> >& l_tubes)
  {
    assert(!l_tubes.empty());
    list<pair<int,TubeVector> > l_clustered;
    list<pair<int,TubeVector> >::iterator it1, it2;
    int size= l_tubes.size();
    for(it1 = l_tubes.begin(); it1 != l_tubes.end(); ++it1)     
      {
	bool clustering = false;
	for(it2 = l_clustered.begin(); it2 != l_clustered.end(); ++it2)
	  {
	    if(
	       !((it1->second & it2->second).is_empty())
	       )
	      {
		it2->second = it2->second | it1->second;
		clustering = true;
	      }
	  }
	if(!clustering)
	  l_clustered.push_back(*it1);
      }
    
    l_tubes = l_clustered;
  }


 /* Clustering algorithm for a list of tubes : modify the list by merging tubes that 
    have a non empty intersection 
    fixed point algorithm until no tube can be merged */

  void Solver::clustering(list<TubeVector>& l_tubes)  {
    //  assert(!l_tubes.empty());
    list<TubeVector> l_clustered;
    list<TubeVector>::iterator it1, it2;
    int nn=0;
    for(it1 = l_tubes.begin(); it1 != l_tubes.end(); ++it1) {
      bool clustering = false;
      int nn1=0;
      for(it2 = l_clustered.begin(); it2 != l_clustered.end(); ++it2) {
		if(!((*it1 & *it2).is_empty()))
	//	if (! empty_intersection(*it1,*it2)) // TO DO : test this algo 
	  {
	    *it2 = (*it1 | *it2);
	    clustering = true;
	    //	cout << " tube " << nn+1 << " in cluster " << nn1+1 << endl;
	    //	cout << ", ti↦" << (*it2)(it2->domain().lb()) << endl;
	    if (m_trace) cout << " tube " << nn+1 << " in cluster " << nn1+1 << endl;
	    break;
	  }
	  nn1++;
      }
      if(!clustering){
	if (m_trace) {cout << "new cluster tube " << nn+1 << endl;
	  cout << ", ti↦" << (*it1)(it1->domain().lb()) << endl;}
	l_clustered.push_back(*it1);
      }
      nn++;
    }
    l_tubes = l_clustered;
  }

  /*   more efficient algorithm in case of tubes with different slicings : no need to compute  the tube intersection  TODO test this algorithm */

  bool Solver::empty_intersection( TubeVector & tubevector1, TubeVector & tubevector2)
  {

      for (int k=0; k< tubevector1.size() ; k++)
         
	 {Slice* s2=tubevector2[k].first_slice();
	   Slice* s1=tubevector1[k].first_slice();
	   if ((s2->input_gate()&s1->codomain()).is_empty()
	       ||
	       (s1->input_gate()&s2->codomain()).is_empty())
	     {cout << "first slice " << endl; 
	       return true;}
	     
	   for (const Slice* s= tubevector1[k].first_slice(); s!=NULL; s=s->next_slice()){

	     while (s2->domain().ub() <= s->domain().lb()){
          
	       s2=s2->next_slice();}
             
	     if ((s2->input_gate()&s->codomain()).is_empty())
	       {cout << " middle s2 " << endl; 
				 cout << *s << endl;
		                 cout << *s2 << endl; 
		 return true;}
	     if (s->domain().ub() <= s2->domain().lb())
	      if ((s->output_gate() & s2->codomain()).is_empty())
		{cout << " middle s "<< endl;  return true;}
	   
	     s2=tubevector2[k].last_slice();
	     s1=tubevector1[k].last_slice();
	     if ((s2->output_gate()&s1->codomain()).is_empty()
		 ||
		 (s1->output_gate()&s2->codomain()).is_empty())
	       {cout << " last slice " << endl; 
		 return true;}
	   
	   }
	 }
  return false;
  }


	
	  
  VIBesFigTubeVector* Solver::figure()
  {
    return m_fig;
  }

  bool Solver::stopping_condition_met(const TubeVector& x)
  {
    assert(x.size() == m_max_thickness.size());
    Vector x_max_thickness = x.max_diam();
    for(int i = 0 ; i < x.size() ; i++)
      {
	//      double tmaxgate;
	//	if(x[i].max_gate_diam(tmaxgate) > m_max_thickness[i])
	if(x_max_thickness[i] > m_max_thickness[i])
        return false;
      }
    /*
    for(int i = 0 ; i < x.size() ; i++)
      cout << "thickness " << i << "  " <<  x_max_thickness[i] << endl;
    */
    return true;
  }

  bool Solver::fixed_point_reached(double volume_before, double volume_after, float fxpt_ratio)
  {
    if (fxpt_ratio > 1) return false;
    if(fxpt_ratio == 0. || volume_after == volume_before)
      return true;

    assert(Interval(0.,1.).contains(fxpt_ratio));
    int n = m_max_thickness.size();
    return (std::pow(volume_after, 1./n) / std::pow(volume_before, 1./n)) >= fxpt_ratio;
  }
  /*
  void Solver::propagation(TubeVector &x, void (*ctc_func)(TubeVector&), float propa_fxpt_ratio)
  {
    assert(Interval(0.,1.).contains(propa_fxpt_ratio));

    if(propa_fxpt_ratio == 0.)
      return;
    
    bool emptiness;
    double volume_before_ctc;
    int nb_iter=0;
    do
    {
      volume_before_ctc = x.volume();
      //      cout << " volume before ctc " << volume_before_ctc << endl;
      ctc_func(x);
      //      cout << "volume after ctc " << x.volume() << endl;
      emptiness = x.is_empty();
      nb_iter++;
    } 
    while(!emptiness
	  //          && !stopping_condition_met(x)
         && !fixed_point_reached(volume_before_ctc, x.volume(), propa_fxpt_ratio));

  }
  */

  void Solver::propagation(TubeVector &x, tubex::Function* f, void (*ctc_func)(TubeVector&), float propa_fxpt_ratio, bool incremental, double t0, TPropagation propa )
  {


    if (f){
      CtcPicard ctc_picard;
  
      ctc_picard.preserve_slicing(true);
      if (x.volume()>= DBL_MAX) 
	ctc_picard.contract(*f, x, FORWARD | BACKWARD);
    }


    if(propa_fxpt_ratio == 0.)
      return;

    

    double volume_before_ctc;

    Ctc* ctc_dyncid;
    CtcIntegration* ctc_integration;
    if (f){
      if (m_contraction_mode==1) 
	ctc_dyncid= new CtcDynCid(*f);     
      else
	ctc_dyncid =  new CtcDynCidGuess(*f);
      ctc_dyncid->set_fast_mode(true);
      ctc_integration = new CtcIntegration (*f,ctc_dyncid);
    }
    do
      {

	volume_before_ctc = x.volume();
        if (ctc_func) ctc_func(x);  // Other constraints contraction
	if (f){                     // ODE contraction
	  if (m_contraction_mode==0){   // CtcDeriv

    
	  CtcDeriv ctc;

      //      cout << " volume before ctc " << volume_before_ctc << endl;
	  ctc.set_fast_mode(true);

	  ctc.contract(x, f->eval_vector(x),FORWARD | BACKWARD);

	  //      cout << "volume after ctc " << x.volume() << endl;

	  }
    
	  else{                           // CtcIntegration
	    




	    TubeVector v = f->eval_vector(x);

	
	    if (incremental && ctc_func==NULL){ // incrementality if first call, incremental and  
	      ctc_integration->set_incremental_mode(true);
	      ctc_integration->contract(x,v,t0,FORWARD) ;
	      ctc_integration->contract(x,v,t0,BACKWARD) ;
	    }
	    else{
	      ctc_integration->set_incremental_mode(false);
	      ctc_integration->contract(x,v,x[0].domain().lb(),FORWARD) ;
	      ctc_integration->contract(x,v,x[0].domain().ub(),BACKWARD) ;

	    }

	  }
	}
      }
    while(!(x.is_empty())
	  && !fixed_point_reached(volume_before_ctc, x.volume(), propa_fxpt_ratio));
    if (ctc_dyncid) delete ctc_dyncid;
  }





  






  //  void Solver::var3b(TubeVector &x, void (*ctc_func)(TubeVector&))
  void Solver::var3b(TubeVector &x, tubex::Function * f,void (*ctc_func)(TubeVector&) )
  {
    if(m_var3b_fxpt_ratio == 0. || x.volume() >= DBL_MAX)
      return;

    double t_bisection;

    if (m_var3b_timept==1)
      t_bisection=x[0].domain().ub();
    else if (m_var3b_timept==-1)
      t_bisection=x[0].domain().lb();
    else  if (m_var3b_timept==2){
      if (rand()%2)
	t_bisection=x[0].domain().lb();
      else
	t_bisection=x[0].domain().ub();
    }
    else
      x.max_gate_diam(t_bisection);  
    //    cout << " t_bisection var3b " << t_bisection << endl;
    for(int k=0; k<x.size() ; k++)
      {

	double rate=m_var3b_bisection_minrate;

	while (rate < m_var3b_bisection_maxrate){
	  try
	    {pair<TubeVector,TubeVector> p_x = x.bisect(t_bisection,k,rate);

	     TubeVector branch_x= p_x.first;
	     //	     propagation(branch_x, ctc_func, m_propa_fxpt_ratio);
	     propagation(branch_x, f, ctc_func, m_propa_fxpt_ratio, true, t_bisection);
	     if (branch_x.is_empty())
	       x=p_x.second;
	     else {x = p_x.second | branch_x ; break;}
	     rate= m_var3b_bisection_ratefactor*rate;
	    }
	  
	  catch (Exception& )
	    {break;}

	}
	//	propagation(x,ctc_func, m_propa_fxpt_ratio);
	propagation(x,f, ctc_func, m_propa_fxpt_ratio, true, t_bisection);
	rate = 1 - m_var3b_bisection_minrate;
       
	while (rate > 1-m_var3b_bisection_maxrate ){
	  try{
	    pair<TubeVector,TubeVector> p_x = x.bisect(t_bisection,k,rate);
	     TubeVector branch_x= p_x.second;
	     //	     propagation(branch_x, ctc_func, m_propa_fxpt_ratio);
	     propagation(branch_x, f, ctc_func, m_propa_fxpt_ratio, true, t_bisection);
	     if (branch_x.is_empty())
	       x=p_x.first;
	     else {x = p_x.first | branch_x ; break;}
	     rate=1-m_var3b_bisection_ratefactor*(1-rate);
	  }
	  catch (Exception& )
	    {break;}

	}
	//	propagation(x,ctc_func, m_propa_fxpt_ratio);
	propagation(x,f , ctc_func, m_propa_fxpt_ratio, true, t_bisection);
      }


  }



 const BoolInterval Solver::solutions_contain(const list<TubeVector>& l_solutions, const TrajectoryVector& truth)
  {
    assert(!l_solutions.empty());

    BoolInterval result = NO;
    list<TubeVector>::const_iterator it;
    for(it = l_solutions.begin() ; it != l_solutions.end() ; ++it)
    {
      assert(truth.size() == it->size());
      BoolInterval b = it->contains(truth);
      if(b == YES) return YES;
      else if(b == MAYBE) result = MAYBE;
    }
    
    return result;
  }
}
