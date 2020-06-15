
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
    m_refining_fxpt_ratio = refining_fxpt_ratio;
  }

  void Solver::set_propa_fxpt_ratio(float propa_fxpt_ratio)
  {
    m_propa_fxpt_ratio = propa_fxpt_ratio;
  }

  void Solver::set_var3b_fxpt_ratio(float var3b_fxpt_ratio)
  {
    m_var3b_fxpt_ratio = var3b_fxpt_ratio;
  }

  
   void Solver::set_var3b_propa_fxpt_ratio(float var3b_propa_fxpt_ratio)
  {
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

 

  double Solver::one_finite_gate(TubeVector &x){
    bool finite=true;
    for (int i=0; i< x.size() ; i++)
      if (x[i].first_slice()->input_gate().diam() >= DBL_MAX)
	{finite=false;break;}
    if (finite==true)
      return x[0].first_slice()->domain().lb();
    else{
      finite=true; 
      for (int i=0; i< x.size() ; i++)
	if (x[i].last_slice()->output_gate().diam() >= DBL_MAX)
	{finite=false;break;}
    }
    if (finite==true)
      return x[0].last_slice()->domain().ub();
    else{
      double t0= x[0].domain().lb()-1;
      for (int i=0 ;i< x.size(); i++){
	for ( Slice*s= x[i].first_slice(); s!=NULL; s=s->next_slice())
	  if (s->input_gate().diam()<DBL_MAX && s->input_gate().diam()>0  )
	    {return s->domain().lb();}
      }
    }
    throw tubex::Exception (" solver " , "solver unable to bisect ");

  }
	       
    

  bool Solver::refining(TubeVector& x)
  {
    int nb_slices=x[0].nb_slices();
    if (nb_slices >= m_max_slices)  // no refining if max_slices is already reached
      return false;

    //   cout << " volume before refining " << x.volume() << endl;

    if  (m_refining_mode == 1)
      {  // one slice is refined
      
      // double t_refining = x[0].wider_slice()->domain().mid() // the widest slice            
      double t_refining= x.steepest_slice()->domain().mid();    
      x.sample(t_refining);
      // cout << "refining point " << t_refining << endl;
      return true;
    }
    else if 
      (m_refining_mode==0)
      { // all slices are refined 

	
	for (int i=0; i<x.size();i++)
	  {
	    int k=0;
	    for ( Slice*s= x[i].first_slice(); s!=NULL; s=s->next_slice()){
	      if (k+nb_slices >= m_max_slices) break;
	      x[i].sample(s->domain().mid(),s);
	      s=s->next_slice();
	      k++;
	    }
	  }
	
	/*
	int k=0;
	for (const Slice*s= x[0].first_slice(); s!=NULL; s=s->next_slice()){
	  if (k+nb_slices >= m_max_slices) break;
	  x.sample(s->domain().mid());
	  s=s->next_slice();
	  k++;
	}
	*/
	return true;
      }
    else if (m_refining_mode== 2 || m_refining_mode== 3)
      return refining_with_threshold(x, nb_slices);

   
  }

  // the refining is focused on "slices " with a larger than average (or median) max difference (in all dimensions)  between input and output gates
  bool Solver::refining_with_threshold (TubeVector & x ,int nb_slices){
    
      vector<double> t_refining;
      vector<double>  slice_step;

      double step_threshold;
      if (m_refining_mode==2) step_threshold= average_refining_threshold(x, slice_step, t_refining);
      else if (m_refining_mode==3) step_threshold= median_refining_threshold(x, slice_step, t_refining);
      //cout << " step threshold " << step_threshold << " nb_slices " << nb_slices << endl;
      
      for (int i=0; i<x.size();i++)
	  {
	    int k=0; int new_slices=0;
	    for ( Slice*s= x[i].first_slice(); s!=NULL; s=s->next_slice()){
	      if (new_slices +nb_slices >= m_max_slices) break;
	      if (slice_step[k] >= step_threshold){
		x[i].sample(s->domain().mid(),s);
		s=s->next_slice();
		new_slices++;
	      }
	      k++;
	    }
	  }

        
      if (nb_slices < m_max_slices && x[0].nb_slices() == nb_slices)  // patch to avoid infinite loops (the selected slices could not be bisected)
	for (int k=0; k<t_refining.size(); k++){
	  //cout << " sample " << k << endl;
	  x.sample(t_refining[k]);
          nb_slices++;
	  if (nb_slices >= m_max_slices) break;
	}
      //      cout << " after sample " << nb_slices << "  " << x[0].nb_slices() << endl;
      return true;
  }

  



  double Solver::median_refining_threshold (const TubeVector &x, vector<double> & slice_step, vector<double>& t_refining) {

	double step_threshold;
	vector<double> stepmed;
	const Slice* s[x.size()];
	for(int k=0; k< x.size(); k++)
	  s[k]=x[k].first_slice();

	for (const Slice*slice=s[0]; slice!=NULL; slice=slice->next_slice()){

	  t_refining.push_back(slice->domain().mid());
	  double step_max= fabs(slice->output_gate().mid() - slice->input_gate().mid());

	  for (int k=1; k< x.size(); k++){

	    step_max=std::max(step_max,fabs(s[k]->output_gate().mid() - s[k]->input_gate().mid()));

	    s[k]=s[k]->next_slice();
	  }
	  
	  slice_step.push_back(step_max);
	  if (step_max < DBL_MAX)  // to not take into account infinite gates in average computation
	    {
	      stepmed.push_back(step_max); // storage for computing the median value
	    }

	}
	
	sort(stepmed.begin(),stepmed.end());
	step_threshold =stepmed[stepmed.size()/2];
	return step_threshold;
  }

 

  double Solver::average_refining_threshold(const TubeVector &x, vector<double> & slice_step, vector<double>& t_refining) {
        int nbsteps=0;
	double step_threshold;

	const Slice* s[x.size()];
	for(int k=0; k< x.size(); k++)
	  s[k]=x[k].first_slice();

	for (const Slice*slice=s[0]; slice!=NULL; slice=slice->next_slice()){
	  t_refining.push_back(slice->domain().mid());
	  double step_max= fabs(slice->output_gate().mid() - slice->input_gate().mid());

	  for (int k=1; k< x.size(); k++){

	    step_max=std::max(step_max,fabs(s[k]->output_gate().mid() - s[k]->input_gate().mid()));

	    s[k]=s[k]->next_slice();
	  }

	  slice_step.push_back(step_max);
	  if (step_max < DBL_MAX)  // to not take into account infinite gates in average computation
	    {
		nbsteps++;
		step_threshold=(step_threshold*(nbsteps-1)+step_max)/nbsteps;
	    }

	}
	return step_threshold;
  }


  
  const list<TubeVector> Solver::solve(const TubeVector& x0,  tubex::Fnc& f, void (*ctc_func)(TubeVector&, double& t0, bool incremental)) { return (solve(x0,&f,ctc_func));}
  

  const list<TubeVector> Solver::solve(const TubeVector& x0, void (*ctc_func)(TubeVector&,double& t0, bool incremental)) {return (solve(x0,NULL, ctc_func));}

 
  const list<TubeVector> Solver::solve(const TubeVector& x0, tubex::Fnc* f, void (*ctc_func)(TubeVector&, double& t0, bool incremental))
  {
    bisections=0;
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
      if (level==0) // && x.volume()> 1.e100)
	propagation(x, f, ctc_func, m_propa_fxpt_ratio, false, t_bisect);
      else
	propagation(x, f, ctc_func, m_propa_fxpt_ratio, true, t_bisect);
      //      cout << " after propagation " << x << endl;
      emptiness = x.is_empty();
      double volume_before_var3b;
      if (!emptiness && m_var3b_fxpt_ratio>=0.0){
	do
	  { 
	    volume_before_var3b=x.volume();
	    //	    cout << " volume before var3b "  << x.volume() << endl;
	    var3b(x, f, ctc_func);
	    //	    cout << " volume after var3b "  << x.volume() << endl;
	    emptiness = x.is_empty();
	  }
	while (!emptiness  
	       && !(stopping_condition_met(x))
	       &&  m_var3b_fxpt_ratio>0.0 
	       && !fixed_point_reached(volume_before_var3b, x.volume(), m_var3b_fxpt_ratio));
      }

      if (! emptiness)
	do
      {
        volume_before_refining = x.volume();
	cout << " volume before refining " <<  volume_before_refining << endl;
        // 1. Refining

	if(m_refining_fxpt_ratio >= 0.0)
	  if (! refining(x))
	    {
	      if (m_trace)
		cout << " end refining " <<  volume_before_refining << " after  " << x.volume() << endl; 
	      break;}
	if (m_trace) cout << " nb_slices after refining step " << x[0].nb_slices() << endl;
	// 2. Propagations up to the fixed point

	propagation(x, f, ctc_func, m_propa_fxpt_ratio, false, x[0].domain().lb());
	// 3.      
	emptiness = x.is_empty();
	double volume_before_var3b;
	if (!emptiness && m_var3b_fxpt_ratio >= 0.0){

	  do
	    { 
	      volume_before_var3b=x.volume();
	      // cout << " volume before var3b "  << x.volume() << endl;
	      var3b(x, f, ctc_func);
	      // cout << " volume after var3b "  << x.volume() << endl;
	      emptiness = x.is_empty();
	    }
	  while(!emptiness  
		&& !(stopping_condition_met(x))
		&& !fixed_point_reached(volume_before_var3b, x.volume(), m_var3b_fxpt_ratio));
	  }
	if (m_trace) cout << " volume after refining " <<  x.volume() << endl;
      }
      
      while(!emptiness
	    && !(stopping_condition_met(x))
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
            if (m_trace) cout << "Bisection... (level " << level << ")" << endl;
	    //	    if (f) bisection_guess (x,*f);  //TODO use bisection_guess
	    double t_bisection;
	    if (x.volume() < DBL_MAX){
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

	   
	    else  // infinite tube : 
	      t_bisection=one_finite_gate(x);
	    }
	    bisections++;
	    level++;
            try{
	      pair<TubeVector,TubeVector> p_x = x.bisect(t_bisection);
	  
	      //   s.push_back(make_pair(level, p_x.second));   // breadth first variant
	      //   s.push_back (make_pair(level, p_x.first));
	    
	      s.push_front(make_pair(make_pair(level,t_bisection), p_x.second));       // depth first variant
	      s.push_front(make_pair(make_pair(level,t_bisection), p_x.first));
              if (m_trace)	      cout << " t_bisection " << t_bisection << " x volume " << x.volume() << " nb_slices " << x.nb_slices()  << endl;
	      //	      cout << "last slice v0 f" << * (p_x.first[0].last_slice())  << endl;
	      //	      cout << "last slice v1 f" << * (p_x.first[1].last_slice())  << endl;
	      //	      cout << "last slice v0 s" << * (p_x.second[0].last_slice())  << endl;
	      //	      cout << "last slice v1 s" << * (p_x.second[1].last_slice())  << endl;
	    }
	    catch (Exception &)   // when the bisection time was not bisectable, change to largest gate
	      {	 
		cout << " exception " << endl;
		x.max_gate_diam(t_bisection);
		//                t_bisection=x.largest_slice()->domain().mid();
	        pair<TubeVector,TubeVector> p_x = x.bisect(t_bisection);

             if (m_trace)	      cout << " t_bisection " << t_bisection << " x volume " << x.volume() << " nb_slices " << x.nb_slices()  << endl;
	     //	      cout << "last slice v0 f" << * (p_x.first[0].slice(t_bisection))  << endl;
	     //	      cout << "last slice v1 f" << * (p_x.first[1].slice(t_bisection))  << endl;
	     //	      cout << "last slice v0 s" << * (p_x.second[0].slice(t_bisection))  << endl;
	     //	      cout << "last slice v1 s" << * (p_x.second[1].slice(t_bisection))  << endl;


	    //   s.push_back(make_pair(level, p_x.second));   // breadth first variant
	    //   s.push_back (make_pair(level, p_x.first));
	    
		s.push_front(make_pair(make_pair(level,t_bisection), p_x.second));       // depth first variant
		s.push_front(make_pair(make_pair(level,t_bisection), p_x.first));
	      }


	    //  cout << " t_bisection " << t_bisection << " x volume " << x.volume() << " nb_slices " << x.nb_slices() <<   endl;
             
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
    if (m_trace) cout << "bisections " << bisections << endl;
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
    assert(!l_tubes.empty());
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
    //    assert(x.size() == m_max_thickness.size());
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

    int n = m_max_thickness.size();
    return (std::pow(volume_after, 1./n) / std::pow(volume_before, 1./n)) >= fxpt_ratio;
  }

  void Solver::picard_contraction (TubeVector &x, tubex::Fnc& f){
    if (x.volume()>= DBL_MAX){
      //      cout << " volume before picard " << x.volume() << endl;
      CtcPicard ctc_picard;
      ctc_picard.preserve_slicing(true);
      ctc_picard.contract(f, x, FORWARD | BACKWARD);
      //   cout << " volume after picard " << x.volume() << endl;
    }
  }

  void Solver::deriv_contraction (TubeVector &x, tubex::Fnc& f){
    CtcDeriv ctc;
    //	    cout << " volume before ctc deriv " << volume_before_ctc << endl;
    ctc.set_fast_mode(true);
    ctc.contract(x, f.eval_vector(x),FORWARD | BACKWARD);
    //	   cout << "volume after ctc deriv " << x.volume() << endl;
  }

  void Solver::integration_contraction(TubeVector &x, tubex::Fnc& f, double t0, bool incremental){
    
    Ctc* ctc_dyncid;
    CtcIntegration* ctc_integration;
    if (m_contraction_mode==0)
      ctc_dyncid =  new CtcDynBasic(f);
    else if (m_contraction_mode==1) 
      ctc_dyncid =  new CtcDynCid(f);
    else if (m_contraction_mode==2){
      // ctc_dyncid =  new CtcDynCidGuess(f,0.);
      ctc_dyncid =  new CtcDynCidGuess(f);
      //      (dynamic_cast <CtcDynCidGuess*> (ctc_dyncid))->set_variant(1);
      //      (dynamic_cast <CtcDynCidGuess*> (ctc_dyncid))->set_dpolicy(2);
          }
  
    ctc_dyncid->set_fast_mode(true);

    ctc_integration = new CtcIntegration (f,ctc_dyncid);
    if(x.volume() >= DBL_MAX) ctc_integration->set_picard_mode(true);

    incremental=false ; // stronger contraction without incrementality ; comment this line for incrementality
    TubeVector v = f.eval_vector(x);
    if (incremental) // incrementality if incremental and no other function
      {
	//	ctc_integration->set_incremental_mode(incremental);
	ctc_integration->set_incremental_mode(false);
	/*
	if (x.volume() < 27.39 && x.volume() > 27.36) 
	  ctc_integration->set_incremental_mode(false);
	*/
	//	cout << " before " << x.volume() << " t0 " << t0 << endl;
	//	ctc_integration->contract(x,v,x[0].domain().lb(),FORWARD);
	//	v = f.eval_vector(x);
	ctc_integration->contract(x,v,t0,FORWARD) ;
	//	v = f.eval_vector(x);
	//	cout << " after t0 forward " <<  x.volume() << endl;
	ctc_integration->contract(x,v,x[0].domain().ub(),BACKWARD);
	//	v = f.eval_vector(x);
	//	cout << " after tf backward  " << x.volume()<<  endl;
	ctc_integration->contract(x,v,t0,BACKWARD) ;
	//	v = f.eval_vector(x);
	ctc_integration->contract(x,v,x[0].domain().lb(),FORWARD);
	//	ctc_integration->contract(x,v,t0,FORWARD) ;
	//	cout << " after t0 backward " << x.volume()<< endl;

	//	cout << " end propagation " <<  x.volume() << endl;
	//	ctc_integration->contract(x,v,t0,FORWARD) ;
	//	ctc_integration->contract(x,v,x[0].domain().ub(),BACKWARD);
	//	ctc_integration->contract(x,v,t0,BACKWARD);
      }
    else
      {
      ctc_integration->set_incremental_mode(false);
      //      cout << " x before " << x << endl;
      ctc_integration->contract(x,v,x[0].domain().lb(),FORWARD);
      //      cout << " x after forward " << x << endl;
      ctc_integration->contract(x,v,x[0].domain().ub(),BACKWARD);
      //      cout << " x after backward " << x << endl;
      }
    delete ctc_dyncid; delete ctc_integration;
  }



  void Solver::propagation(TubeVector &x, tubex::Fnc* f, void (*ctc_func)(TubeVector&, double& t0, bool incremental), float propa_fxpt_ratio, bool incremental, double t0 )
  {
    if  (propa_fxpt_ratio <0.0) return;
    double volume_before_ctc;
    bool first_iteration=true;
    do
      {
	volume_before_ctc = x.volume();
        if (ctc_func) {ctc_func(x, t0, incremental);  // Other constraints contraction
	  incremental=false;}
	if (f){                     // ODE contraction
	  
	  if (m_contraction_mode==4){   // CtcPicard + CtcDeriv
	    picard_contraction(x,*f);
	    deriv_contraction(x,*f);
	  }
    
	  else if (m_contraction_mode <=2 ){                           // CtcIntegration 
	    if (first_iteration)
	      integration_contraction(x,*f,t0,incremental);
	    else
	      integration_contraction(x,*f,t0,false);
	  }
	  //	  cout << " tube after contraction " << x << " volume after " << x.volume() << endl;
	  
	}
	first_iteration=false;
      }

    while(!(x.is_empty()) 
	  && propa_fxpt_ratio >0
	  && !fixed_point_reached(volume_before_ctc, x.volume(), propa_fxpt_ratio));
    //    cout << " tube after propagation " <<  x << " volume after " << x.volume() << endl;
  }





  void Solver::var3b(TubeVector &x, tubex::Fnc * f,void (*ctc_func) (TubeVector&,double& t0,bool incremental))
  {
    //    cout << " volume before var3b " << x.volume() << endl;
    if(m_var3b_fxpt_ratio < 0. || x.volume() >= DBL_MAX)
      return;

    int contraction_mode = m_contraction_mode;
    m_contraction_mode=4;  // var3B calls CtcDeriv as internal contractor 
    // incremental contractors using CtcIntegration are too weak and non incremental contractors as
    // CtcIntegration with CtcDyncid or CtcDyncidGuess are too costly
    //    cout << " var3b " << x << endl;
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
	     //	     cout << " before var3b prop1 " << branch_x << endl;
	     propagation(branch_x, f, ctc_func, m_var3b_propa_fxpt_ratio, true, t_bisection);
	     //	     cout << " after var3b prop1 " << endl;
	     if (branch_x.is_empty())
	       x=p_x.second;
	     else {x = p_x.second | branch_x ; break;}

	     rate= m_var3b_bisection_ratefactor*rate;

	    }
	  catch (Exception& )
	    {break;}
	}

	propagation(x,f, ctc_func, m_var3b_propa_fxpt_ratio, true, t_bisection);

	rate = 1 - m_var3b_bisection_minrate;
       
	while (rate > 1-m_var3b_bisection_maxrate ){
	  try{
	    pair<TubeVector,TubeVector> p_x = x.bisect(t_bisection,k,rate);
	     TubeVector branch_x= p_x.second;
	     //	     cout << " before var3b prop2 " << endl;
	     propagation(branch_x, f, ctc_func, m_var3b_propa_fxpt_ratio, true, t_bisection);
	     //	     cout << " after var3b prop2" << endl;
	     if (branch_x.is_empty())
	       x=p_x.first;
	     else {x = p_x.first | branch_x ; break;}
	     rate=1-m_var3b_bisection_ratefactor*(1-rate);

	  }
	  catch (Exception& )
	    {break;}
	}

	propagation(x,f , ctc_func, m_var3b_propa_fxpt_ratio, true , t_bisection);

       
      }
    m_contraction_mode=contraction_mode;
    //    cout << " volume after var3b " << x.volume() << endl;
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
