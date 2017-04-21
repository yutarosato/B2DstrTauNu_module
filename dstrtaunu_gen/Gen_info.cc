#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Gen_info( std::map<int,int> child_id_map_d  [][2], std::multiset<int> n_particle_set_d  [][2],
			    std::map<int,int> child_id_map_tau[][2], std::multiset<int> n_particle_set_tau[][2],
			    int gen_b_decay_info[][30],
			    int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
			    const int fl_message )
  {
    //std::cout << "Gen start!" << std::endl;

    // [ Check DATA/MC ]
    if( mc !=1 ){
      std::cout << "Data is Real Data. Gen Function skip. "
		<< std::endl;
      return;
    }

    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    if( fl_message==2 ) genMgr.dump();

    // [ Check event_type ]
    event_type = 0; // 1(uds), 2(charm), 3(unmixed), -3(mixed), 4(charged), 5(Y(nS))
    Gen_hepevt& Upsilon4S( genMgr( Panther_ID(1) ) );
    Gen_hepevt& B1       ( genMgr( Panther_ID(2) ) );
    Gen_hepevt& B2       ( genMgr( Panther_ID(3) ) );
    if     ( Upsilon4S.idhep() != Upsilon4S_LUND ) event_type = 1; // ROOT PARTICLE is Upsilon(4S) ? -> [continuum]
    else if( Upsilon4S.daLast() - Upsilon4S.daFirst() == 1 ){ // Upsilon(4S) decays two body ?
      if     (     B1.idhep()*B2.idhep() ==    B0_LUND*    B0_LUND ) event_type = -3; // -> [ mixed ]
      if     (     B1.idhep()*B2.idhep() ==    B0_LUND*antiB0_LUND ) event_type =  3; // -> [unmixed]
      else if(     B1.idhep()*B2.idhep() == Bplus_LUND*Bminus_LUND ) event_type =  4; // -> [charged]
    }else event_type = 5; // -> [Y(nS)]

    if( event_type==1 ){ // continuum
      int tmp_ndid1, tmp_ndid2;
      int nd = search_daughter_charm( Upsilon4S, tmp_ndid1, tmp_ndid2, fl_message ); // count # of D in virtual photon
      if( nd ) event_type = 2; // charm
    }else if( event_type==3 || event_type==-3 || event_type==4 ){
      // [ judgement of sig/tag side ]
      int semilept_id1 = 0;
      int semilept_id2 = 0;
      int semilept1 = semilept_dec( B1, semilept_id1, false );
      int semilept2 = semilept_dec( B2, semilept_id2, false );
      if( semilept1==3 ){ // B1 is signal side
	gen_b_decay_info[0][0] = 2;
	gen_b_decay_info[1][0] = 3;
      }else if( semilept2==3 ){ // B2 is signal side
	gen_b_decay_info[0][0] = 3;
	gen_b_decay_info[1][0] = 2;
      }else{ // no signal events
	gen_b_decay_info[0][0] = 2;
	gen_b_decay_info[1][0] = 3;
      }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // [ Check B decay ]
      Gen_hepevt B[2];
      for( int i=0; i<2; i++ ){
	B[i] = Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[i][0]) ) ); // B[0] : Bsig, B[1] : Btag
	gen_b_decay_info[i][1] = semilept_dec( B[i], gen_b_decay_info[i][2], false );      // Check semileptonic B decay
	gen_b_decay_info[i][3] = ( B[i].daLast()==0 ? 0 : B[i].daLast() - B[i].daFirst() + 1 ); // Check B children
	
	// [ Check D particle in B decay ]
	for( int j=0; j< gen_b_decay_info[i][3]; j++ ){
	if( genMgr( Panther_ID( B[i].daFirst() + j ) ).idhep()==Gamma_LUND ) gen_b_decay_info[i][4]++; // count # of gamma directly from B
	if( D_set.count( genMgr( Panther_ID( B[i].daFirst() + j ) ).idhep() ) ){
	  if( gen_b_decay_info[i][5]==0 ){
	    gen_b_decay_info[i][6] = genMgr( Panther_ID( B[i].daFirst() + j ) ).idhep();  // 1st-rootD->idhep()
	    gen_b_decay_info[i][7] = genMgr( Panther_ID( B[i].daFirst() + j ) ).get_ID(); // 1st-rootD->get_ID()
	  }else{
	    gen_b_decay_info[i][8] = genMgr( Panther_ID( B[i].daFirst() + j ) ).idhep();  // 2nd-rootD->idhep()
	    gen_b_decay_info[i][9] = genMgr( Panther_ID( B[i].daFirst() + j ) ).get_ID(); // 2nd-rootD->get_ID()
	  }
	  gen_b_decay_info[i][5]++; // count # of rootD
	}
	}
	
	// [ Check D particle ]
	gen_b_decay_info[i][10] = search_daughter_D( B[i], gen_b_decay_info[i][11], gen_b_decay_info[i][12], fl_message ); // count # of D in B decay
	
	// [ Flag og D/D*/D**/DD ]
	if     ( gen_b_decay_info[i][5]> 1 ) gen_b_decay_info[i][13] = 4; // two D
	else if( gen_b_decay_info[i][5]==0 ) gen_b_decay_info[i][13] = 5; //  no D
	else if( abs(gen_b_decay_info[i][6])==D0_LUND    || abs(gen_b_decay_info[i][6])==Dplus_LUND ) gen_b_decay_info[i][13] = 1; // D
	else if( abs(gen_b_decay_info[i][6])==Dstr0_LUND || abs(gen_b_decay_info[i][6])==Dstrp_LUND ) gen_b_decay_info[i][13] = 2; // D*
	else                                                                                          gen_b_decay_info[i][13] = 3; // D**
	
	if( gen_b_decay_info[i][13]==2 ){ // if D* exist, accompany particle in D* decay is checked.
	  Gen_hepevt& acc( genMgr( Panther_ID(gen_b_decay_info[i][7]) ) );
	  gen_b_decay_info[i][14] = genMgr( Panther_ID( acc.daFirst()+1 ) ).idhep();
	  gen_b_decay_info[i][15] = genMgr( Panther_ID( acc.daFirst()+1 ) ).get_ID();
	}
	
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// search Charmonium from B decay
	gen_b_decay_info[i][16] = search_daughter_CC( B[i], false );
	
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// [ Check D decay ]
	if( gen_b_decay_info[i][10]>0 ){
	  Gen_hepevt& D( genMgr( Panther_ID(gen_b_decay_info[i][11]) ) );
	  gen_b_decay_info[i][17] = semilept_dec( D, gen_b_decay_info[i][18], false );	
	  find_fin_child( D, child_id_map_d[i][0], n_particle_set_d[i][0], fl_message );
	  if( make_mode_digit( n_particle_set_d[i][0], gen_d_mode_info[i][0] ) ) std::cerr << "[WARNING] n_rest(1st-D) exists" << std::endl;
      }
	if( gen_b_decay_info[i][10]>1 ){
	  Gen_hepevt& D( genMgr( Panther_ID(gen_b_decay_info[i][12]) ) );
	  gen_b_decay_info[i][19] = semilept_dec( D, gen_b_decay_info[i][20], fl_message );	
	  find_fin_child( D, child_id_map_d[i][1], n_particle_set_d[i][1], fl_message );
	  if( make_mode_digit( n_particle_set_d[i][1], gen_d_mode_info[i][1] ) ) std::cerr << "[WARNING] n_rest(2nd-D) exists" << std::endl;
	}
	
	if( fl_message ){
	  std::cout << "------------------------------------------------------------------------------"
		  << std::endl;
	  std::cout << " # of 1stD rec-child : " << child_id_map_d[i][0].size() << " [ ";
	  for( std::map<int, int>::iterator m = child_id_map_d[i][0].begin(); m != child_id_map_d[i][0].end(); m++ ){
	    std::cout << m->second << "(" << m->first << ") ";
	  }
	  std::cout << " ] " << std::endl;
	  std::cout << " # of 1stD all-child : " <<  n_particle_set_d[i][0].size() << " [ ";
	  int count=0;
	  for( std::multiset<int>::iterator m = n_particle_set_d[i][0].begin(); m != n_particle_set_d[i][0].end(); m++ ){
	    std::cout << *m << " ";
	  }
	  std::cout << " ]" << std::endl;
	}
	
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// [ Check tau decay ]
	if( gen_b_decay_info[i][1]==3 ){
	  Gen_hepevt& tau( genMgr( Panther_ID(gen_b_decay_info[i][2]) ) );
	  find_fin_child( tau, child_id_map_tau[i][0], n_particle_set_tau[i][0], fl_message );
	  if( make_mode_digit( n_particle_set_tau[i][0], gen_tau_mode_info[i] ) )  std::cerr << "[WARNING] n_rest(tau) exists" << std::endl;
	  search_prompt_lepton( tau, gen_b_decay_info[i][21], gen_b_decay_info[i][22], false );
	}
      }
    }


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // [ Event number ]
    double exprunNo = expNo*10000+runNo;
    Gen_dist->column( "exprun",    exprunNo   );
    Gen_dist->column( "event",     evtNo      );
    Gen_dist->column( "evt_type",  event_type );

    Gen_dist->column( "nb1",       gen_b_decay_info[0][ 3] ); Gen_dist->column( "nb2",       gen_b_decay_info[1][ 3] );
    Gen_dist->column( "nb1g",      gen_b_decay_info[0][ 4] ); Gen_dist->column( "nb2g",      gen_b_decay_info[1][ 4] );
    Gen_dist->column( "semi1",     gen_b_decay_info[0][ 1] ); Gen_dist->column( "semi2",     gen_b_decay_info[1][ 1] );
    Gen_dist->column( "nrootd1",   gen_b_decay_info[0][ 5] ); Gen_dist->column( "nrootd2",   gen_b_decay_info[1][ 5] );
    Gen_dist->column( "rootdf1",   gen_b_decay_info[0][ 6] ); Gen_dist->column( "rootdf2",   gen_b_decay_info[1][ 6] ); 
    Gen_dist->column( "rootds1",   gen_b_decay_info[0][ 8] ); Gen_dist->column( "rootds2",   gen_b_decay_info[1][ 8] ); // buf fixed @ 20140308
    Gen_dist->column( "nd1",       gen_b_decay_info[0][10] ); Gen_dist->column( "nd2",       gen_b_decay_info[1][10] );
    Gen_dist->column( "gm_ddst1",  gen_b_decay_info[0][13] ); Gen_dist->column( "gm_ddst2",  gen_b_decay_info[1][13] );
    Gen_dist->column( "dst1_acc",  gen_b_decay_info[0][14] ); Gen_dist->column( "dst2_acc",  gen_b_decay_info[1][14] );
    Gen_dist->column( "cc1",       gen_b_decay_info[0][16] ); Gen_dist->column( "cc2",       gen_b_decay_info[1][16] );

    if( gen_b_decay_info[0][15] ) Gen_dist->column( "acc1_E", genMgr( Panther_ID(gen_b_decay_info[0][15]) ).E() );
    else                          Gen_dist->column( "acc1_E", -10 );
    if( gen_b_decay_info[1][15] ) Gen_dist->column( "acc2_E", genMgr( Panther_ID(gen_b_decay_info[1][15]) ).E() );
    else                          Gen_dist->column( "acc2_E", -10 );

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Gen_dist->column( "pi0_E",   -10 ); // initialize
    Gen_dist->column( "pi0_cos", -10 ); // initialize
    Gen_dist->column( "pi0_ghE", -10 ); // initialize
    Gen_dist->column( "pi0_glE", -10 ); // initialize
    if( gen_b_decay_info[0][11] ){
      int nchild = genMgr( Panther_ID(gen_b_decay_info[0][11]) ).daLast()==0 ? 0 : genMgr( Panther_ID(gen_b_decay_info[0][11]) ).daLast() - genMgr( Panther_ID(gen_b_decay_info[0][11]) ).daFirst() + 1;
      for( int i=0; i<nchild; i++ ){
	Gen_hepevt& pi0_cand( genMgr( Panther_ID(genMgr(Panther_ID(gen_b_decay_info[0][11])).daFirst() + i ) ) );
	if( pi0_cand.idhep()==PI0_LUND ){
	  Gen_dist->column( "pi0_E",  pi0_cand.E() );
	  int pi0_nchild = pi0_cand.daLast()==0 ? 0 : pi0_cand.daLast() - pi0_cand.daFirst() + 1;
	  if( pi0_nchild==2 ){ // for pi0 -> two gamma
	    Gen_hepevt& acc_pi0_gam1( genMgr( Panther_ID(pi0_cand.daFirst()+0    ) ) ); 
	    Gen_hepevt& acc_pi0_gam2( genMgr( Panther_ID(pi0_cand.daFirst()+1    ) ) );
	    Hep3Vector par1_mom(acc_pi0_gam1.PX(), acc_pi0_gam1.PY(), acc_pi0_gam1.PZ());
	    Hep3Vector par2_mom(acc_pi0_gam2.PX(), acc_pi0_gam2.PY(), acc_pi0_gam2.PZ());
	    double cos = (par1_mom * par2_mom)/(par1_mom.mag()*par2_mom.mag());
	    Gen_dist->column( "pi0_cos", cos );
	    if( acc_pi0_gam1.E() > acc_pi0_gam2.E() ){
	      Gen_dist->column( "pi0_ghE", acc_pi0_gam1.E() );
	      Gen_dist->column( "pi0_glE", acc_pi0_gam2.E() );
	    }else{
	      Gen_dist->column( "pi0_glE", acc_pi0_gam1.E() );
	      Gen_dist->column( "pi0_ghE", acc_pi0_gam2.E() );
	    }
	  }
	}
      }
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Gen_dist->column( "acc1_ghE", -10 ); // initialize
    Gen_dist->column( "acc1_glE", -10 ); // initialize
    Gen_dist->column( "acc2_ghE", -10 ); // initialize
    Gen_dist->column( "acc2_glE", -10 ); // initialize
    Gen_dist->column( "acc1_cos", -10 ); // initialize
    Gen_dist->column( "acc2_cos", -10 ); // initialize

    if( gen_b_decay_info[0][14]==PI0_LUND ){
      int nchild = genMgr( Panther_ID(gen_b_decay_info[0][15]) ).daLast()==0 ? 0 : genMgr( Panther_ID(gen_b_decay_info[0][15]) ).daLast() - genMgr( Panther_ID(gen_b_decay_info[0][15]) ).daFirst() + 1;
      if( nchild==2 ){ // for pi0 -> two gamma
	Gen_hepevt& acc_pi0     ( genMgr( Panther_ID(gen_b_decay_info[0][15]) ) );
	Gen_hepevt& acc_pi0_gam1( genMgr( Panther_ID(acc_pi0.daFirst()+0    ) ) ); 
	Gen_hepevt& acc_pi0_gam2( genMgr( Panther_ID(acc_pi0.daFirst()+1    ) ) );
	Hep3Vector par1_mom(acc_pi0_gam1.PX(), acc_pi0_gam1.PY(), acc_pi0_gam1.PZ());
	Hep3Vector par2_mom(acc_pi0_gam2.PX(), acc_pi0_gam2.PY(), acc_pi0_gam2.PZ());
	double cos = (par1_mom * par2_mom)/(par1_mom.mag()*par2_mom.mag());
	Gen_dist->column( "acc1_cos", cos );
	if( acc_pi0_gam1.E() > acc_pi0_gam2.E() ){
	  Gen_dist->column( "acc1_ghE", acc_pi0_gam1.E() );
	  Gen_dist->column( "acc1_glE", acc_pi0_gam2.E() );
	}else{
	  Gen_dist->column( "acc1_glE", acc_pi0_gam1.E() );
	  Gen_dist->column( "acc1_ghE", acc_pi0_gam2.E() );
	}
      }
    }
    
    if( gen_b_decay_info[1][14]==PI0_LUND ){
      int nchild = genMgr( Panther_ID(gen_b_decay_info[1][15]) ).daLast()==0 ? 0 : genMgr( Panther_ID(gen_b_decay_info[1][15]) ).daLast() - genMgr( Panther_ID(gen_b_decay_info[1][15]) ).daFirst() + 1;
      if( nchild==2 ){ // for pi0 -> two gamma
	Gen_hepevt& acc_pi0     ( genMgr( Panther_ID(gen_b_decay_info[1][15]) ) );
	Gen_hepevt& acc_pi0_gam1( genMgr( Panther_ID(acc_pi0.daFirst()+0    ) ) ); 
	Gen_hepevt& acc_pi0_gam2( genMgr( Panther_ID(acc_pi0.daFirst()+1    ) ) );
	Hep3Vector par1_mom(acc_pi0_gam1.PX(), acc_pi0_gam1.PY(), acc_pi0_gam1.PZ());
	Hep3Vector par2_mom(acc_pi0_gam2.PX(), acc_pi0_gam2.PY(), acc_pi0_gam2.PZ());
	double cos = (par1_mom * par2_mom)/(par1_mom.mag()*par2_mom.mag());
	Gen_dist->column( "acc2_cos", cos );
	if( acc_pi0_gam1.E() > acc_pi0_gam2.E() ){
	  Gen_dist->column( "acc2_ghE", acc_pi0_gam1.E() );
	  Gen_dist->column( "acc2_glE", acc_pi0_gam2.E() );
	}else{
	  Gen_dist->column( "acc2_glE", acc_pi0_gam1.E() );
	  Gen_dist->column( "acc2_ghE", acc_pi0_gam2.E() );
	}
      }
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Gen_dist->column( "gm_df1hr",  gen_d_mode_info [0][0][0]  ); Gen_dist->column( "gm_df2hr",  gen_d_mode_info [1][0][0] );
    Gen_dist->column( "gm_ds1hr",  gen_d_mode_info [0][1][0]  ); Gen_dist->column( "gm_ds2hr",  gen_d_mode_info [1][1][0] );
    Gen_dist->column( "gm_df1hnr", gen_d_mode_info [0][0][1]  ); Gen_dist->column( "gm_df2hnr", gen_d_mode_info [1][0][1] );
    Gen_dist->column( "gm_ds1hnr", gen_d_mode_info [0][1][1]  ); Gen_dist->column( "gm_ds2hnr", gen_d_mode_info [1][1][1] );
    Gen_dist->column( "gm_df1lep", gen_d_mode_info [0][0][2]  ); Gen_dist->column( "gm_df2lep", gen_d_mode_info [1][0][2] );
    Gen_dist->column( "gm_ds1lep", gen_d_mode_info [0][1][2]  ); Gen_dist->column( "gm_ds2lep", gen_d_mode_info [1][1][2] );
    Gen_dist->column( "gm_df1nu",  gen_d_mode_info [0][0][3]  ); Gen_dist->column( "gm_df2nu",  gen_d_mode_info [1][0][3] );
    Gen_dist->column( "gm_ds1nu",  gen_d_mode_info [0][1][3]  ); Gen_dist->column( "gm_ds2nu",  gen_d_mode_info [1][1][3] );
    Gen_dist->column( "gm_df1gam", gen_d_mode_info [0][0][4]  ); Gen_dist->column( "gm_df2gam", gen_d_mode_info [1][0][4] );
    Gen_dist->column( "gm_ds1gam", gen_d_mode_info [0][1][4]  ); Gen_dist->column( "gm_ds2gam", gen_d_mode_info [1][1][4] );

    Gen_dist->column( "gm_t1hr",  gen_tau_mode_info [0][0]  ); Gen_dist->column( "gm_t2hr",  gen_tau_mode_info [1][0] );
    Gen_dist->column( "gm_t1hnr", gen_tau_mode_info [0][1]  ); Gen_dist->column( "gm_t2hnr", gen_tau_mode_info [1][1] );
    Gen_dist->column( "gm_t1lep", gen_tau_mode_info [0][2]  ); Gen_dist->column( "gm_t2lep", gen_tau_mode_info [1][2] );
    Gen_dist->column( "gm_t1nu",  gen_tau_mode_info [0][3]  ); Gen_dist->column( "gm_t2nu",  gen_tau_mode_info [1][3] );
    Gen_dist->column( "gm_t1gam", gen_tau_mode_info [0][4]  ); Gen_dist->column( "gm_t2gam", gen_tau_mode_info [1][4] );
    Gen_dist->column( "gm_t1pro", gen_tau_mode_info [0][5]  ); Gen_dist->column( "gm_t2pro", gen_tau_mode_info [1][5] );

    Gen_dist->column( "gm_dd",   make_flag_dd  (gen_b_decay_info[0][13], gen_b_decay_info[1][13]) );
    Gen_dist->column( "gm_semi", make_flag_semi(gen_b_decay_info[0][ 1], gen_b_decay_info[1][ 1]) );

    Gen_dist->dumpData();
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
