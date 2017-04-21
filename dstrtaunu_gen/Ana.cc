#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Ana( std::vector<Particle>& dtau_list, 
		       int gen_b_decay_info[][30], 
		       int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
		       std::vector<Particle>& trk_list, std::vector<Particle>& gam_list, std::vector<Particle>& pi0_list, std::vector<Particle>& ks_list, 
		       const bool fl_message )
  {

    if( fl_message ) std::cout << "AnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAna" << std::endl;
    for( std::vector<Particle>::iterator dtau1 = dtau_list.begin(); dtau1 != dtau_list.end(); dtau1++ ){
      for( std::vector<Particle>::iterator dtau2 = dtau1+1; dtau2 != dtau_list.end(); dtau2++ ){

	// [ number ]
	double exprunNo = expNo*10000+runNo;
	Rec_dist->column( "exprun",    exprunNo   );
	Rec_dist->column( "event",     evtNo      );
	Rec_dist->column( "evt_type",  event_type );


	// [ check (DD/D*D/D*D*) ]
	// (if D*D, recD1 is D* and recD2 is D)
	int fl_dstr1  = 0; // 1(D*), 0(D)
	int fl_dstr2  = 0;
	int fl_switch = 0; // 1(DD*), other(0)
	if     ( abs(dtau1->child(0).lund())==Dstrp_LUND || abs(dtau1->child(0).lund())==Dstr0_LUND ) fl_dstr1 = 1;
	else if( abs(dtau1->child(0).lund())==Dplus_LUND || abs(dtau1->child(0).lund())==D0_LUND    ) fl_dstr1 = 0;
	else std::cerr << "[ABORT] Wrong D(*) LUND" << std::endl, abort();
	if     ( abs(dtau2->child(0).lund())==Dstrp_LUND || abs(dtau2->child(0).lund())==Dstr0_LUND ) fl_dstr2 = 1;
	else if( abs(dtau2->child(0).lund())==Dplus_LUND || abs(dtau2->child(0).lund())==D0_LUND    ) fl_dstr2 = 0;
	else std::cerr << "[ABORT] Wrong D(*) LUND" << std::endl, abort();
	if( fl_dstr1==0 && fl_dstr2==1 ) fl_switch = 1;
	int rm_dd = 0;
	if     (  fl_dstr1==1 && fl_dstr2==1  ) rm_dd = 3; // D*D*
	else if( (fl_dstr1==1 && fl_dstr2==0) ||
		 (fl_dstr1==0 && fl_dstr2==1) ) rm_dd = 2; // D*D
	else if(  fl_dstr1==0 && fl_dstr2==0  ) rm_dd = 1; // DD
	Rec_dist->column( "rm_dd",  rm_dd );
	
	   
	Particle& Dtau1 = (fl_switch ? *dtau2 : *dtau1);
	Particle& Dstr1 = Dtau1.child(0); // if Dtaunu mode, "Dstr1" is D (not D*)
	UserInfo& info_Dtau1 = dynamic_cast<UserInfo&>( Dtau1.userInfo() );
	UserInfo& info_Dstr1 = dynamic_cast<UserInfo&>( Dstr1.userInfo() );
	
	Particle& D1    = ( info_Dstr1.dm()==-1 ? Dtau1.child(0) : Dstr1.child(0) );
	Particle& lep1  = Dtau1.child(1);
	UserInfo& info_D1    = dynamic_cast<UserInfo&>( D1.userInfo()    );
	UserInfo& info_lep1  = dynamic_cast<UserInfo&>( lep1.userInfo()  );

	Particle& Dtau2 = (fl_switch ? *dtau1 : *dtau2);
	Particle& Dstr2 = Dtau2.child(0); // if Dtaunu mode, "Dstr1" is D (not D*)
	UserInfo& info_Dtau2 = dynamic_cast<UserInfo&>( Dtau2.userInfo() );
	UserInfo& info_Dstr2 = dynamic_cast<UserInfo&>( Dstr2.userInfo() );
	
	Particle& D2    = ( info_Dstr2.dm()==-1 ? Dtau2.child(0) : Dstr2.child(0) );
	Particle& lep2  = Dtau2.child(1);
	UserInfo& info_D2    = dynamic_cast<UserInfo&>( D2.userInfo()    );
	UserInfo& info_lep2  = dynamic_cast<UserInfo&>( lep2.userInfo()  );


	if( fl_message ) std::cout << "Dtau_CNTID(1st) = " << std::setw(3) << std::right << info_Dtau1.cntid() << ", "
				   << "Dtau_CNTID(2nd) = " << std::setw(3) << std::right << info_Dtau2.cntid();

	// charge check
	if( Dtau1.charge()+Dtau2.charge() != 0 ){
	  if( fl_message ) std::cout << " -> not allowed charge : " << Dtau1.charge() << " + " << Dtau2.charge() << std::endl;
	  continue;
	}
	Rec_dist->column( "rm_bb",   (Dtau1.charge() ? 4 : 3) ); // 3(B0B0bar) 4(B+B-)
	
	// Check duplication
	if( !check_dupli_daughter(dtau1,dtau2) ){
	  if( fl_message ) std::cout << " -> duplication " << std::endl;
	  continue;
	}

	// flavor check
	if( info_Dtau1.flavor() && info_Dtau2.flavor() && info_Dtau1.flavor()*info_Dtau2.flavor()==1 ){
	  if( fl_message ) std::cout << " -> not allowed flavor : " << info_Dtau1.flavor() << ", " << info_Dtau2.flavor() << std::endl;
	  continue;
	}


	int rm_bfl = 0; // flavor is determined
	if     ( !info_Dtau1.flavor() && !info_Dtau1.flavor() ) rm_bfl = 2; // both B-flavor is not determined
	else if( !info_Dtau1.flavor() || !info_Dtau1.flavor() ) rm_bfl = 1; // only one B-flavor is determined
	Rec_dist->column( "rm_bfl", rm_bfl );
	
	// [rec-mode]
	Rec_dist->column( "eecl",   calEclEnergyWithMatch2GoodGamma(Dtau1, Dtau2          ) );
	Rec_dist->column( "remtrk", cnt_remain_trk                 (Dtau1, Dtau2          ) );
	Rec_dist->column( "rempi0", cnt_remain_pi0                 (Dtau1, Dtau2, pi0_list) );
	Rec_dist->column( "remks",  cnt_remain_ks                  (Dtau1, Dtau2, ks_list ) );
	Rec_dist->column( "remgam", cnt_remain_gamma               (Dtau1, Dtau2, gam_list) );
	cal_Mmiss_Evis( trk_list, gam_list );

	// [rec-mode]	
	Rec_dist->column( "rm_l1",   info_Dtau1.rec_mode() );
	Rec_dist->column( "rm_d1",   info_D1.rec_mode()    );
	Rec_dist->column( "cosdl1",  info_Dtau1.cos()      );
	Rec_dist->column( "d1_m",    info_D1.m_org()       );

	Rec_dist->column( "rm_l2",  info_Dtau2.rec_mode() );
	Rec_dist->column( "rm_d2",  info_D2.rec_mode()    );
	Rec_dist->column( "cosdl2", info_Dtau2.cos()      );
	Rec_dist->column( "d2_m",   info_D2.m_org()       );

	Rec_dist->column( "d1self", info_D1.self() );
	Rec_dist->column( "d2self", info_D2.self() );

	// [momentum]
	Rec_dist->column( "l1pcm",  info_lep1.Vcm().vect().mag() );
	Rec_dist->column( "l2pcm",  info_lep2.Vcm().vect().mag() );
	Rec_dist->column( "d1pcm", info_Dstr1.Vcm().vect().mag() );
	Rec_dist->column( "d2pcm", info_Dstr2.Vcm().vect().mag() );

	Rec_dist->column( "l1p",   lep1.p().vect().mag() );
	Rec_dist->column( "l2p",   lep2.p().vect().mag() );
	Rec_dist->column( "d1p",  Dstr1.p().vect().mag() );
	Rec_dist->column( "d2p",  Dstr2.p().vect().mag() );
	Rec_dist->column( "l1pt",  lep1.p().perp()       );
	Rec_dist->column( "l2pt",  lep2.p().perp()       );


	if     ( info_Dtau1.rec_mode() == 10 ) Rec_dist->column( "l1pid",    info_lep1.eidProb()        );
	else if( info_Dtau1.rec_mode() ==  1 ) Rec_dist->column( "l1pid",    info_lep1.muonLikelihood() );
	if     ( info_Dtau2.rec_mode() == 10 ) Rec_dist->column( "l2pid",    info_lep2.eidProb()        );
	else if( info_Dtau2.rec_mode() ==  1 ) Rec_dist->column( "l2pid",    info_lep2.muonLikelihood() );

	
	if( info_Dstr1.dm()==-1 ){ // B -> D lnu
	  Rec_dist->column( "dm1",       -1 );
	  Rec_dist->column( "dst1_m",    -1 );
	  Rec_dist->column( "dst1self",  -1 );
	  Rec_dist->column( "rm_dst1",    0 );
	}else{ // B -> D* lnu
	  Rec_dist->column( "dm1",       info_Dstr1.dm()   );
	  Rec_dist->column( "dst1_m",    Dstr1.mass()      );
	  Rec_dist->column( "dst1self",  info_Dstr1.self() );
	  Rec_dist->column( "rm_dst1",   info_Dstr1.rec_mode() );
	}

	if( info_Dstr2.dm()==-1 ){ // B -> D lnu
	  Rec_dist->column( "dm2",       -1 );
	  Rec_dist->column( "dst2_m",    -1 );
	  Rec_dist->column( "dst2self",  -1 );
	  Rec_dist->column( "rm_dst2",    0 );

	}else{ // B -> D* lnu
	  Rec_dist->column( "dm2",      info_Dstr2.dm()   );
	  Rec_dist->column( "dst2_m",   Dstr2.mass()      );
	  Rec_dist->column( "dst2self", info_Dstr2.self() );
	  Rec_dist->column( "rm_dst2",  info_Dstr2.rec_mode() );
	}

	if( info_D1.self()==1 ){
	  int d_org = which_B( D1.genHepevt(),   gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  int l_org = which_B( lep1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  Rec_dist->column( "rmdmo1", D1.genHepevt().mother().idhep() );
	  if( d_org==l_org ){
	    if( (gen_b_decay_info[d_org-1][1]==1 || gen_b_decay_info[d_org-1][1]==2 ) && lep1.genHepevt().get_ID()==gen_b_decay_info[d_org-1][ 2] ) Rec_dist->column( "l1self", 2 ); // tag-side lepton
	    else if(                                gen_b_decay_info[d_org-1][1]==3   && lep1.genHepevt().get_ID()==gen_b_decay_info[d_org-1][22] ) Rec_dist->column( "l1self", 1 ); // sig-side lepton
	    else Rec_dist->column( "l1self", -3 ); // random lepton from self-D side
	  }else{
	    Rec_dist->column( "l1self", -2 ); // random lepton from opposite side to self-D
	  }
	}else{ // fake-D
	  Rec_dist->column( "l1self", -1 );
	  Rec_dist->column( "rmdmo1",  0 );
	}

	if( info_D2.self()==1 ){
	  int d_org = which_B( D2.genHepevt(),   gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  int l_org = which_B( lep2.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  Rec_dist->column( "rmdmo2", D2.genHepevt().mother().idhep() );
	  if( d_org==l_org ){
	    if( (gen_b_decay_info[d_org-1][1]==1 || gen_b_decay_info[d_org-1][1]==2 ) && lep2.genHepevt().get_ID()==gen_b_decay_info[d_org-1][ 2] ) Rec_dist->column( "l2self", 2 );
	    else if(                                gen_b_decay_info[d_org-1][1]==3   && lep2.genHepevt().get_ID()==gen_b_decay_info[d_org-1][22] ) Rec_dist->column( "l2self", 1 );
	    else Rec_dist->column( "l2self", -3 );
	  }else{
	    Rec_dist->column( "l2self", -2 );
	  }

	}else{
	  Rec_dist->column( "l2self", -1 );
	  Rec_dist->column( "rmdmo2",  0 );
	}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int recfB_id = 0;
	int recsB_id = 0;
	int recgen   = 0;
	if( info_D1.self()==1 ){
	  int d_org = which_B( D1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  if( d_org==1 ){
	    recfB_id = 0;
	    recsB_id = 1;
	    recgen   = 1;
	  }else if( d_org==2 ){
	    recfB_id =  1;
	    recsB_id =  0;
	    recgen   = -1;
	  }
	}else if( info_D2.self()==1 ){
	  int d_org = which_B( D2.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  if( d_org==1 ){
	    recfB_id =  1;
	    recsB_id =  0;
	    recgen   = -1;
	  }else if( d_org==2 ){
	    recfB_id = 1;
	    recsB_id = 0;
	    recgen   = 1;
	  }
	}else{
	  recfB_id = 0;
	  recsB_id = 1;
	  recgen   = 0;
	}


	// calculation B directon
	Hep3Vector B_plus;
	Hep3Vector B_minus;
	double rm_bdir = calBdirection( Dtau1, Dtau2, B_plus, B_minus );
	Rec_dist->column( "rm_bdir",  rm_bdir );

	// [gen]
	Rec_dist->column( "recgen",  recgen );
	Rec_dist->column( "nb1",       gen_b_decay_info[0][ 3] ); Rec_dist->column( "nb2",       gen_b_decay_info[1][ 3] );
	Rec_dist->column( "semi1",     gen_b_decay_info[0][ 1] ); Rec_dist->column( "semi2",     gen_b_decay_info[1][ 1] );
	Rec_dist->column( "nrootd1",   gen_b_decay_info[0][ 5] ); Rec_dist->column( "nrootd2",   gen_b_decay_info[1][ 5] );
	Rec_dist->column( "rootdf1",   gen_b_decay_info[0][ 6] ); Rec_dist->column( "rootdf2",   gen_b_decay_info[1][ 6] ); 
	Rec_dist->column( "rootds1",   gen_b_decay_info[0][ 7] ); Rec_dist->column( "rootds2",   gen_b_decay_info[1][ 7] );
	Rec_dist->column( "nd1",       gen_b_decay_info[0][10] ); Rec_dist->column( "nd2",       gen_b_decay_info[1][10] );
	Rec_dist->column( "gm_ddst1",  gen_b_decay_info[0][13] ); Rec_dist->column( "gm_ddst2",  gen_b_decay_info[1][13] );
	Rec_dist->column( "dst1_acc",  gen_b_decay_info[0][14] ); Rec_dist->column( "dst2_acc",  gen_b_decay_info[1][14] );
	Rec_dist->column( "cc1",       gen_b_decay_info[0][16] ); Rec_dist->column( "cc2",       gen_b_decay_info[1][16] );
	// [gen-MODE]	
	Rec_dist->column( "gm_df1hr",  gen_d_mode_info [0][0][0]  ); Rec_dist->column( "gm_df2hr",  gen_d_mode_info [1][0][0] );
	Rec_dist->column( "gm_ds1hr",  gen_d_mode_info [0][1][0]  ); Rec_dist->column( "gm_ds2hr",  gen_d_mode_info [1][1][0] );
	Rec_dist->column( "gm_df1hnr", gen_d_mode_info [0][0][1]  ); Rec_dist->column( "gm_df2hnr", gen_d_mode_info [1][0][1] );
	Rec_dist->column( "gm_ds1hnr", gen_d_mode_info [0][1][1]  ); Rec_dist->column( "gm_ds2hnr", gen_d_mode_info [1][1][1] );
	Rec_dist->column( "gm_df1lep", gen_d_mode_info [0][0][2]  ); Rec_dist->column( "gm_df2lep", gen_d_mode_info [1][0][2] );
	Rec_dist->column( "gm_ds1lep", gen_d_mode_info [0][1][2]  ); Rec_dist->column( "gm_ds2lep", gen_d_mode_info [1][1][2] );
	Rec_dist->column( "gm_df1nu",  gen_d_mode_info [0][0][3]  ); Rec_dist->column( "gm_df2nu",  gen_d_mode_info [1][0][3] );
	Rec_dist->column( "gm_ds1nu",  gen_d_mode_info [0][1][3]  ); Rec_dist->column( "gm_ds2nu",  gen_d_mode_info [1][1][3] );
	Rec_dist->column( "gm_df1gam", gen_d_mode_info [0][0][4]  ); Rec_dist->column( "gm_df2gam", gen_d_mode_info [1][0][4] );
	Rec_dist->column( "gm_ds1gam", gen_d_mode_info [0][1][4]  ); Rec_dist->column( "gm_ds2gam", gen_d_mode_info [1][1][4] );
	
	Rec_dist->column( "gm_t1hr",  gen_tau_mode_info [0][0]  ); Rec_dist->column( "gm_t2hr",  gen_tau_mode_info [1][0] );
	Rec_dist->column( "gm_t1hnr", gen_tau_mode_info [0][1]  ); Rec_dist->column( "gm_t2hnr", gen_tau_mode_info [1][1] );
	Rec_dist->column( "gm_t1lep", gen_tau_mode_info [0][2]  ); Rec_dist->column( "gm_t2lep", gen_tau_mode_info [1][2] );
	Rec_dist->column( "gm_t1nu",  gen_tau_mode_info [0][3]  ); Rec_dist->column( "gm_t2nu",  gen_tau_mode_info [1][3] );
	Rec_dist->column( "gm_t1gam", gen_tau_mode_info [0][4]  ); Rec_dist->column( "gm_t2gam", gen_tau_mode_info [1][4] );
	Rec_dist->column( "gm_t1pro", gen_tau_mode_info [0][5]  ); Rec_dist->column( "gm_t2pro", gen_tau_mode_info [1][5] );

	Rec_dist->column( "gm_dd",   make_flag_dd  (gen_b_decay_info[0][13], gen_b_decay_info[1][13]) );
	Rec_dist->column( "gm_semi", make_flag_semi(gen_b_decay_info[0][ 1], gen_b_decay_info[1][ 1]) );

	Rec_dist->dumpData();
	if( fl_message ) std::cout << " -> dump !" << std::endl;
      }
    }
    if( fl_message ) std::cout << "AnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAna" << std::endl;

  }

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
