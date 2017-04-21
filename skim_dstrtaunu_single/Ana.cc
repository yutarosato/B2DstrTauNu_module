#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  int DSTRTAUNU::Ana( std::vector<Particle>& dtau_list, 
		      int gen_b_decay_info[][40], 
		      int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
		      std::vector<Particle>& trk_list, std::vector<Particle>& gam_list, std::vector<Particle>& pi0_list, std::vector<Particle>& ks_list, 
		      std::vector<Particle>& pi0_list_test0, std::vector<Particle>& pi0_list_test1, std::vector<Particle>& pi0_list_test2, std::vector<Particle>& pi0_list_test3,
		      std::vector<Particle>& pi0_list_test4, std::vector<Particle>& pi0_list_test5, std::vector<Particle>& pi0_list_test6, std::vector<Particle>& pi0_list_test7, std::vector<Particle>& pi0_list_test8,
		      BelleTuple* dist, const bool fl_message )
  {
    int cnt_dump = 0;
    if( fl_message ) std::cout << "AnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAna" << std::endl;
    for( std::vector<Particle>::iterator dtau1 = dtau_list.begin(); dtau1 != dtau_list.end(); dtau1++ ){
      for( std::vector<Particle>::iterator dtau2 = dtau1+1; dtau2 != dtau_list.end(); dtau2++ ){


	// ++++++++++++++++++++++++++++++++++++++++ [ number ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	double exprunNo = expNo*10000+runNo;
	if( runNo==0 ) exprunNo = expNo*10000+runNoforsigMC; // for signal MC

	dist->column( "exprun",    exprunNo   );
	dist->column( "event",     evtNo      );
	dist->column( "evt_type",  event_type );

	// +++++++++++++++++++++++++++++++++ [ make Particl object ] ++++++++++++++++++++++++++++++++++++++++++++++++++
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


	// ++++++++++++++++++++++++++++++++++++++++ [ check ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// charge check
	if( Dtau1.charge()+Dtau2.charge() != 0 ){
	  if( fl_message ) std::cout << " -> not allowed charge : " << Dtau1.charge() << " + " << Dtau2.charge() << std::endl;
	  // continue; // removed @20140919
	}

	// Check duplication
	if( !check_dupli_daughter(dtau1,dtau2) ){
	  if( fl_message ) std::cout << " -> duplication " << std::endl;
	  continue;
	}

	// flavor check
	if( info_Dtau1.flavor() && info_Dtau2.flavor() && info_Dtau1.flavor()*info_Dtau2.flavor()==1 ){
	  if( fl_message ) std::cout << " -> not allowed flavor : " << info_Dtau1.flavor() << ", " << info_Dtau2.flavor() << std::endl;
	  //continue; // removed @20140919
	}

	// ++++++++++++++++++++++++++++++++++++ [ rec-mode(BB) ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++

	int rm_bb = (Dtau1.charge() ? 4 : 3); // 3(B0B0bar) 4(B+B-)
	dist->column( "rm_bb", rm_bb );
	dist->column( "rm_dd", rm_dd ); // 1(DD), 2(D*D), 3(D*D*)

	int rm_bfl = 0; // flavor is determined
	if     ( !info_Dtau1.flavor() && !info_Dtau2.flavor() ) rm_bfl = 2; // both B-flavor is not determined
	else if( !info_Dtau1.flavor() || !info_Dtau2.flavor() ) rm_bfl = 1; // only one B-flavor is determined
	dist->column( "rm_bfl", rm_bfl );

	// ++++++++++++++++++++++++++++++++++++ [ rec-mode(B) ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++
	dist->column( "rm_l1",       info_Dtau1.rec_mode() );
	dist->column( "rm_d1",       info_D1.rec_mode()    );
	dist->column( "rm_d1lund",   D1.lund()             );
	dist->column( "rm_dst1lund", info_Dstr1.dm()==-1 ? 0 : Dstr1.lund() );
	dist->column( "dfl1",        info_Dtau1.flavor()   );
	dist->column( "dtchg1",      Dtau1.charge()        );

	dist->column( "rm_l2",       info_Dtau2.rec_mode() );
	dist->column( "rm_d2",       info_D2.rec_mode()    );
	dist->column( "rm_d2lund",   D2.lund()             );
	dist->column( "rm_dst2lund", info_Dstr2.dm()==-1 ? 0 : Dstr2.lund() );
	dist->column( "dfl2",        info_Dtau2.flavor() );
	dist->column( "dtchg2",      Dtau2.charge() );

	if( info_Dstr1.dm()==-1 ) dist->column( "rm_dst1",   0                     ); // B -> D  lnu
	else                      dist->column( "rm_dst1",   info_Dstr1.rec_mode() ); // B -> D* lnu
	if( info_Dstr2.dm()==-1 ) dist->column( "rm_dst2",   0                     ); // B -> D  lnu
	else                      dist->column( "rm_dst2",   info_Dstr2.rec_mode() ); // B -> D* lnu

	// ++++++++++++++++++++++++++++++++++ [ event variable ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double eecl_detail [39]={0};
	int    neecl_detail[39]={0};

	dist->column( "eecl",    calEclEnergyWithMatch2GoodGamma(Dtau1, Dtau2, !true, eecl_detail, neecl_detail, NULL, gen_b_decay_info[0][0], gen_b_decay_info[1][0] ) );
	dist->column( "eeclth1", calEclEnergyWithMatch2GoodGamma(Dtau1, Dtau2, !true,        NULL,         NULL, NULL, gen_b_decay_info[0][0], gen_b_decay_info[1][0], 6.0, 0.94, 0.10, 0.05, 0.15 ) ); // same wit nominal
	dist->column( "eeclth2", calEclEnergyWithMatch2GoodGamma(Dtau1, Dtau2, !true,        NULL,         NULL, NULL, gen_b_decay_info[0][0], gen_b_decay_info[1][0], 6.0, 0.94, 0.05, 0.05, 0.05 ) );
	dist->column( "eeclth3", calEclEnergyWithMatch2GoodGamma(Dtau1, Dtau2, !true,        NULL,         NULL, NULL, gen_b_decay_info[0][0], gen_b_decay_info[1][0], 6.0, 0.94, 0.10, 0.05, 0.10 ) );

	dist->column( "neecl",   neecl_detail[0]+neecl_detail[1]+neecl_detail[2] );
	
	dist->column( "eecl_b1", eecl_detail[ 0] ); dist->column( "neecl_b1", neecl_detail[ 0] );
	dist->column( "eecl_b2", eecl_detail[ 1] ); dist->column( "neecl_b2", neecl_detail[ 1] );
	dist->column( "eecl_bg", eecl_detail[ 2] ); dist->column( "neecl_bg", neecl_detail[ 2] );
	dist->column( "eecl_ef", eecl_detail[ 3] ); dist->column( "neecl_ef", neecl_detail[ 3] );
	dist->column( "eecl_b",  eecl_detail[ 4] ); dist->column( "neecl_b",  neecl_detail[ 4] );
	dist->column( "eecl_eb", eecl_detail[ 5] ); dist->column( "neecl_eb", neecl_detail[ 5] );

	dist->column( "eecl_neu", eecl_detail[ 8]+eecl_detail[ 9]+eecl_detail[10]+eecl_detail[11]+eecl_detail[12]+eecl_detail[13] );
	dist->column( "eecl_pi0", eecl_detail[14]+eecl_detail[15]+eecl_detail[16]+eecl_detail[17]+eecl_detail[18]+eecl_detail[19] );
	dist->column( "eecl_trk", eecl_detail[24]+eecl_detail[26]+eecl_detail[28]+eecl_detail[30]+eecl_detail[32]+eecl_detail[34] );
	dist->column( "eecl_had", eecl_detail[25]+eecl_detail[27]+eecl_detail[29]+eecl_detail[31]+eecl_detail[33]+eecl_detail[35] );

	dist->column( "neecl_neu", neecl_detail[ 8]+neecl_detail[ 9]+neecl_detail[10]+neecl_detail[11]+neecl_detail[12]+neecl_detail[13] );
	dist->column( "neecl_pi0", neecl_detail[14]+neecl_detail[15]+neecl_detail[16]+neecl_detail[17]+neecl_detail[18]+neecl_detail[19] );
	dist->column( "neecl_trk", neecl_detail[24]+neecl_detail[26]+neecl_detail[28]+neecl_detail[30]+neecl_detail[32]+neecl_detail[34] );
	dist->column( "neecl_had", neecl_detail[25]+neecl_detail[27]+neecl_detail[29]+neecl_detail[31]+neecl_detail[33]+neecl_detail[35] );

	dist->column( "eecl_0",  eecl_detail[ 6] ); dist->column( "neecl_0",  neecl_detail[ 6] ); // other?
	dist->column( "eecl_1",  eecl_detail[ 7] ); dist->column( "neecl_1",  neecl_detail[ 7] ); // beam B.G.
	dist->column( "eecl_2",  eecl_detail[ 8] ); dist->column( "neecl_2",  neecl_detail[ 8] ); // Lambda
	dist->column( "eecl_3",  eecl_detail[ 9] ); dist->column( "neecl_3",  neecl_detail[ 9] ); // Ks->pi+pi-
	dist->column( "eecl_4",  eecl_detail[10] ); dist->column( "neecl_4",  neecl_detail[10] ); // Ks->pi0pi0
	dist->column( "eecl_5",  eecl_detail[11] ); dist->column( "neecl_5",  neecl_detail[11] ); // KL
	dist->column( "eecl_6",  eecl_detail[12] ); dist->column( "neecl_6",  neecl_detail[12] ); // neutron
	dist->column( "eecl_7",  eecl_detail[13] ); dist->column( "neecl_7",  neecl_detail[13] ); // neutral baryon
	dist->column( "eecl_8",  eecl_detail[14] ); dist->column( "neecl_8",  neecl_detail[14] ); // gamma from pi0 decay from D
	dist->column( "eecl_9",  eecl_detail[15] ); dist->column( "neecl_9",  neecl_detail[15] ); // gamma from pi0 decay from D*
	dist->column( "eecl_10", eecl_detail[16] ); dist->column( "neecl_10", neecl_detail[16] ); // gamma from pi0 decay from D**
	dist->column( "eecl_11", eecl_detail[17] ); dist->column( "neecl_11", neecl_detail[17] ); // gamma from pi0 decay from Ds(*)
	dist->column( "eecl_12", eecl_detail[18] ); dist->column( "neecl_12", neecl_detail[18] ); // gamma from pi0 decay from B
	dist->column( "eecl_13", eecl_detail[19] ); dist->column( "neecl_13", neecl_detail[19] ); // gamma from pi0 decay from others
	dist->column( "eecl_14", eecl_detail[20] ); dist->column( "neecl_14", neecl_detail[20] ); // gamma from eta decay
	dist->column( "eecl_15", eecl_detail[21] ); dist->column( "neecl_15", neecl_detail[21] ); // gamma from D*
	dist->column( "eecl_16", eecl_detail[22] ); dist->column( "neecl_16", neecl_detail[22] ); // brems from semileptonic B decay
	dist->column( "eecl_17", eecl_detail[23] ); dist->column( "neecl_17", neecl_detail[23] ); // gamma from others
	dist->column( "eecl_18", eecl_detail[24] ); dist->column( "neecl_18", neecl_detail[24] ); // K track
	dist->column( "eecl_19", eecl_detail[25] ); dist->column( "neecl_19", neecl_detail[25] ); // K shower
	dist->column( "eecl_20", eecl_detail[26] ); dist->column( "neecl_20", neecl_detail[26] ); // pi track
	dist->column( "eecl_21", eecl_detail[27] ); dist->column( "neecl_21", neecl_detail[27] ); // pi shower
	dist->column( "eecl_22", eecl_detail[28] ); dist->column( "neecl_22", neecl_detail[28] ); // electron track
	dist->column( "eecl_23", eecl_detail[29] ); dist->column( "neecl_23", neecl_detail[29] ); // electron shower
	dist->column( "eecl_24", eecl_detail[30] ); dist->column( "neecl_24", neecl_detail[30] ); // muon track
	dist->column( "eecl_25", eecl_detail[31] ); dist->column( "neecl_25", neecl_detail[31] ); // muon shower
	dist->column( "eecl_26", eecl_detail[32] ); dist->column( "neecl_26", neecl_detail[32] ); // proton track
	dist->column( "eecl_27", eecl_detail[33] ); dist->column( "neecl_27", neecl_detail[33] ); // proton shower
	dist->column( "eecl_28", eecl_detail[34] ); dist->column( "neecl_28", neecl_detail[34] ); // charged baryon track
	dist->column( "eecl_29", eecl_detail[35] ); dist->column( "neecl_29", neecl_detail[35] ); // charged baryon shower
	dist->column( "eecl_30", eecl_detail[36] ); dist->column( "neecl_30", neecl_detail[36] ); // gamma from pi0 decay from D** (forward)
	dist->column( "eecl_31", eecl_detail[37] ); dist->column( "neecl_31", neecl_detail[37] ); // gamma from pi0 decay from D** (barrel)
	dist->column( "eecl_32", eecl_detail[38] ); dist->column( "neecl_32", neecl_detail[38] ); // gamma from pi0 decay from D** (backward)

	cal_Mmiss_Evis( dist, trk_list, gam_list, "mmiss",  "evis"  ); // mmiss  evis
	cal_Mmiss_Evis( dist, Dtau1,    Dtau2,    "mmiss2", "evis2" ); // mmiss2 evis2

	// +++++++++++++++++++++++++++++++ [ remaining particles ] ++++++++++++++++++++++++++++++++++++++++++++++++++++
	int remtrk = cnt_remain_trk  ( Dtau1, Dtau2           );
	int rempi0 = cnt_remain_pi0  ( Dtau1, Dtau2, pi0_list );
	int remks  = cnt_remain_ks   ( Dtau1, Dtau2, ks_list  );
	int remgam = cnt_remain_gamma( Dtau1, Dtau2, gam_list );
	if( remtrk || rempi0 || remks ){
	  if( fl_message ) std::cout << " -> remaining particles exist : "
				     << "Ntrk = " << remtrk << ", "
				     << "Npi0 = " << rempi0 << ", "
				     << "Nks = "  << remks
				     << std::endl;
	  //continue; // removed @20140919
	  //if( remtrk || remks ) continue; // added @20140919
	  if( remtrk + remks > 2 ) continue; // modified @20150130
	}
	  
	dist->column( "remtrk", remtrk );
	dist->column( "rempi0", rempi0 );
	dist->column( "remks",  remks  );
	dist->column( "remgam", remgam );
	// +++++++++++++++++++++++++++++++++++++++++++++
	// pi0 study
	dist->column( "rempi0_0", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test0) );
	dist->column( "rempi0_1", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test1) );
	dist->column( "rempi0_2", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test2) );
	dist->column( "rempi0_3", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test3) );
	dist->column( "rempi0_4", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test4) );
	dist->column( "rempi0_5", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test5) );
	dist->column( "rempi0_6", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test6) );
	dist->column( "rempi0_7", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test7) );
	dist->column( "rempi0_8", cnt_remain_pi0(Dtau1, Dtau2, pi0_list_test8) );

	// +++++++++++++++++++++++++++++++++++++ [ momentum ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	dist->column( "l1pcm",   info_lep1.Vcm().vect().mag() );
	dist->column( "l2pcm",   info_lep2.Vcm().vect().mag() );
	dist->column( "l1p",     lep1.p().vect().mag()        );
	dist->column( "l2p",     lep2.p().vect().mag()        );
	dist->column( "l1pt",    lep1.p().perp()              );
	dist->column( "l2pt",    lep2.p().perp()              );
	dist->column( "l1pc",    lep1.p().cosTheta()          ); // @20141202
	dist->column( "l2pc",    lep2.p().cosTheta()          ); // @20141202

	if( info_Dstr1.dm()==-1 ){ // B -> D  lnu
	  dist->column( "dst1pcm", -1 );
	  dist->column( "dst1p",   -1 );
	  dist->column( "acc1p",   -1 );
	  dist->column( "acc1pc",  -1 ); // @ 20141202
	  dist->column( "acc1m",   -1 ); // @ 20141105
	}else{ // B -> D* lnu
	  dist->column( "dst1pcm", info_Dstr1.Vcm().vect().mag()   );
	  dist->column( "dst1p",   Dstr1.p().vect().mag()          );
	  dist->column( "acc1p",   Dstr1.child(1).p().vect().mag() );
	  dist->column( "acc1pc",  Dstr1.child(1).p().cosTheta()   ); // @20141202
	  UserInfo& info_accpi0 = dynamic_cast<UserInfo&>( Dstr1.child(1).userInfo() );
	  dist->column( "acc1m",   info_accpi0.m_org()             ); // @ 20141105
	}
	
	if( info_Dstr2.dm()==-1 ){ // B -> D  lnu
	  dist->column( "dst2pcm", -1 );
	  dist->column( "dst2p",   -1 );
	  dist->column( "acc2p",   -1 );
	  dist->column( "acc2pc",  -1 ); // @ 20141202
	  dist->column( "acc2m",   -1 ); // @ 20141105
	}else{ // B -> D* lnu
	  dist->column( "dst2pcm", info_Dstr2.Vcm().vect().mag()   );
	  dist->column( "dst2p",   Dstr2.p().vect().mag()          );
	  dist->column( "acc2p",   Dstr2.child(1).p().vect().mag() );
	  dist->column( "acc2pc",  Dstr2.child(1).p().cosTheta()   ); // @20141202
	  UserInfo& info_accpi0 = dynamic_cast<UserInfo&>( Dstr2.child(1).userInfo() );
	  dist->column( "acc2m",   info_accpi0.m_org()             ); // @ 20141105
	}
	
	dist->column( "d1pcm",   info_D1.Vcm().vect().mag() );
	dist->column( "d2pcm",   info_D2.Vcm().vect().mag() );
	dist->column( "d1p",     Dstr1.p().vect().mag()     );
	dist->column( "d2p",     Dstr2.p().vect().mag()     );

	// ++++++++++++++++++++++++++++++++++++++ [ D(*)tau ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	dist->column( "cosdl1", info_Dtau1.cos()      );
	dist->column( "cosdl2", info_Dtau2.cos()      );

	double cosdll = ( info_Dtau1.cos() > info_Dtau2.cos() ? info_Dtau2.cos() : info_Dtau1.cos() ); // low  value
	double cosdlh = ( info_Dtau1.cos() > info_Dtau2.cos() ? info_Dtau1.cos() : info_Dtau2.cos() ); // high value
	dist->column( "cosdll", cosdll );
	dist->column( "cosdlh", cosdlh );

	//dist->column( "kfl1cl",   info_Dtau1.kf_cl()    );
	//dist->column( "kfl1chi2", info_Dtau1.kf_chisq() );
	//dist->column( "kfl1ndf",  info_Dtau1.kf_ndf()   );
	
	//dist->column( "kfl2cl",   info_Dtau2.kf_cl()    );
	//dist->column( "kfl2chi2", info_Dtau2.kf_chisq() );
	//dist->column( "kfl2ndf",  info_Dtau2.kf_ndf()   );
	
	// calculation B directon
	double rm_bdir = calBdirection( Dtau1, Dtau2 );
	dist->column( "rm_bdir",  rm_bdir );

	// ++++++++++++++++++++++++++++ [D* and accomapany particle] ++++++++++++++++++++++++++++++++++++++++++++++++++

	int dst1self = 0;
	int dst2self = 0;
	int acc1self = 0; // accompany particle (e.g. slow pion)
	int acc2self = 0;
	int dst1org  = -1;
	int dst2org  = -1;
	int acc1org  = -1;
	int acc2org  = -1;
	int acc1mo   = 0;
	int acc2mo   = 0;
	int acc1gmo  = 0;
	int acc2gmo  = 0;


	if( info_Dstr1.dm()==-1 ){ // B -> D lnu
	  dist->column( "dm1",       -1 );
	  dist->column( "dst1_morg", -1 );
	  dist->column( "dst1_m",    -1 );
	  dist->column( "kfs1cl",    -1 );
	  dist->column( "kfs1chi2",  -1 );
	  dist->column( "kfs1ndf",   -1 );
	  dst1self = -1;
	  acc1self = -1;
	}else{ // B -> D* lnu
	  dist->column( "dm1",       info_Dstr1.dm()       );
	  dist->column( "dst1_morg", info_Dstr1.m_org()    );
	  dist->column( "dst1_m",    Dstr1.mass()          );
	  dist->column( "kfs1cl",    info_Dstr1.kf_cl()    );
	  dist->column( "kfs1chi2",  info_Dstr1.kf_chisq() );
	  dist->column( "kfs1ndf",   info_Dstr1.kf_ndf()   );
	  dst1self = info_Dstr1.self();
	  dst1org = which_B( Dstr1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  UserInfo& info_acc1 = dynamic_cast<UserInfo&>( Dstr1.child(1).userInfo() );
	  acc1self = info_acc1.self();
	  acc1org = which_B( Dstr1.child(1).genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  if( Dstr1.child(1).genHepevt() && Dstr1.child(1).genHepevt().mother()                                                 ) acc1mo  = Dstr1.child(1).genHepevt().mother().idhep();
	  if( Dstr1.child(1).genHepevt() && Dstr1.child(1).genHepevt().mother() && Dstr1.child(1).genHepevt().mother().mother() ) acc1gmo = Dstr1.child(1).genHepevt().mother().mother().idhep();
	}

	if( info_Dstr2.dm()==-1 ){ // B -> D lnu
	  dist->column( "dm2",       -1 );
	  dist->column( "dst2_morg", -1 );
	  dist->column( "dst2_m",    -1 );
	  dist->column( "kfs2cl",    -1 );
	  dist->column( "kfs2chi2",  -1 );
	  dist->column( "kfs2ndf",   -1 );
	  dst2self = -1;
	  acc2self = -1;
	}else{ // B -> D* lnu
	  dist->column( "dm2",       info_Dstr2.dm()       );
	  dist->column( "dst2_morg", info_Dstr2.m_org()    );
	  dist->column( "dst2_m",    Dstr2.mass()          );
	  dist->column( "kfs2cl",    info_Dstr2.kf_cl()    );
	  dist->column( "kfs2chi2",  info_Dstr2.kf_chisq() );
	  dist->column( "kfs2ndf",   info_Dstr2.kf_ndf()   );
	  dst2self = info_Dstr2.self();
	  dst2org = which_B( Dstr2.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  UserInfo& info_acc2 = dynamic_cast<UserInfo&>( Dstr2.child(1).userInfo() );
	  acc2self = info_acc2.self();
	  acc2org = which_B( Dstr2.child(1).genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	  if( Dstr2.child(1).genHepevt() && Dstr2.child(1).genHepevt().mother()                                                 ) acc2mo  = Dstr2.child(1).genHepevt().mother().idhep();
	  if( Dstr2.child(1).genHepevt() && Dstr2.child(1).genHepevt().mother() && Dstr2.child(1).genHepevt().mother().mother() ) acc2gmo = Dstr2.child(1).genHepevt().mother().mother().idhep();
	}

	dist->column( "dst1self",  dst1self );
	dist->column( "dst1org",   dst1org  );
	dist->column( "acc1self",  acc1self );
	dist->column( "acc1org",   acc1org  );
	dist->column( "acc1mo",    acc1mo   );
	dist->column( "acc1gmo",   acc1gmo  );
	dist->column( "dst2self",  dst2self ); // bug fixed @ 20140306
	dist->column( "dst2org",   dst2org  );
	dist->column( "acc2self",  acc2self );
	dist->column( "acc2org",   acc2org  );
	dist->column( "acc2mo",    acc2mo   );
	dist->column( "acc2gmo",   acc2gmo  );


	// +++++++++++++++++++++++++++++++++++++++++++++++++++

	if( info_Dstr1.rec_mode()==PI0_LUND ){
	  UserInfo& info_pi0 = dynamic_cast<UserInfo&>( Dstr1.child(1).userInfo() );
	  dist->column( "dst1_pi0cos", info_pi0.cos() );
	}else{
	  dist->column( "dst1_pi0cos", 2 );
	}
	if( info_Dstr2.rec_mode()==PI0_LUND ){
	  UserInfo& info_pi0 = dynamic_cast<UserInfo&>( Dstr2.child(1).userInfo() );
	  dist->column( "dst2_pi0cos", info_pi0.cos() );
	}else{
	  dist->column( "dst2_pi0cos", 2 );
	}

	// +++++++++++++++++++++++++++++++++++++++++ [D] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	dist->column( "d1_morg", info_D1.m_org() );
	dist->column( "d2_morg", info_D2.m_org() );
	dist->column( "d1_m",    D1.mass()       );
	dist->column( "d2_m",    D2.mass()       );

	dist->column( "kfd1cl",   info_D1.kf_cl()    );
	dist->column( "kfd1chi2", info_D1.kf_chisq() );
	dist->column( "kfd1ndf",  info_D1.kf_ndf()   );
	dist->column( "kfd2cl",   info_D2.kf_cl()    );
	dist->column( "kfd2chi2", info_D2.kf_chisq() );
	dist->column( "kfd2ndf",  info_D2.kf_ndf()   );

	int d1self = info_D1.self();
	int d2self = info_D2.self();
	dist->column( "d1self", d1self );
	dist->column( "d2self", d2self );

	if( D1.genHepevt() ) dist->column( "d1moid", D1.genHepevt().mother().idhep() );
	else                 dist->column( "d1moid", 0                               );
	if( D2.genHepevt() ) dist->column( "d2moid", D2.genHepevt().mother().idhep() );
	else                 dist->column( "d2moid", 0                               );

	int d1org = which_B( D1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	int d2org = which_B( D2.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );

	dist->column( "d1org", d1org );
	dist->column( "d2org", d2org );

	// ++++++++++++++++++++++++++++++++++++++ [lepton] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int l1org = which_B( lep1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	int l2org = which_B( lep2.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );

	int l1self;

	if( lep1.lund()==info_lep1.id() ){ // true lepton
	  if     ( (gen_b_decay_info[0][1]==1 || gen_b_decay_info[0][1]==2) && lep1.genHepevt().get_ID()==gen_b_decay_info[0][ 2] ) l1self = 2; // tag-side true lepton
	  else if( (gen_b_decay_info[1][1]==1 || gen_b_decay_info[1][1]==2) && lep1.genHepevt().get_ID()==gen_b_decay_info[1][ 2] ) l1self = 2; // tag-side true lepton
	  else if(                               gen_b_decay_info[0][1]==3  && lep1.genHepevt().get_ID()==gen_b_decay_info[0][22] ) l1self = 1; // sig-side true lepton
	  else if(                               gen_b_decay_info[1][1]==3  && lep1.genHepevt().get_ID()==gen_b_decay_info[1][22] ) l1self = 1; // sig-side true lepton
	  else if( info_D1.self()==1 ){ // true D
	    if( d1org==l1org )                                                                                                      l1self = 3; // random   true lepton from same     D side
	    else                                                                                                                    l1self = 5; // random   true lepton from opposite D side
	  }else                                                                                                                     l1self = 7; // random   true lepton (fake D)
	}else{ // fake lepton
	  if( info_D1.self()==1 ){ // true D
	    if( d1org==l1org ) l1self = -3; // same D side
	    else               l1self = -5; // opposite D side
	  }else                l1self = -7; // fake D
	}

	int l2self;

	if( lep2.lund()==info_lep2.id() ){ // true lepton
	  if     ( (gen_b_decay_info[0][1]==1 || gen_b_decay_info[0][1]==2) && lep2.genHepevt().get_ID()==gen_b_decay_info[0][ 2] ) l2self = 2; // tag-side true lepton
	  else if( (gen_b_decay_info[1][1]==1 || gen_b_decay_info[1][1]==2) && lep2.genHepevt().get_ID()==gen_b_decay_info[1][ 2] ) l2self = 2; // tag-side true lepton
	  else if(                               gen_b_decay_info[0][1]==3  && lep2.genHepevt().get_ID()==gen_b_decay_info[0][22] ) l2self = 1; // sig-side true lepton
	  else if(                               gen_b_decay_info[1][1]==3  && lep2.genHepevt().get_ID()==gen_b_decay_info[1][22] ) l2self = 1; // sig-side true lepton
	  else if( info_D2.self()==1 ){ // true D
	    if( d2org==l2org )                                                                                                      l2self = 3; // random   true lepton from same     D side
	    else                                                                                                                    l2self = 5; // random   true lepton from opposite D side
	  }else                                                                                                                     l2self = 7; // random   true lepton (fake D)
	}else{ // fake lepton
	  if( info_D2.self()==1 ){ // true D
	    if( d2org==l2org ) l2self = -3; // same D side
	    else               l2self = -5; // opposite D side
	  }else                l2self = -7; // fake D
	}
	
	dist->column( "l1self",   l1self               );
	dist->column( "l2self",   l2self               );
	dist->column( "l1selfid", info_lep1.id()       );
	dist->column( "l2selfid", info_lep2.id()       );
	dist->column( "l1moid",   info_lep1.motherid() );
	dist->column( "l2moid",   info_lep2.motherid() );
	dist->column( "l1org",    l1org                );
	dist->column( "l2org",    l2org                );



	if     ( info_Dtau1.rec_mode() == 10 ) dist->column( "l1pid",    info_lep1.eidProb()        ); // 10(e)
	else if( info_Dtau1.rec_mode() ==  1 ) dist->column( "l1pid",    info_lep1.muonLikelihood() ); //  1(mu)
	if     ( info_Dtau2.rec_mode() == 10 ) dist->column( "l2pid",    info_lep2.eidProb()        ); // 10(e)
	else if( info_Dtau2.rec_mode() ==  1 ) dist->column( "l2pid",    info_lep2.muonLikelihood() ); //  1(mu)


	// +++++++++++++++++++++++++++++++ [D's daughter particles] +++++++++++++++++++++++++++++++++++++++++++++++++++
	dist->column( "d1nch", D1.nChildren() );
	for( int i=0; i<4; i++ ){
	  int reclund =  -1;
	  int genlund =  -1;
	  int moid    =  -1;
	  int gmoid   =  -1;
	  int org     =  -1;
	  double mom  =  -1; // @20141202
	  double cos  = -10; // @20141202
	  if( i < D1.nChildren() ){
	    Particle& Ddaughter = D1.child(i);
	    reclund = Ddaughter.lund();
	    if( Ddaughter.genHepevt() ){
	      genlund = Ddaughter.genHepevt().idhep();
	      if( Ddaughter.genHepevt().mother() ) moid = Ddaughter.genHepevt().mother().idhep();
	      else                                 moid = 0;
	      if( Ddaughter.genHepevt().mother() && Ddaughter.genHepevt().mother().mother() ) gmoid = Ddaughter.genHepevt().mother().mother().idhep();
	      else                                                                            gmoid = 0;
	      org = which_B( Ddaughter.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	    }
	    mom = Ddaughter.p().vect().mag();
	    cos = Ddaughter.p().cosTheta();
	  }

	  if      ( i==0 ){ dist->column( "d1ch0recid", reclund ); dist->column( "d1ch0selfid", genlund ); dist->column( "d1ch0mo", moid ); dist->column( "d1ch0gmo", gmoid ); dist->column( "d1ch0org", org ); dist->column( "d1ch0p", mom ); dist->column( "d1ch0pc", cos );
	  }else if( i==1 ){ dist->column( "d1ch1recid", reclund ); dist->column( "d1ch1selfid", genlund ); dist->column( "d1ch1mo", moid ); dist->column( "d1ch1gmo", gmoid ); dist->column( "d1ch1org", org ); dist->column( "d1ch1p", mom ); dist->column( "d1ch1pc", cos );
	  }else if( i==2 ){ dist->column( "d1ch2recid", reclund ); dist->column( "d1ch2selfid", genlund ); dist->column( "d1ch2mo", moid ); dist->column( "d1ch2gmo", gmoid ); dist->column( "d1ch2org", org ); dist->column( "d1ch2p", mom ); dist->column( "d1ch2pc", cos );
	  }else if( i==3 ){ dist->column( "d1ch3recid", reclund ); dist->column( "d1ch3selfid", genlund ); dist->column( "d1ch3mo", moid ); dist->column( "d1ch3gmo", gmoid ); dist->column( "d1ch3org", org ); dist->column( "d1ch3p", mom ); dist->column( "d1ch3pc", cos );
	  }
	}

	// ++++++++++++++

	dist->column( "d2nch", D2.nChildren() );
	for( int i=0; i<4; i++ ){
	  int reclund =  -1;
	  int genlund =  -1;
	  int moid    =  -1;
	  int gmoid   =  -1;
	  int org     =  -1;
	  double mom  =  -1; // @20141202
	  double cos  = -10; // @20141202
	  if( i < D2.nChildren() ){
	    Particle& Ddaughter = D2.child(i);
	    reclund = Ddaughter.lund();
	    if( Ddaughter.genHepevt() ){
	      genlund = Ddaughter.genHepevt().idhep();
	      if( Ddaughter.genHepevt().mother() ) moid = Ddaughter.genHepevt().mother().idhep();
	      else                                 moid = 0;
	      if( Ddaughter.genHepevt().mother() && Ddaughter.genHepevt().mother().mother() ) gmoid = Ddaughter.genHepevt().mother().mother().idhep();
	      else                                                                            gmoid = 0;
	      org = which_B( Ddaughter.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	    }
	    mom = Ddaughter.p().vect().mag();
	    cos = Ddaughter.p().cosTheta();
	  }

	  if      ( i==0 ){ dist->column( "d2ch0recid", reclund ); dist->column( "d2ch0selfid", genlund ); dist->column( "d2ch0mo", moid ); dist->column( "d2ch0gmo", gmoid ); dist->column( "d2ch0org", org ); dist->column( "d2ch0p", mom ); dist->column( "d2ch0pc", cos );
	  }else if( i==1 ){ dist->column( "d2ch1recid", reclund ); dist->column( "d2ch1selfid", genlund ); dist->column( "d2ch1mo", moid ); dist->column( "d2ch1gmo", gmoid ); dist->column( "d2ch1org", org ); dist->column( "d2ch1p", mom ); dist->column( "d2ch1pc", cos );
	  }else if( i==2 ){ dist->column( "d2ch2recid", reclund ); dist->column( "d2ch2selfid", genlund ); dist->column( "d2ch2mo", moid ); dist->column( "d2ch2gmo", gmoid ); dist->column( "d2ch2org", org ); dist->column( "d2ch2p", mom ); dist->column( "d2ch2pc", cos );
	  }else if( i==3 ){ dist->column( "d2ch3recid", reclund ); dist->column( "d2ch3selfid", genlund ); dist->column( "d2ch3mo", moid ); dist->column( "d2ch3gmo", gmoid ); dist->column( "d2ch3org", org ); dist->column( "d2ch3p", mom ); dist->column( "d2ch3pc", cos );
	  }
	}

	// +++++++++++++++++++++++++++++++++++++++++++++++++++

	double k1pidmax = -1;
	double k1pidmin = 2;
	for( int i=0; i<D1.nChildren(); i++ ){
	  if( abs(D1.child(i).lund())==Kplus_LUND ){
	    UserInfo& info_k = dynamic_cast<UserInfo&>( D1.child(i).userInfo() );
	    if( k1pidmax < info_k.selKPI() ) k1pidmax = info_k.selKPI();
	    if( k1pidmin > info_k.selKPI() ) k1pidmin = info_k.selKPI();
	  }
	}
	if( k1pidmax==-1 ) k1pidmax = 2;
	dist->column( "k1pidmin", k1pidmin );
	dist->column( "k1pidmax", k1pidmax );

	double k2pidmax = -1;
	double k2pidmin = 2;
	for( int i=0; i<D2.nChildren(); i++ ){
	  if( abs(D2.child(i).lund())==Kplus_LUND ){
	    UserInfo& info_k = dynamic_cast<UserInfo&>( D2.child(i).userInfo() );
	    if( k2pidmax < info_k.selKPI() ) k2pidmax = info_k.selKPI();
	    if( k2pidmin > info_k.selKPI() ) k2pidmin = info_k.selKPI();
	  }
	}
	if( k2pidmax==-1 ) k2pidmax = 2;
	dist->column( "k2pidmin", k2pidmin );
	dist->column( "k2pidmax", k2pidmax );

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
	// additional pi0 information
	if( info_D1.rec_mode() > 999 ){
	  for( int i=0; i<D1.nChildren(); i++ ){
	    if( D1.child(i).lund()==PI0_LUND ){
	      UserInfo& info_pi0 = dynamic_cast<UserInfo&>( D1.child(i).userInfo() );
	      dist->column( "d1_pi0cos", info_pi0.cos()   );
	      dist->column( "pi01m",     info_pi0.m_org() );
	      dist->column( "pi01self",  info_pi0.self()  );
	    }
	  }
	}else{
	  dist->column( "d1_pi0cos",  2 );
	  dist->column( "pi01m",     -1 );
	  dist->column( "pi01self",   1 );
	}

	if( info_D2.rec_mode() > 999 ){
	  for( int i=0; i<D2.nChildren(); i++ ){
	    if( D2.child(i).lund()==PI0_LUND ){
	      UserInfo& info_pi0 = dynamic_cast<UserInfo&>( D2.child(i).userInfo() );
	      dist->column( "d2_pi0cos", info_pi0.cos()   );
	      dist->column( "pi02m",     info_pi0.m_org() );
	      dist->column( "pi02self",  info_pi0.self()  );
	    }
	  }
	}else{
	  dist->column( "d2_pi0cos",  2 );
	  dist->column( "pi02m",     -1 );
	  dist->column( "pi02self",   1 );
	}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
	double d1_pi0_mass = -1;
	for( int i=0; i<D1.nChildren(); i++ ){
	  if( D1.child(i).lund()==PI0_LUND ){
	    UserInfo& info_pi0 = dynamic_cast<UserInfo&>( D1.child(i).userInfo() );
	    d1_pi0_mass = info_pi0.m_org();
	  }
	}
	dist->column( "pi01m", d1_pi0_mass );

	double d2_pi0_mass = -1;
	for( int i=0; i<D2.nChildren(); i++ ){
	  if( D2.child(i).lund()==PI0_LUND ){
	    UserInfo& info_pi0 = dynamic_cast<UserInfo&>( D2.child(i).userInfo() );
	    d2_pi0_mass = info_pi0.m_org();
	  }
	}
	dist->column( "pi02m", d2_pi0_mass );

	// ++++++++++++++++++++++++++++++++ [bremsstrahlung study] ++++++++++++++++++++++++++++++++++++++++++++++++++++
	if( mc==1 ){
	  // <bremsstrahlung from generator inforamation>
	  dist->column( "l1geg", radgam_energy(lep1.genHepevt()) );
	  dist->column( "l2geg", radgam_energy(lep2.genHepevt()) );
	  dist->column( "l1gng", radgam_cnt   (lep1.genHepevt()) );
	  dist->column( "l2gng", radgam_cnt   (lep2.genHepevt()) );
	  // +++++++++++++++++++++++++++++++
	  Hep3Vector lep1_mom_gen( lep1.genHepevt().PX(), lep1.genHepevt().PY(), lep1.genHepevt().PZ() );
	  Hep3Vector lep2_mom_gen( lep2.genHepevt().PX(), lep2.genHepevt().PY(), lep2.genHepevt().PZ() );
	  dist->column( "l1ag", lep1.p3().angle(lep1_mom_gen) ); // deviation angle between rec and gen
	  dist->column( "l2ag", lep2.p3().angle(lep2_mom_gen) );
	  // +++++++++++++++++++++++++++++++
	  // <reconstructed bremsstrahlung gamma>
	  dist->column( "l1gnr", info_lep1.radgam_n() );
	  dist->column( "l2gnr", info_lep2.radgam_n() );
	  // +++++++++++++++++++++++++++++++
	  // <truely reconstructed bremsstrahlung gamma>
	  int    cnt_true_radgam1   = 0;
	  int    cnt_true_radgam2   = 0;
	  double cnt_true_radgam1_e = 0;
	  double cnt_true_radgam2_e = 0;
	  for( int i=0; i<info_lep1.radgam_n(); i++ ){
	    if( Dtau1.child(2+i).genHepevt() && lep1.genHepevt().mother().get_ID() == Dtau1.child(2+i).genHepevt().mother().get_ID() ){
	      cnt_true_radgam1   += 1;
	      cnt_true_radgam1_e += Dtau1.child(2+i).mdstEcl().energy();
	    }
	  }
	  for( int i=0; i<info_lep2.radgam_n(); i++ ){
	    if( Dtau2.child(2+i).genHepevt() && lep2.genHepevt().mother().get_ID() == Dtau2.child(2+i).genHepevt().mother().get_ID() ){
	      cnt_true_radgam2   += 1;
	      cnt_true_radgam2_e += Dtau2.child(2+i).mdstEcl().energy();
	    }
	  }
	  dist->column( "l1gntr", cnt_true_radgam1   );
	  dist->column( "l2gntr", cnt_true_radgam2   );
	  dist->column( "l1getr", cnt_true_radgam1_e );
	  dist->column( "l2getr", cnt_true_radgam2_e );
	  // +++++++++++++++++++++++++++++++
	  //<detected bremsstrahlung gamma>
	  int    cnt_det_radgam1   = 0;
	  int    cnt_det_radgam2   = 0;
	  double cnt_det_radgam1_e = 0;
	  double cnt_det_radgam2_e = 0;
	  
	  dist->column( "rgd1e", -10 );
	  dist->column( "rgd1a", -10 );
	  dist->column( "rgd1r",   0 );
	  dist->column( "rgd2e", -10 );
	  dist->column( "rgd2a", -10 );
	  dist->column( "rgd2r",   0 );
	  dist->column( "rgd3e", -10 );
	  dist->column( "rgd3a", -10 );
	  dist->column( "rgd3r",   0 );
	  dist->column( "rgd4e", -10 );
	  dist->column( "rgd4a", -10 );
	  dist->column( "rgd4r",   0 );
	  
	  
	  for( std::vector<Particle>::iterator gam = gam_list.begin(); gam != gam_list.end(); gam++ ){
	    if( !gam->genHepevt() ) continue;
	    // ++++++++++++++++++++++++++++++++++++
	    double tmp_angle = -10;
	    if( gam->genHepevt().mother().get_ID() == lep1.genHepevt().mother().get_ID() ){
	      cnt_det_radgam1   += 1;
	      cnt_det_radgam1_e += gam->mdstEcl().energy();
	      tmp_angle = lep1.p3().angle( gam->p3() );
	    }else if( gam->genHepevt().mother().get_ID() == lep2.genHepevt().mother().get_ID() ){
	      cnt_det_radgam2   += 1;
	      cnt_det_radgam2_e += gam->mdstEcl().energy();
	      tmp_angle = lep2.p3().angle( gam->p3() );
	    }
	    // ++++++++++++++++++++++++++++++++++++
	    int tmp_rec = 0;
	    for( int i=0; i<info_lep1.radgam_n(); i++ ){
	      if( Dtau1.child(2+i).genHepevt() && Dtau1.child(2+i).genHepevt().get_ID() == gam->genHepevt().get_ID() ){
		tmp_rec = 1;
	      }
	    }
	    for( int i=0; i<info_lep2.radgam_n(); i++ ){
	      if( Dtau2.child(2+i).genHepevt() && Dtau2.child(2+i).genHepevt().get_ID() == gam->genHepevt().get_ID() ){
		tmp_rec = 1;
	      }
	    }
	    // ++++++++++++++++++++++++++++++++++++
	    if( gam->genHepevt().mother().get_ID() == lep1.genHepevt().mother().get_ID() ||
		gam->genHepevt().mother().get_ID() == lep2.genHepevt().mother().get_ID() ){
	      int cnt = cnt_det_radgam1 + cnt_det_radgam2;
	      
	      if( cnt==1 ){
		dist->column( "rgd1e", gam->mdstEcl().energy() );
		dist->column( "rgd1a", tmp_angle               );
		dist->column( "rgd1r", tmp_rec                 );
	      }else if( cnt==2 ){
		dist->column( "rgd2e", gam->mdstEcl().energy() );
		dist->column( "rgd2a", tmp_angle               );
		dist->column( "rgd2r", tmp_rec                 );
	      }else if( cnt==3 ){
		dist->column( "rgd3e", gam->mdstEcl().energy() );
		dist->column( "rgd3a", tmp_angle               );
		dist->column( "rgd3r", tmp_rec                 );
	      }else{
		dist->column( "rgd4e", gam->mdstEcl().energy() );
		dist->column( "rgd4a", tmp_angle               );
		dist->column( "rgd4r", tmp_rec                 );
	      }
	    }
	  }
	  
	  dist->column( "l1ged", cnt_det_radgam1_e );
	  dist->column( "l2ged", cnt_det_radgam2_e );
	  dist->column( "l1gnd", cnt_det_radgam1   );
	  dist->column( "l2gnd", cnt_det_radgam2   );
	}
	// +++++++++++++++++++++++++++++++
	// make more bremsstrahlung recovery
	dist->column( "cosdla1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.05, MUminus_LUND) );
	dist->column( "cosdla2",  cosBDl_with_more_bremsstrahlung       (Dtau2,        gam_list, 0.05, MUminus_LUND) );
	dist->column( "rm_bdira", calBdirection_with_more_bremsstrahlung(Dtau1, Dtau2, gam_list, 0.05, MUminus_LUND) );
	dist->column( "cosdlb1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.10, MUminus_LUND) );
	dist->column( "cosdlb2",  cosBDl_with_more_bremsstrahlung       (Dtau2,        gam_list, 0.10, MUminus_LUND) );
	dist->column( "rm_bdirb", calBdirection_with_more_bremsstrahlung(Dtau1, Dtau2, gam_list, 0.10, MUminus_LUND) );
	dist->column( "cosdlc1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.10              ) );
	dist->column( "cosdlc2",  cosBDl_with_more_bremsstrahlung       (Dtau2,        gam_list, 0.10              ) );
	dist->column( "rm_bdirc", calBdirection_with_more_bremsstrahlung(Dtau1, Dtau2, gam_list, 0.10              ) );
	dist->column( "cosdld1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.20, MUminus_LUND) );
	dist->column( "cosdld2",  cosBDl_with_more_bremsstrahlung       (Dtau2,        gam_list, 0.20, MUminus_LUND) );
	dist->column( "rm_bdird", calBdirection_with_more_bremsstrahlung(Dtau1, Dtau2, gam_list, 0.20, MUminus_LUND) );
	dist->column( "cosdle1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.20              ) );
	dist->column( "cosdle2",  cosBDl_with_more_bremsstrahlung       (Dtau2,        gam_list, 0.20              ) );
	dist->column( "rm_bdire", calBdirection_with_more_bremsstrahlung(Dtau1, Dtau2, gam_list, 0.20              ) );

	// ++++++++++++++++++++++++++++++++++++++ [D** study] +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ @ 20140703
	dist->column( "rm_ddst_angle", cal_nearestgammadirection( Dstr1, Dstr2, Dtau1, Dtau2, gam_list ) ); // not useful
	Particle* veto_pi0 = NULL;
	double cosdl1_ctrl = info_Dtau1.cos();
	double cosdl2_ctrl = info_Dtau2.cos();
	if( info_Dstr1.dm()>0 && info_Dstr1.rec_mode()==PI0_LUND ){
	  veto_pi0    = &Dstr1.child(1); // if B1 side is D*l and slow pion is neutral, ....
	  cosdl1_ctrl = info_Dtau1.cos_ctrl();
	}else if( info_Dstr2.dm()>0 && info_Dstr2.rec_mode()==PI0_LUND ){
	  veto_pi0 = &Dstr2.child(1); // if B2 side is D*l and slow pion is neutral, ....
	  cosdl2_ctrl = info_Dtau2.cos_ctrl();
	}
	
	cal_Mmiss_Evis( dist, trk_list, gam_list, "mmiss_ctrl",  "evis_ctrl",  veto_pi0 );
	cal_Mmiss_Evis( dist, Dtau1,    Dtau2,    "mmiss2_ctrl", "evis2_ctrl", veto_pi0 ); // (strictly speaking) slow pi should be estimated from D* - D due to kinematic fitting.
	
	dist->column( "eecl_ctrl", calEclEnergyWithMatch2GoodGamma(Dtau1, Dtau2, !true, NULL, NULL, veto_pi0, gen_b_decay_info[0][0], gen_b_decay_info[1][0] ) );
	
	dist->column( "cosdl1_ctrl", cosdl1_ctrl );
	dist->column( "cosdl2_ctrl", cosdl2_ctrl );

	dist->column( "cosdll_ctrl", ( cosdl1_ctrl > cosdl2_ctrl ? cosdl2_ctrl : cosdl1_ctrl ) ); // low  value
	dist->column( "cosdlh_ctrl", ( cosdl1_ctrl > cosdl2_ctrl ? cosdl1_ctrl : cosdl2_ctrl ) ); // high value

	// +++++++++++++++++++++++++++++++++ [ gen-rec matching ] +++++++++++++++++++++++++++++++++++++++++++++++++++++
	int recfB_id = 0;
	int recsB_id = 0;
	int recgen   = 0;
	if( info_D1.self()==1 ){ // true D1
	  if( d1org==1 ){
	    recfB_id = 0;
	    recsB_id = 1;
	    recgen   = 1;
	  }else if( d1org==2 ){
	    recfB_id =  1;
	    recsB_id =  0;
	    recgen   = -1;
	  }
	}else if( info_D2.self()==1 ){ // true D2
	  if( d2org==1 ){
	    recfB_id =  1;
	    recsB_id =  0;
	    recgen   = -1;
	  }else if( d2org==2 ){
	    recfB_id = 0;
	    recsB_id = 1;
	    recgen   = 1;
	  }
	}else{ // fake D1 and D2
	  recfB_id = 0;
	  recsB_id = 1;
	  recgen   = 0;
	}

	dist->column( "recgen",  recgen ); // 1(rec1=gen1,rec2=gen2), -1(rec1=gen2,rec2=gen1), 0(rec1=gen1,rec2=gen2)

	// +++++++++++++++++++++++++++++++++++ [ gen-mode(BB) ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int gm_dd   = make_flag_dd  ( gen_b_decay_info[0][13], gen_b_decay_info[1][13] );
	int gm_semi = make_flag_semi( gen_b_decay_info[0][ 1], gen_b_decay_info[1][ 1] );
	dist->column( "gm_dd",   gm_dd   );
	dist->column( "gm_semi", gm_semi );

	// +++++++++++++++++++++++++++++++++++ [ gen-mode(B) ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int nb1  = gen_b_decay_info[0][3];
	int nb2  = gen_b_decay_info[1][3];
	int nb1g = gen_b_decay_info[0][4];
	int nb2g = gen_b_decay_info[1][4];
	dist->column( "nb1",       nb1                     ); dist->column( "nb2",       nb2                     );
    	dist->column( "nb1g",      nb1g                    ); dist->column( "nb2g",      nb2g                    );
	dist->column( "nb1pip",    gen_b_decay_info[0][29] ); dist->column( "nb2pip",    gen_b_decay_info[1][29] );
	dist->column( "nb1pi0",    gen_b_decay_info[0][30] ); dist->column( "nb2pi0",    gen_b_decay_info[1][30] );
	dist->column( "semi1",     gen_b_decay_info[0][ 1] ); dist->column( "semi2",     gen_b_decay_info[1][ 1] );
	dist->column( "nrootd1",   gen_b_decay_info[0][ 5] ); dist->column( "nrootd2",   gen_b_decay_info[1][ 5] );
	dist->column( "rootdf1",   gen_b_decay_info[0][ 6] ); dist->column( "rootdf2",   gen_b_decay_info[1][ 6] ); 
	dist->column( "rootds1",   gen_b_decay_info[0][ 8] ); dist->column( "rootds2",   gen_b_decay_info[1][ 8] ); // fixed @ 20150123
	dist->column( "nd1",       gen_b_decay_info[0][10] ); dist->column( "nd2",       gen_b_decay_info[1][10] );
	int gm_ddst1 = gen_b_decay_info[0][13];
	int gm_ddst2 = gen_b_decay_info[1][13];
	dist->column( "gm_ddst1",  gm_ddst1                ); dist->column( "gm_ddst2",  gm_ddst2                );
	dist->column( "dst1_acc",  gen_b_decay_info[0][14] ); dist->column( "dst2_acc",  gen_b_decay_info[1][14] );
	dist->column( "dst1_d",    gen_b_decay_info[0][27] ); dist->column( "dst2_d",    gen_b_decay_info[1][27] ); // 20141002
	dist->column( "cc1",       gen_b_decay_info[0][16] ); dist->column( "cc2",       gen_b_decay_info[1][16] );

	dist->column( "rootd1nc",  gen_b_decay_info[0][31] ); dist->column( "rootd2nc",  gen_b_decay_info[1][31] ); // 20141028
	dist->column( "fldstst1",  gen_b_decay_info[0][32] ); dist->column( "fldstst2",  gen_b_decay_info[1][32] ); // 20141028
	dist->column( "incl1gmc",  gen_b_decay_info[0][33] ); dist->column( "incl2gmc",  gen_b_decay_info[1][33] ); // 20141028
	dist->column( "nddst1pp",  gen_b_decay_info[0][34] ); dist->column( "nddst2pp",  gen_b_decay_info[1][34] ); // 20141028
	dist->column( "nddst1p0",  gen_b_decay_info[0][35] ); dist->column( "nddst2p0",  gen_b_decay_info[1][35] ); // 20141028
	dist->column( "flmcddst",  flag_DststMC );

	dist->column( "taulep1",   gen_b_decay_info[0][21] ); dist->column( "taulep2",   gen_b_decay_info[1][21] ); // 20140925
	dist->column( "ddst1_d",   gen_b_decay_info[0][23] ); dist->column( "ddst2_d",   gen_b_decay_info[1][23] ); // 20140703
	dist->column( "ddst1_acc", gen_b_decay_info[0][25] ); dist->column( "ddst2_acc", gen_b_decay_info[1][25] ); // 20140703

	// +++++++++++++++++++++++++++++++++++ [ gen-mode(D) ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	dist->column( "gm_df1hr",  gen_d_mode_info [0][0][0]  ); dist->column( "gm_df2hr",  gen_d_mode_info [1][0][0] );
	dist->column( "gm_ds1hr",  gen_d_mode_info [0][1][0]  ); dist->column( "gm_ds2hr",  gen_d_mode_info [1][1][0] );
	dist->column( "gm_df1hnr", gen_d_mode_info [0][0][1]  ); dist->column( "gm_df2hnr", gen_d_mode_info [1][0][1] );
	dist->column( "gm_ds1hnr", gen_d_mode_info [0][1][1]  ); dist->column( "gm_ds2hnr", gen_d_mode_info [1][1][1] );
	dist->column( "gm_df1lep", gen_d_mode_info [0][0][2]  ); dist->column( "gm_df2lep", gen_d_mode_info [1][0][2] );
	dist->column( "gm_ds1lep", gen_d_mode_info [0][1][2]  ); dist->column( "gm_ds2lep", gen_d_mode_info [1][1][2] );
	dist->column( "gm_df1nu",  gen_d_mode_info [0][0][3]  ); dist->column( "gm_df2nu",  gen_d_mode_info [1][0][3] );
	dist->column( "gm_ds1nu",  gen_d_mode_info [0][1][3]  ); dist->column( "gm_ds2nu",  gen_d_mode_info [1][1][3] );
	dist->column( "gm_df1gam", gen_d_mode_info [0][0][4]  ); dist->column( "gm_df2gam", gen_d_mode_info [1][0][4] );
	dist->column( "gm_ds1gam", gen_d_mode_info [0][1][4]  ); dist->column( "gm_ds2gam", gen_d_mode_info [1][1][4] );
	
	dist->column( "gm_t1hr",  gen_tau_mode_info [0][0]  ); dist->column( "gm_t2hr",  gen_tau_mode_info [1][0] );
	dist->column( "gm_t1hnr", gen_tau_mode_info [0][1]  ); dist->column( "gm_t2hnr", gen_tau_mode_info [1][1] );
	dist->column( "gm_t1lep", gen_tau_mode_info [0][2]  ); dist->column( "gm_t2lep", gen_tau_mode_info [1][2] );
	dist->column( "gm_t1nu",  gen_tau_mode_info [0][3]  ); dist->column( "gm_t2nu",  gen_tau_mode_info [1][3] );
	dist->column( "gm_t1gam", gen_tau_mode_info [0][4]  ); dist->column( "gm_t2gam", gen_tau_mode_info [1][4] );
	dist->column( "gm_t1pro", gen_tau_mode_info [0][5]  ); dist->column( "gm_t2pro", gen_tau_mode_info [1][5] );

	// ++++++++++++++++++++++++++++++++++++++++ [    D** mass(gen)     ] ++++++++++++++++++++++++++++++++++++++++++
	Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
	if( gen_b_decay_info[0][13]==3 ) dist->column( "ddst1m", genMgr( Panther_ID(gen_b_decay_info[0][7]) ).M() );
	else                             dist->column( "ddst1m", -10                                              );
	if( gen_b_decay_info[1][13]==3 ) dist->column( "ddst2m", genMgr( Panther_ID(gen_b_decay_info[1][7]) ).M() );
	else                             dist->column( "ddst2m", -10                                              );

	// ++++++++++++++++++++++++++++++++++++++++ [kinematics for weight ] ++++++++++++++++++++++++++++++++++++++++++
	// for MC model correction, calculate q^2, p_l*, w, cos(theta_l) in B->D(**) l nu.
	for( int i=0; i<2; i++ ){
	  double q2     = -10;
	  double plep   = -10;
	  double w      = -10;
	  double coslep = -10;
	  if( (gen_b_decay_info[i][1]==1 || gen_b_decay_info[i][1]==2) &&
	      (
	       ( gen_b_decay_info[i][13]==1 ) || // B -> D  l nu
	       ( gen_b_decay_info[i][13]==2 ) || // B -> D* l nu
	       ( gen_b_decay_info[i][13]==3 && gen_b_decay_info[i][32]==1 ) )
	      ){
	    calKinematics( gen_b_decay_info[i][0], gen_b_decay_info[i][7], gen_b_decay_info[i][2], q2, plep, w, coslep );
	  }
	  if( fl_message ) std::cout << "Kinematics Calculation " << i << " : "
				     << "q2 = "     << q2   << ", "
				     << "plep = "   << plep << ", "
				     << "w = "      << w    << ", "
				     << "coslep = " << coslep
				     << std::endl;
	  if( i==0 ){
	    dist->column( "k1_q2",  q2     );
	    dist->column( "k1_pl",  plep   );
	    dist->column( "k1_w",   w      );
	    dist->column( "k1_cos", coslep );
	  }else{
	    dist->column( "k2_q2",  q2     );
	    dist->column( "k2_pl",  plep   );
	    dist->column( "k2_w",   w      );
	    dist->column( "k2_cos", coslep );
	  }
	}
	
	// ++++++++++++++++++++++++++++++++++++++++ [ self ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int self = -10;
	/*
	if     ( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) && (nb1-nb1g==3&&nb2-nb2g==3) && ( dst1self!=0&&dst2self!=0 && d1self==1&&d2self==1 ) && ( l1self>0&&l2self>0 ) && ( l1self*l2self==1 || l1self*l2self==2 ) ) self =  1; // sig
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) && (nb1-nb1g==3&&nb2-nb2g==3) && ( dst1self!=0&&dst2self!=0 && d1self==1&&d2self==1 ) && ( l1self>0&&l2self>0 ) && ( l1self*l2self==4                     ) ) self =  2; // tag
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) && (nb1-nb1g==3&&nb2-nb2g==3)                                                                                                                               ) self =  0; // fake
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) &&!(nb1-nb1g==3&&nb2-nb2g==3)                                                                                                                               ) self = -1; // non-resonant semileptonic
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi==0 )                                                                                                                                                             ) self = -2; // hadronic decay
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi==0 )                                                                                                                                                             ) self = -3; // other DD(hadronic decay)
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi!=0 && gm_dd!=0   )                                                                                                                                               ) self = -4; // other DD(D(*)D(*))
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi!=0 && gm_dd==0   ) && ( gm_ddst1==4||gm_ddst2==4 )                                                                                                               ) self = -5; // other DD(double D)
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi!=0 && gm_dd==0   ) && ( gm_ddst1==3||gm_ddst2==3 )                                                                                                               ) self = -6; // other DD(D**)
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi!=0 && gm_dd==0   )                                                                                                                                               ) self = -7; // other DD(other)
	else if( ( abs(event_type)!=rm_bb&&abs(event_type)>2 )                                                                                                                                                                      ) self = -8; // other BB(charged)
	else if( ( abs(event_type)!=rm_bb&&abs(event_type)<3 )                                                                                                                                                                      ) self = -9; // continuum
	else  std::cout << "[self check] (exp,run,event) = "
			<< expNo << ", "
			<< runNo << ", "
			<< evtNo << std::endl;
	*/
	if     ( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) && (nb1-nb1g==3&&nb2-nb2g==3) && ( dst1self!=0&&dst2self!=0 && d1self==1&&d2self==1 ) && ( l1self>0&&l2self>0 ) && ( l1self*l2self==1 || l1self*l2self==2 ) ) self =  1; // sig
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) && (nb1-nb1g==3&&nb2-nb2g==3) && ( dst1self!=0&&dst2self!=0 && d1self==1&&d2self==1 ) && ( l1self>0&&l2self>0 ) && ( l1self*l2self==4                     ) ) self =  2; // tag
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) && (nb1-nb1g==3&&nb2-nb2g==3)                                                                                                                               ) self =  0; // fake
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi!=0 ) &&!(nb1-nb1g==3&&nb2-nb2g==3)                                                                                                                               ) self = -1; // non-resonant semileptonic
	else if( ( abs(event_type)==rm_bb&&gm_dd==rm_dd && gm_semi==0 )                                                                                                                                                             ) self = -2; // hadronic decay
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi==0 )                                                                                                                                                             ) self = -3; // other DD(hadronic decay) -> [double DD]
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi!=0 && gm_dd!=0   )                                                                                                                                               ) self = -4; // other DD(D(*)D(*))
	else if( ( abs(event_type)==rm_bb&&gm_dd!=rm_dd && gm_semi!=0 && gm_dd==0   )                                                                                                                                               ) self = -5; // other DD(other) -> [D**]
	else if( ( abs(event_type)!=rm_bb&&abs(event_type)>2 )                                                                                                                                                                      ) self = -6; // other BB(charged)
	else if( ( abs(event_type)!=rm_bb&&abs(event_type)<3 )                                                                                                                                                                      ) self = -7; // continuum
	else  std::cout << "[self check] (exp,run,event) = "
			<< expNo << ", "
			<< runNo << ", "
			<< evtNo << std::endl;
	dist->column( "self", self );
	/*
	if( info_D1.rec_mode()!=1201 && info_D1.rec_mode()!=1102 && // veto dirty modes
	    info_D2.rec_mode()!=1201 && info_D2.rec_mode()!=1102 &&
	    !(info_Dstr1.rec_mode()==111 && info_Dstr2.rec_mode()==111) && // veto double slow pi0
	    rempi0==0 && // no remainint pi0
	    Dtau1.charge()+Dtau2.charge() == 0 && // charge check
	    !(info_Dtau1.flavor() && info_Dtau2.flavor() && info_Dtau1.flavor()*info_Dtau2.flavor()==1) && // flavor check
	    (
	     (info_Dtau1.cos() > info_Dtau2.cos() && info_Dstr1.Vcm().vect().mag() < 2.5 && info_Dstr2.Vcm().vect().mag() < 2.0) ||
	     (info_Dtau1.cos() < info_Dtau2.cos() && info_Dstr1.Vcm().vect().mag() < 2.0 && info_Dstr2.Vcm().vect().mag() < 2.5)
	     ) && // D* momentum cut
	    ( -2 < cosdlh && cosdlh < 1.5 ) // cos(theta)B-Dl cut
	    // slightly loose mass cut
	    ) ckEclEnergyWithMatch2GoodGamma( Dtau1, Dtau2, self, neecl_detail[0]+neecl_detail[1]+neecl_detail[2] );
	*/
	dist->dumpData();
	cnt_dump++;
	if( fl_message ) std::cout << " -> dump !" << std::endl;
      }
    }
    if( fl_message ) std::cout << "AnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAnaAna" << std::endl;
    return cnt_dump;
  }

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
