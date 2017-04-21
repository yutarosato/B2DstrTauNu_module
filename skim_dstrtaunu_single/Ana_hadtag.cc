#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  int DSTRTAUNU::Ana_hadtag( std::vector<Particle>& dtau_list,
			     int gen_b_decay_info[][40], 
			     int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
			     std::vector<Particle>& trk_list, std::vector<Particle>& gam_list, std::vector<Particle>& lep_list,
			     std::vector<Particle>& pi0_list, std::vector<Particle>& ks_list, 
			     std::vector<Particle>& pi0_list_test0, std::vector<Particle>& pi0_list_test1, std::vector<Particle>& pi0_list_test2, std::vector<Particle>& pi0_list_test3,
			     std::vector<Particle>& pi0_list_test4, std::vector<Particle>& pi0_list_test5, std::vector<Particle>& pi0_list_test6, std::vector<Particle>& pi0_list_test7, std::vector<Particle>& pi0_list_test8,
			     BelleTuple* dist, const bool fl_message )
  {
    int cnt_dump = 0;
    if( fl_message ) std::cout << "ANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANA" << std::endl;
    Ekpfullrecon_Manager&  frecMgr = Ekpfullrecon_Manager::get_manager();
    Brecon_header_Manager& bhMgr   = Brecon_header_Manager::get_manager();
    AnaBrecon              brecon;
    if( bhMgr.count()!=1 ) return 0; // 1 if btag candidate exists.
    for( Ekpfullrecon_Manager::iterator btag_it = frecMgr.begin(); btag_it != frecMgr.end(); btag_it++ ){
      Ekpfullrecon &btag = *btag_it;
      int    tagID       = btag.tag_id(); // LUND code of B
      int    decayMode   = btag.decay();  // B decay mode ( e.g. D0 pi- pi- pi+  : 522054 )
      int    DDecayMode1 = btag.Ddec(0);  // D decay mode
      int    DDecayMode2 = btag.Ddec(1);
      int    DDecayMode3 = btag.Ddec(2);
      int    DDecayMode4 = btag.Ddec(3);
      int    MCinfo      = btag.MCinfo(); // for NeuroBayes trainings. (not for end-user)
      double deltaE      = btag.DeltaE();
      double Mbc         = btag.Mbc();
      double nboutDef    = btag.NBout();
      int    bestDef     = btag.NBRank();
      double nboutCont   = btag.cont_NBout();
      int    bestCont    = btag.cont_NBRank();
      int    nFS         = btag.nFS();
      if( !(bestDef==1 && abs(tagID)==B0_LUND && nboutDef>0.005) ) continue; // only select best B0 candidate.
      Particle& bcand = const_cast<Particle&>( brecon.getParticle((int)btag.get_ID()) );
      setUserInfo( bcand );
      UserInfo& info = dynamic_cast<UserInfo&>( bcand.userInfo() );
      info.self( check_selfR2(bcand) ); // check if candidate is true or false : true(1), false(0)
      HepLorentzVector b_4Vcm = bcand.p();
      b_4Vcm.boost( cmboost );
      info.Vcm( b_4Vcm );
      
      int nroot   [5] = {0};
      int rootlund[4] = {0};
      int digit_D [4] = {0};
      int digit_B     = 0;
      FullReconMode( decayMode, DDecayMode1, DDecayMode2, DDecayMode3, DDecayMode4,
		     digit_B, digit_D, rootlund, nroot );
      if( fl_message ){
	std::cout << "++++++++++++++++++++++++++++++++++++++++++"
		  << "[ Btag]"
		  << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "B_LUND(" << tagID << "), "
		  << "BdecayMode(" << std::setw(7) << std::right << decayMode   << "), "
		  << "DDecayMode("
		  << std::setw(8) << std::right << DDecayMode1 << ", "
		  << std::setw(8) << std::right << DDecayMode2 << ", "
		  << std::setw(8) << std::right << DDecayMode3 << ", "
		  << std::setw(8) << std::right << DDecayMode4 << "), "
		  << "MCinfo("    << std::setw(2) << std::right << MCinfo << "), "
		  << "(dE,Mbc)=("
		  << std::setw(8) << std::right << deltaE << ", "
		  << std::setw(8) << std::right << Mbc    << "), "
		  << "(NBDef,nbNBCont)=("
		  << std::setw(10) << std::right << nboutDef  << ", "
		  << std::setw(10) << std::right << nboutCont << "), "
		  << "(RankDef,RankCont)=(" 
		  << std::setw(2) << std::right << bestDef  << ", "
		  << std::setw(2) << std::right << bestCont << "), "
		  << "nFS(" << std::setw(2) << std::right << nFS << ")"
		  << std::endl;
	std::cout << "selfR2 = " << info.self() << std::endl;
	display_rec_particle( bcand, fl_message );
	std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      }      
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for( std::vector<Particle>::iterator dtau = dtau_list.begin(); dtau != dtau_list.end(); dtau++ ){
	Particle& Dtau      = *dtau;
	UserInfo& info_Dtau = dynamic_cast<UserInfo&>( Dtau.userInfo() );

	if( fl_message ) std::cout << "Dtau_CNTID = " << std::setw(3) << std::right << info_Dtau.cntid();
	// charge check
	if( Dtau.charge()+bcand.charge() != 0 ){
	  if( fl_message ) std::cout << " -> not allowed charge : " << Dtau.charge() << " + " << bcand.charge() << std::endl;
	  continue;
	}
	
	// Check duplication
	if( !check_dupli_daughter( Dtau, bcand) ){
	  if( fl_message ) std::cout << " -> duplication " << std::endl;
	  continue;
	}
	
	// flavor check
	//if( info_Dtau1.flavor() && info_Dtau2.flavor() && info_Dtau1.flavor()*info_Dtau2.flavor()==1 ){
	//if( fl_message ) std::cout << " -> not allowed flavor : " << info_Dtau1.flavor() << ", " << info_Dtau2.flavor() << std::endl;
	//continue; // removed @20140919
	//}
	
	Ana_single_sub( dtau, gen_b_decay_info, gen_d_mode_info, gen_tau_mode_info, gam_list, lep_list, dist, false, fl_message );
	dist->column( "tmbc",       Mbc         );
	dist->column( "tde",        deltaE      );
	dist->column( "tmcinfo",    MCinfo      );
	dist->column( "tself",      info.self() );
	dist->column( "tblund",     tagID       );
	dist->column( "tbdecay",    decayMode   );
	dist->column( "tdecay1",    DDecayMode1 );
	dist->column( "tdecay2",    DDecayMode2 );
	dist->column( "tdecay3",    DDecayMode3 );
	dist->column( "tdecay4",    DDecayMode4 );
	dist->column( "tnfs",       nFS         );
	dist->column( "tnboutdef",  nboutDef    );
	dist->column( "tbestdef",   bestDef     );
	dist->column( "tnboutcont", nboutCont   );
	dist->column( "tbestcont",  bestCont    );

	dist->column( "tdgb",       digit_B     );
	dist->column( "tdgdst1",    digit_D[0]  );
	dist->column( "tdgd1",      digit_D[1]  );
	dist->column( "tdgdst2",    digit_D[2]  );
	dist->column( "tdgd2",      digit_D[3]  );
	dist->column( "tdg",        digit_B + digit_D[0] + digit_D[1] + digit_D[2] + digit_D[3] );
	dist->column( "tndst",      nroot[0]    );
	dist->column( "tnd",        nroot[1]    );
	dist->column( "tndsst",     nroot[2]    );
	dist->column( "tnds",       nroot[3]    );
	dist->column( "tnjpsi",     nroot[4]    );
	dist->column( "tdst1lund",  rootlund[0] );
	dist->column( "td1lund",    rootlund[1] );
	dist->column( "tdst2lund",  rootlund[2] );
	dist->column( "td2lund",    rootlund[3] );
	// ++++++++++++++++++++++++++++++++++ [ event variable ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double eecl_detail [39]={0};
	int    neecl_detail[39]={0};
	dist->column( "eecl",    calEclEnergyWithMatch2GoodGamma(Dtau, bcand, !true, eecl_detail, neecl_detail, NULL, gen_b_decay_info[0][0], gen_b_decay_info[1][0] ) );
	dist->column( "neecl",   neecl_detail[0]+neecl_detail[1]+neecl_detail[2] );
	cal_Mmiss_Evis( dist, trk_list, gam_list, "mmiss",  "evis"  ); // mmiss  evis
	cal_Mmiss_Evis( dist, Dtau,     bcand,    "mmiss2", "evis2" ); // mmiss2 evis2

	// +++++++++++++++++++++++++++++++ [ remaining particles ] ++++++++++++++++++++++++++++++++++++++++++++++++++++
	int remtrk = cnt_remain_trk  ( Dtau, bcand           );
	int rempi0 = cnt_remain_pi0  ( Dtau, bcand, pi0_list );
	int remks  = cnt_remain_ks   ( Dtau, bcand, ks_list  );
	int remgam = cnt_remain_gamma( Dtau, bcand, gam_list );
	if( remtrk || rempi0 || remks ){
	  if( fl_message ) std::cout << " -> remaining particles exist : "
				     << "Ntrk = " << remtrk << ", "
				     << "Npi0 = " << rempi0 << ", "
				     << "Nks = "  << remks
				     << std::endl;
	  if( remtrk || remks ) continue;
	}
	  
	dist->column( "remtrk", remtrk );
	dist->column( "rempi0", rempi0 );
	dist->column( "remks",  remks  );
	dist->column( "remgam", remgam );
	// +++++++++++++++++++++++++++++++++++++++++++++
	// pi0 study
	dist->column( "rempi0_0", cnt_remain_pi0(Dtau, bcand, pi0_list_test0) );
	dist->column( "rempi0_1", cnt_remain_pi0(Dtau, bcand, pi0_list_test1) );
	dist->column( "rempi0_2", cnt_remain_pi0(Dtau, bcand, pi0_list_test2) );
	dist->column( "rempi0_3", cnt_remain_pi0(Dtau, bcand, pi0_list_test3) );
	dist->column( "rempi0_4", cnt_remain_pi0(Dtau, bcand, pi0_list_test4) );
	dist->column( "rempi0_5", cnt_remain_pi0(Dtau, bcand, pi0_list_test5) );
	dist->column( "rempi0_6", cnt_remain_pi0(Dtau, bcand, pi0_list_test6) );
	dist->column( "rempi0_7", cnt_remain_pi0(Dtau, bcand, pi0_list_test7) );
	dist->column( "rempi0_8", cnt_remain_pi0(Dtau, bcand, pi0_list_test8) );


	dist->dumpData();
	cnt_dump++;
	if( fl_message ) std::cout << " -> dump !" << std::endl;
      }
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
    if( fl_message ) std::cout << "ANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANAANA" << std::endl;
    return cnt_dump;
  }

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
