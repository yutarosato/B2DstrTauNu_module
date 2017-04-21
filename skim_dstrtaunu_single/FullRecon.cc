#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif


  int DSTRTAUNU::FullReconN( const int lund, const double nb_th, const bool fl_nb_cont )
  {
    Ekpfullrecon_Manager&  frecMgr = Ekpfullrecon_Manager::get_manager();
    Brecon_header_Manager& bhMgr   = Brecon_header_Manager::get_manager();
    AnaBrecon              brecon;
    if( bhMgr.count()!=1 ) return 0; // 1 if btag candidate exists.
    int cnt_cand = 0;
    for( Ekpfullrecon_Manager::iterator btag_it = frecMgr.begin(); btag_it != frecMgr.end(); btag_it++ ){
      Ekpfullrecon &btag = *btag_it;
      double nb = (fl_nb_cont ? btag.cont_NBout() : btag.NBout() );
      if( abs(btag.tag_id())==B0_LUND && nb>nb_th ) cnt_cand++;
    }
    return cnt_cand++;
  }

  void DSTRTAUNU::FullRecon( bool fl_message, const bool fl_dump )
  {
    
    Ekpfullrecon_Manager&  frecMgr = Ekpfullrecon_Manager::get_manager();
    Brecon_header_Manager& bhMgr   = Brecon_header_Manager::get_manager();
    AnaBrecon              brecon;
    bool init_fl_message = fl_message;
    if( fl_message ) std::cout << "FULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECO" << std::endl;
    if( bhMgr.count()!=1 ) return; // 1 if btag candidate exists.
    int nbtag = 0;
    for( Ekpfullrecon_Manager::iterator btag_it = frecMgr.begin(); btag_it != frecMgr.end(); btag_it++ ) nbtag++;
    if( fl_message ) std::cout << nbtag << " candidates are found." << std::endl;
    int cntnbtag = 0;
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
      fl_message = (nboutDef>0.005 ? 1&&init_fl_message : 0);

      Particle& bcand = const_cast<Particle&>( brecon.getParticle((int)btag.get_ID()) );
      setUserInfo( bcand );
      UserInfo& info = dynamic_cast<UserInfo&>( bcand.userInfo() );
      int self = check_selfR2(bcand); // check if candidate is true or false : true(1), false(0)
      info.self( self );
      
      int nroot   [5] = {0};
      int rootlund[4] = {0};
      int digit_D [4] = {0};
      int digit_B     = 0;
      FullReconMode( decayMode, DDecayMode1, DDecayMode2, DDecayMode3, DDecayMode4,
		     digit_B, digit_D, rootlund, nroot );
      
      if( fl_message ){
	std::cout << "++++++++++++++++++++++++++++++++++++++++++"
		  << "[ Btag : " << ++cntnbtag << " ]"
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
	std::cout << "selfR2 = " << self << std::endl;
	display_rec_particle( bcand, fl_message );
	std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      }
      // fill ntuple
      ///*
      hadtag_dist->column( "exp",        expNo       );
      hadtag_dist->column( "run",        runNo       );
      hadtag_dist->column( "evt",        evtNo       );
      hadtag_dist->column( "nbtag",      nbtag       );
      hadtag_dist->column( "tmbc",       Mbc         );
      hadtag_dist->column( "tde",        deltaE      );
      hadtag_dist->column( "tmcinfo",    MCinfo      );
      hadtag_dist->column( "tself",      self        );
      hadtag_dist->column( "tblund",     tagID       );
      hadtag_dist->column( "tbdecay",    decayMode   );
      hadtag_dist->column( "tdecay1",    DDecayMode1 );
      hadtag_dist->column( "tdecay2",    DDecayMode2 );
      hadtag_dist->column( "tdecay3",    DDecayMode3 );
      hadtag_dist->column( "tdecay4",    DDecayMode4 );
      hadtag_dist->column( "tnfs",       nFS         );
      hadtag_dist->column( "tnboutdef",  nboutDef    );
      hadtag_dist->column( "tbestdef",   bestDef     );
      hadtag_dist->column( "tnboutcont", nboutCont   );
      hadtag_dist->column( "tbestcont",  bestCont    );
      
      hadtag_dist->column( "tdgb",       digit_B     );
      hadtag_dist->column( "tdgdst1",    digit_D[0]  );
      hadtag_dist->column( "tdgd1",      digit_D[1]  );
      hadtag_dist->column( "tdgdst2",    digit_D[2]  );
      hadtag_dist->column( "tdgd2",      digit_D[3]  );
      hadtag_dist->column( "tdg",        digit_B + digit_D[0] + digit_D[1] + digit_D[2] + digit_D[3] );
      hadtag_dist->column( "tndst",      nroot[0]    );
      hadtag_dist->column( "tnd",        nroot[1]    );
      hadtag_dist->column( "tndsst",     nroot[2]    );
      hadtag_dist->column( "tnds",       nroot[3]    );
      hadtag_dist->column( "tnjpsi",     nroot[4]    );
      hadtag_dist->column( "tdst1lund",  rootlund[0] );
      hadtag_dist->column( "td1lund",    rootlund[1] );
      hadtag_dist->column( "tdst2lund",  rootlund[2] );
      hadtag_dist->column( "td2lund",    rootlund[3] );
      
      if( fl_dump ) hadtag_dist->dumpData();
      //*/
    }
    
    if( fl_message ) std::cout << "FULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECOFULLRECO"   << std::endl;
    return;
  }
  
  int DSTRTAUNU::FullReconMode( const int decay, const int decay1, const int decay2, const int decay3, const int decay4,
				int& digit_B, int digit_D[4], int rootlund[4], int nroot[5]
				){
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // nroot   [5] : [D*, D, Ds*, Ds, Jpsi]
    // rootlund[4] : [D*, D, Ds*, Ds      ]
    // digit_D [4] : [D*, D, Ds*, Ds      ], definition of digit : (gamma,pi0,pi+,Ks,K+)
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // initialize
    for( int i=0; i<5; i++ ) nroot   [i] = 0;
    for( int i=0; i<4; i++ ) rootlund[i] = 0;
    for( int i=0; i<4; i++ ) digit_D [i] = 0;
    digit_B;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
    int input[4] = {0}; // [1st-D(s)*, 1st-D(s), 2nd-D(s)*, 2nd-D(s)]
    if  ( abs(decay1)/1000==413 || abs(decay1)/1000==423 || abs(decay1)/1000==433 ){ // 1st D(s)*
      input[0] = decay1;
      input[1] = decay2;
      if( abs(decay3)/1000==413 || abs(decay3)/1000==423 || abs(decay3)/1000==433 ){ // 2nd D(s)*
	input[2] = decay3;
	input[3] = decay4;
      }else if( abs(decay3)/1000==421 || abs(decay3)/1000==411 || abs(decay3)/1000==431 || abs(decay3)/1000==432 ){ // 2nd D(s)
	input[3] = decay3;
      }      
    }else if( abs(decay1)/1000==421 || abs(decay1)/1000==411 || abs(decay1)/1000==431 || abs(decay1)/1000==432 ){ // 1st D(s)
      input[1] = decay1;
      if( abs(decay2)/1000==413 || abs(decay2)/1000==423 || abs(decay2)/1000==433 ){ // 2nd D(s)*
	input[2] = decay2;
	input[3] = decay3;
      }else if( abs(decay2)/1000==421 || abs(decay2)/1000==411 || abs(decay2)/1000==431 || abs(decay2)/1000==432 ){ // 2nd D(s)
	input[3] = decay2;
      }      
    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
    digit_B = FullReconMode_B( decay, nroot );
    
    for( int i=0; i<4; i++ ){
      if     ( i%2==0 ) digit_D[i] = FullReconMode_Dstr( input[i], rootlund[i] );
      else if( i%2==1 ) digit_D[i] = FullReconMode_D   ( input[i], rootlund[i] );
    }

    return 0;
  }

  int DSTRTAUNU::FullReconMode_B( const int decay, int nroot[5] ){
    int& nDstr  = nroot[0];
    int& nD     = nroot[1];
    int& nDsstr = nroot[2];
    int& nDs    = nroot[3];
    int& nJpsi  = nroot[4];

    if     ( abs(decay)==511624 ){ nDstr =1;           return  100; } //  1 anti-B0 -> D*+ pi-
    else if( abs(decay)==511735 ){ nDstr =1;           return 1100; } //  2 anti-B0 -> D*+ pi- pi0
    else if( abs(decay)==512046 ){ nDstr =1;           return  300; } //  3 anti-B0 -> D*+ pi- pi+ pi-
    else if( abs(decay)==511622 ){ nD    =1;           return  100; } //  4 anti-B0 -> D+  pi-
    else if( abs(decay)==511733 ){ nD    =1;           return 1100; } //  5 anti-B0 -> D+  pi- pi0 (with D*+ veto)
    else if( abs(decay)==512044 ){ nD    =1;           return  300; } //  6 anti-B0 -> D+  pi- pi+ pi-
    else if( abs(decay)==511846 ){ nDsstr=1; nDstr =1; return    0; } //  7 anti-B0 -> D*+ Ds*-
    else if( abs(decay)==511845 ){ nDstr =1; nDs   =1; return    0; } //  8 anti-B0 -> D*+ Ds-
    else if( abs(decay)==511844 ){ nD    =1; nDsstr=1; return    0; } //  9 anti-B0 -> D+  Ds*-
    else if( abs(decay)==511842 ){ nD    =1; nDs   =1; return    0; } // 10 anti-B0 -> D+  Ds-
    else if( abs(decay)==511753 ){ nJpsi =1;           return   10; } // 11 anti-B0 -> J/psi Ks
    else if( abs(decay)==511975 ){ nJpsi =1;           return  101; } // 12 anti-B0 -> J/psi K- pi+
    else if( abs(decay)==512175 ){ nJpsi =1;           return  210; } // 13 anti-B0 -> J/psi Ks pi+ pi-
    else if( abs(decay)==511532 ){ nD    =1;           return 1000; } // 14 anti-B0 -> D0 pi0
    else if( abs(decay)==512157 ){ nDstr =1;           return 1300; } // 15 anti-B0 -> D*+ pi- pi- pi+ pi0
    else if( abs(decay)==521634 ){ nDstr =1;           return  100; } //  1 B+ -> D*0 pi-
    else if( abs(decay)==521745 ){ nDstr =1;           return 1100; } //  2 B+ -> D*0 pi- pi0
    else if( abs(decay)==522056 ){ nDstr =1;           return  300; } //  3 B+ -> D*0 pi- pi- pi+
    else if( abs(decay)==521632 ){ nD    =1;           return  100; } //  4 B+ -> D0  pi-
    else if( abs(decay)==521743 ){ nD    =1;           return 1100; } //  5 B+ -> D0 pi- pi0 (with D*0 veto)
    else if( abs(decay)==522054 ){ nD    =1;           return  300; } //  6 B+ -> D0 pi- pi- pi+
    else if( abs(decay)==521856 ){ nDstr =1; nDsstr=1; return    0; } //  7 B+ -> D*0 Ds*-
    else if( abs(decay)==521855 ){ nDstr =1; nDs   =1; return    0; } //  8 B+ -> D*0 Ds-
    else if( abs(decay)==521854 ){ nD    =1; nDsstr=1; return    0; } //  9 B+ -> D0  Ds*-
    else if( abs(decay)==521852 ){ nD    =1; nDs   =1; return    0; } // 10 B+ -> D0  Ds-
    else if( abs(decay)==521764 ){ nJpsi =1;           return    1; } // 11 B+ -> J/psi K-
    else if( abs(decay)==522186 ){ nJpsi =1;           return  201; } // 12 B+ -> J/psi K- pi+ pi-
    else if( abs(decay)==521742 ){ nD    =1;           return    1; } // 13 B+ -> D0 K-
    else if( abs(decay)==521833 ){ nD    =1;           return  200; } // 14 B+ -> D+ pi- pi-
    else if( abs(decay)==522167 ){ nDstr =1;           return 1300; } // 15 B+ -> D*0 pi- pi- pi+ pi0
    else if( abs(decay)==521875 ){ nJpsi =1;           return 1001; } // 16 B+ -> J/psi K- pi0
    else if( abs(decay)==521964 ){ nJpsi =1;           return  110; } // 17 B+ -> J/psi Ks pi-
    else std::cout << "[WARNING] Wrong B decay mode : " << decay << std::endl;
  }


  int DSTRTAUNU::FullReconMode_Dstr( const int decay, int& lund ){
    lund = decay/1000;
    if     ( abs(decay)==     0 ) return    0;
    else if( abs(decay)==413632 ) return   100; // D*+  -> D0 pi+
    else if( abs(decay)==413522 ) return  1000; // D*+  -> D+ pi0
    else if( abs(decay)==423532 ) return  1000; // D*0  -> D0 pi0
    else if( abs(decay)==423443 ) return 10000; // D*0  -> D0 gamma
    else if( abs(decay)==433453 ) return 10000; // Ds*+ -> Ds gamma
    else std::cout << "[WARNING] Wrong D(s)* decay mode : " << decay << std::endl;
  }

  int DSTRTAUNU::FullReconMode_D( const int decay, int& lund ){
    lund = decay/1000;
    if     ( lund== 432 ) lund =  431;
    else if( lund==-432 ) lund = -431;
    
    if     ( abs(decay)==     0 ) return    0;
    else if( abs(decay)==421532 ) return  101; // D0 -> K- pi+
    else if( abs(decay)==421643 ) return 1101; // D0 -> K- pi+ pi0
    else if( abs(decay)==421954 ) return  301; // D0 -> K- pi+ pi+ pi-
    else if( abs(decay)==421421 ) return 1010; // D0 -> Ks pi0
    else if( abs(decay)==421732 ) return  210; // D0 -> Ks pi+ pi-
    else if( abs(decay)==421952 ) return   12; // D0 -> Ks K+  K-
    else if( abs(decay)==421642 ) return    2; // D0 -> K+ K-
    else if( abs(decay)==421422 ) return  200; // D0 -> pi+ pi-
    else if( abs(decay)==421533 ) return 1200; // D0 -> pi+ pi- pi0
    else if( abs(decay)==421843 ) return 1210; // D0 -> Ks pi+ pi- pi0
    // D+
    else if( abs(decay)==411743 ) return  201; // D+ -> K- pi+ pi+
    else if( abs(decay)==411854 ) return 1201; // D+ -> K- pi+ pi+ pi0
    else if( abs(decay)==411521 ) return  110; // D+ -> Ks pi+
    else if( abs(decay)==411632 ) return 1110; // D+ -> Ks pi+ pi0
    else if( abs(decay)==411943 ) return  310; // D+ -> Ks pi+ pi+ pi-
    else if( abs(decay)==411853 ) return  102; // D+ -> K+ K- pi+
    else if( abs(decay)==411964 ) return 1102; // D+ -> K+ K- pi+ pi0
    // Ds
    else if( abs(decay)==431853 ) return  102; // Ds -> K+ K- pi+
    else if( abs(decay)==431631 ) return   11; // Ds -> Ks K+
    else if( abs(decay)==431633 ) return  300; // Ds -> pi+ pi+ pi-
    else if( abs(decay)==431964 ) return 1102; // Ds -> K+ K- pi+ pi0
    else if( abs(decay)==432054 ) return  211; // Ds -> Ks K- pi+ pi+
    else if( abs(decay)==432053 ) return  211; // Ds -> Ks K+ pi+ pi-
    else if( abs(decay)==431743 ) return  201; // Ds -> K+ pi+ pi-
    else if( abs(decay)==432275 ) return  302; // Ds -> K+ K- pi+ pi+ pi-
    else std::cout << "[WARNING] Wrong D(s) decay mpode : " << decay << std::endl;
  }

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
