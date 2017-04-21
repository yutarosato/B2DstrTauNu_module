#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Rec_D_3body( std::vector<Particle>& d_list,
			       std::vector<Particle>& kaon_list, // Ks/K+
			       std::vector<Particle>& pion_list, // pi+
			       std::vector<Particle>& pi0_list,  // pi0
			       const bool fl_message
			       )
  {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( !kaon_list.size() ) return;
    if( !pion_list.size() ) return;
    if(  !pi0_list.size() ) return;
    int rec_mode_d = 0;
    int cnt        = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if     ( abs(kaon_list.begin()->lund())==Kplus_LUND && abs(pion_list.begin()->lund())==PIplus_LUND && abs(pi0_list.begin()->lund())==PI0_LUND ) rec_mode_d = 1101; // [charged K  + charged pi + pi0 ]
    else if( abs(kaon_list.begin()->lund())==Ks_LUND    && abs(pion_list.begin()->lund())==PIplus_LUND && abs(pi0_list.begin()->lund())==PI0_LUND ) rec_mode_d = 1110; // [        Ks + charged pi + pi0]
    else std::cerr << "[ABORT] Wrong recD-mode(3body)" << std::endl, abort();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message ) std::cout << "D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)" << std::endl;
    for( std::vector<Particle>::iterator kaon = kaon_list.begin(); kaon != kaon_list.end(); kaon++ ){
      for( std::vector<Particle>::iterator pion = pion_list.begin(); pion != pion_list.end(); pion++ ){
	for( std::vector<Particle>::iterator pi0 = pi0_list.begin(); pi0 != pi0_list.end(); pi0++ ){
	  cnt++;
	  UserInfo& info_kaon = dynamic_cast<UserInfo&>( kaon->userInfo() );
	  UserInfo& info_pion = dynamic_cast<UserInfo&>( pion->userInfo() );
	  UserInfo& info_pi0  = dynamic_cast<UserInfo&>( pi0 ->userInfo() );
	  if( fl_message ){
	    std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt << " : ";
	    rec_message( kaon ); std::cout << ", ";
	    rec_message( pion ); std::cout << ", ";
	    rec_message( pi0 );
	  }
	  
	  int chg  = int( kaon->charge() + pion->charge() + pi0->charge() ); // total charge
	  int kchg = int( kaon->charge()                                  ); // kaon charge
	  
	  // charge check
	  if( abs(chg) > 1 ){
	    if( fl_message ) std::cout << " -> not allowed charge : " << chg << std::endl;
	    continue;
	  }
	  
	  // Check duplication
	  if( !check_dupli_daughter(kaon,pion) ){
	    if( fl_message ) std::cout << " -> duplication " << std::endl;
	    continue;
	  }

	  // flavor-charge check
	  if( kchg && chg*kchg==1 ){
	    if( fl_message ) std::cout << " -> not allowed flavor-charge : " << chg << ", " << kchg << std::endl;
	    continue;
	  }
	  
	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	  int d_lund   = 0;
	  int d_flavor = 0;
	  if( chg==0 ){ // neutral D
	    if( kchg < 0 ){ // D0 -> K- pi+ pi0
	      d_lund   = D0_LUND;
	      d_flavor = 1;
	    }else if( kchg > 0 ){ // anti-D0 -> K+ pi- pi0
	      d_lund = antiD0_LUND;
	      d_flavor = -1;
	    }else{ // none
	      std::cerr << "[WARNING] not supported rec-D mode : " << rec_mode_d << std::endl;
	    }
	  }else if( chg==1 ){ // D+ -> Ks pi+ pi0
	    d_lund   = Dplus_LUND;
	    d_flavor = 1;
	  }else if( chg==-1){ // D- -> Ks pi- pi0
	    d_lund   = Dminus_LUND;
	    d_flavor = -1;
	  }
	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	  
	  // reconstruction
	  Particle D( kaon->p() + pion->p() + pi0->p(), Ptype(d_lund) );
	  
	  D.relation().append( *kaon );
	  D.relation().append( *pion );
	  D.relation().append( *pi0  );
	  setUserInfo( D );
	  UserInfo& info = dynamic_cast<UserInfo&>( D.userInfo() );
	  info.rec_mode( rec_mode_d );
	  info.cntid   ( cnt        );
	  info.flavor  ( d_flavor   );
	  info.m_org   ( D.momentum().mass() );
	  //info.self    ( check_selfR(D)      );
	  info.self    ( check_selfR2(D)      );
	  if( !masscut( D, masscut_D_L, masscut_D_H, chg ? PDG_DplusMass : PDG_D0Mass ) ){
	    if( fl_message ) std::cout << " -> mass cut rejection : " << D.momentum().mass() << std::endl;
	    continue;
	  }
	  
	HepLorentzVector d_4Vcm = D.p();
	d_4Vcm.boost( cmboost );
	info.Vcm( d_4Vcm );

	  //kinematic_fit( D, 4 ); // mass-constrined fit
	  if( fl_message ) std::cout << " -> push_back : "
				     << "LUND="     << D.lund()        << ", "
				     << "rec_mode=" << info.rec_mode() << ", "
				     << "d_flavor=" << info.flavor()   << ", "
				     << "self="     << info.self()
				     << ")"         << std::endl;
	  d_list.push_back( D );
	}
      }
    }
    if( fl_message ) std::cout << "D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)" << std::endl;
  }


  void DSTRTAUNU::Rec_D_3body( std::vector<Particle>& d_list,
			       std::vector<Particle>& k_list, // K+/Ks/pi+/pi0
			       std::vector<Particle>& p_list, // k+/pi+
			       const bool fl_message
			       )
  {

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( !k_list.size() ) return;
    if( !p_list.size() ) return;
    int rec_mode_d = 0;
    int cnt        = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if     ( abs(k_list.begin()->lund())==    Ks_LUND && abs(p_list.begin()->lund())==PIplus_LUND ) rec_mode_d =  210; // [         Ks + charged pi + charged pi ]
    else if( abs(k_list.begin()->lund())==    Ks_LUND && abs(p_list.begin()->lund())== Kplus_LUND ) rec_mode_d =   12; // [         Ks + charged K  + charged K  ]
    else if( abs(k_list.begin()->lund())== Kplus_LUND && abs(p_list.begin()->lund())==PIplus_LUND ) rec_mode_d =  201; // [ charged K  + charged pi + charged pi ]
    else if( abs(k_list.begin()->lund())==PIplus_LUND && abs(p_list.begin()->lund())== Kplus_LUND ) rec_mode_d =  102; // [ charged pi + charged K  + charged K  ]
    else if( abs(k_list.begin()->lund())==   PI0_LUND && abs(p_list.begin()->lund())==PIplus_LUND ) rec_mode_d = 1200; // [ neutral pi + charged pi + charged pi ]
    else std::cerr << "[ABORT] Wrong recD-mode(3body)" << std::endl, abort();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message ) std::cout << "D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)" << std::endl;
    for( std::vector<Particle>::iterator k = k_list.begin(); k != k_list.end(); k++ ){
      for( std::vector<Particle>::iterator p1 = p_list.begin(); p1 != p_list.end(); p1++ ){
	for( std::vector<Particle>::iterator p2 = p1+1; p2 != p_list.end(); p2++ ){
	  cnt++;
	  UserInfo& info_k  = dynamic_cast<UserInfo&>(  k->userInfo() );
	  UserInfo& info_p1 = dynamic_cast<UserInfo&>( p1->userInfo() );
	  UserInfo& info_p2 = dynamic_cast<UserInfo&>( p2->userInfo() );
	  if( fl_message ){
	    std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt << " : ";
	    rec_message( k  ); std::cout << ", ";
	    rec_message( p1 ); std::cout << ", ";
	    rec_message( p2 );
	  }
	  
	  int chg  = int( k->charge() + p1->charge() + p2->charge() );      // total charge
	  int kchg = ( abs(k->lund())==Kplus_LUND ? (int)k->charge() : 0 ); // kaon  charge
	  
	  // charge check
	  if( abs(chg) > 1 ){
	    if( fl_message ) std::cout << " -> not allowed charge : " << chg << std::endl;
	    continue;
	  }
	  
	  // Check duplication
	  if( !check_dupli_daughter(k,p1) || !check_dupli_daughter(k,p2) ){
	    if( fl_message ) std::cout << " -> duplication " << std::endl;
	    continue;
	  }

	  // flavor-charge check
	  if( kchg && chg*kchg==1 ){
	    if( fl_message ) std::cout << " -> not allowed flavor-charge : " << chg << ", " << kchg << std::endl;
	    continue;
	  }
	  
	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	  int d_lund   = 0;
	  int d_flavor = 0;
	  if( chg==0 ){ // neutral D
	    if( k->charge() ){ // none
	      std::cerr << "[WARNING] not supported rec-D mode : " << rec_mode_d << std::endl;
	    }else{
	      d_lund   = D0_LUND; 
	      d_flavor = 0; // unknown flavor
	    }
	  }else if( chg==1 ){ // D+ -> K+ K- pi+ or D+ -> pi+ pi- K+
	    d_lund   = Dplus_LUND;
	    d_flavor = 1;
	  }else if( chg==-1){ // D- -> K+ K- pi- or D- -> pi+ pi- K-
	    d_lund   = Dminus_LUND;
	    d_flavor = -1;
	  }
	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	  
	  // reconstruction
	  Particle D( k->p() + p1->p() + p2->p(), Ptype(d_lund) );
	  
	  D.relation().append( *k  );
	  D.relation().append( *p1 );
	  D.relation().append( *p2 );
	  setUserInfo( D );
	  UserInfo& info = dynamic_cast<UserInfo&>( D.userInfo() );
	  info.rec_mode( rec_mode_d );
	  info.cntid   ( cnt        );
	  info.flavor  ( d_flavor   );
	  info.m_org   ( D.momentum().mass() );
	  //info.self    ( check_selfR(D)      );
	  info.self    ( check_selfR2(D)      );
	  if( !masscut( D, masscut_D_L, masscut_D_H, chg ? PDG_DplusMass : PDG_D0Mass ) ){
	    if( fl_message ) std::cout << " -> mass cut rejection : " << D.momentum().mass() << std::endl;
	    continue;
	  }
	  
	  HepLorentzVector d_4Vcm = D.p();
	  d_4Vcm.boost( cmboost );
	  info.Vcm( d_4Vcm );
	  
	  //kinematic_fit( D, 4 ); // mass-constrined fit
	  if( fl_message ) std::cout << " -> push_back : "
				     << "LUND="     << D.lund()        << ", "
				     << "rec_mode=" << info.rec_mode() << ", "
				     << "d_flavor=" << info.flavor()   << ", "
				     << "self="     << info.self()
				     << ")"         << std::endl;
	  d_list.push_back( D );
	}
      }
    }
    if( fl_message ) std::cout << "D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)D(3body)" << std::endl;

  }

  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
