#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Rec_D_2body( std::vector<Particle>& d_list,
			       std::vector<Particle>& k_list,
			       std::vector<Particle>& p_list,
			       const bool fl_message
			       )
  {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( !k_list.size() ) return;
    if( !p_list.size() ) return;
    int rec_mode_d = 0;
    int cnt        = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if     ( abs(k_list.begin()->lund())==Kplus_LUND  && abs(p_list.begin()->lund())==PIplus_LUND ) rec_mode_d =  101; // [charged K  + charged pi ]
    else if( abs(k_list.begin()->lund())==Kplus_LUND  && abs(p_list.begin()->lund())==   PI0_LUND ) rec_mode_d = 1001; // [charged K  +         pi0] -> it will be rejected by flavor-charged correlation
    else if( abs(k_list.begin()->lund())==   Ks_LUND  && abs(p_list.begin()->lund())==PIplus_LUND ) rec_mode_d =  110; // [        Ks + charged pi ]
    else if( abs(k_list.begin()->lund())==   Ks_LUND  && abs(p_list.begin()->lund())==   PI0_LUND ) rec_mode_d = 1010; // [        Ks +         pi0]
    else     std::cerr << "[ABORT] Wrong recD-mode(2body)" << std::endl, abort();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message ) std::cout << "D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)" << std::endl;
    for( std::vector<Particle>::iterator k = k_list.begin(); k != k_list.end(); k++ ){
      for( std::vector<Particle>::iterator p = p_list.begin(); p != p_list.end(); p++ ){
	cnt++;
	UserInfo& info_k = dynamic_cast<UserInfo&>( k->userInfo() );
	UserInfo& info_p = dynamic_cast<UserInfo&>( p->userInfo() );
	if( fl_message ){
	  std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt << " : ";
	  rec_message( k ); std::cout << ", ";
	  rec_message( p );
	}
	int chg  = int( k->charge() + p->charge() ); // total charge
	int kchg = int( k->charge()               ); // kaon  charge
	
	// charge check
	if( abs(chg) > 1 ){
	  if( fl_message ) std::cout << " -> not allowed charge : " << chg << std::endl;
	  continue;
	}
	
	// Check duplication
	if( !check_dupli_daughter(k,p) ){
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
	  if( kchg < 0 ){ // D0 -> K- pi+
	    d_lund   = D0_LUND;
	    d_flavor = 1;
	  }else if( kchg > 0 ){ // anti-D0 -> K+ pi-
	    d_lund = antiD0_LUND;
	    d_flavor = -1;
	  }else if( kchg == 0 ){ // (anti-)D0 -> Ks pi0 // unknown flavor
	    d_lund   = D0_LUND; // assign D0_LUND in the case of unknown flavor
	    d_flavor = 0;       // unknown flavor
	  }
	}else if( chg==1 ){ // D+
	  d_lund   = Dplus_LUND;
	  d_flavor = 1;
	}else if( chg==-1){ // D-
	  d_lund   = Dminus_LUND;
	  d_flavor = -1;
	}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	
	// reconstruction
	Particle D( k->p() + p->p(), Ptype(d_lund) );
	
	D.relation().append( *k );
	D.relation().append( *p );
	setUserInfo( D );
	UserInfo& info = dynamic_cast<UserInfo&>( D.userInfo() );
	info.rec_mode( rec_mode_d );
	info.cntid   ( cnt        );
	info.flavor  ( d_flavor   );
	//info.self    ( check_selfR(D)      );
	info.self    ( check_selfR2(D)      );
	info.m_org( D.mass() );

	kinematic_fit( D, 1 ); // vertex-constrined fit

	if( !masscut( D, masscut_D_L, masscut_D_H, chg ? PDG_DplusMass : PDG_D0Mass ) ){
	  if( fl_message ) std::cout << " -> mass cut rejection : " << D.momentum().mass() << std::endl;
	  continue;
	}
	HepLorentzVector d_4Vcm = D.p();
	d_4Vcm.boost( cmboost );
	info.Vcm( d_4Vcm );

	if( fl_message ) std::cout << " -> push_back : "
				   << "LUND="     << D.lund()        << ", "
				   << "rec_mode=" << info.rec_mode() << ", "
				   << "d_flavor=" << info.flavor()   << ", "
				   << "self="     << info.self()
				   << ")"         << std::endl;
	d_list.push_back( D );

      }
    }
    if( fl_message ) std::cout << "D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)" << std::endl;
  }



  void DSTRTAUNU::Rec_D_2body( std::vector<Particle>& d_list,
			       std::vector<Particle>& p_list, // K+/pi+
			       const bool fl_message
			       )
  {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( !p_list.size() ) return;
    int rec_mode_d = 0;
    int cnt        = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if     ( abs(p_list.begin()->lund())== Kplus_LUND ) rec_mode_d =    2; // [charged K  + charged K  ]
    else if( abs(p_list.begin()->lund())==PIplus_LUND ) rec_mode_d =  200; // [charged pi + charged pi ]
    else std::cerr << "[ABORT] Wrong recD-mode(2body)" << std::endl, abort();
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message ) std::cout << "D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)" << std::endl;
    for( std::vector<Particle>::iterator p1 = p_list.begin(); p1 != p_list.end(); p1++ ){
      for( std::vector<Particle>::iterator p2 = p1+1; p2 != p_list.end(); p2++ ){
	cnt++;
	UserInfo& info_p1 = dynamic_cast<UserInfo&>( p1->userInfo() );
	UserInfo& info_p2 = dynamic_cast<UserInfo&>( p2->userInfo() );
	if( fl_message ){
	  std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt << " : ";
	  rec_message( p1 ); std::cout << ", ";
	  rec_message( p2 );
	}
	int chg  = int( p1->charge() + p2->charge() ); // total charge
	int kchg = 0;                                  // kaon  charge
	
	// charge check
	if( abs(chg) > 1 ){
	  if( fl_message ) std::cout << " -> not allowed charge : " << chg << std::endl;
	  continue;
	}
	
	// Check duplication
	if( !check_dupli_daughter(p1,p2) ){
	  if( fl_message ) std::cout << " -> duplication " << std::endl;
	  continue;
	}

	// flavor-charge check
	if( kchg && chg*kchg==1 ){
	  if( fl_message ) std::cout << " -> not allowed flavor-charge : " << chg << ", " << kchg << std::endl;
	  continue;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	int d_lund   = D0_LUND;
	int d_flavor = 0; // unknown flavor
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	
	// reconstruction
	Particle D( p1->p() + p2->p(), Ptype(d_lund) );
	
	D.relation().append( *p1 );
	D.relation().append( *p2 );
	setUserInfo( D );
	UserInfo& info = dynamic_cast<UserInfo&>( D.userInfo() );
	info.rec_mode( rec_mode_d );
	info.cntid   ( cnt        );
	info.flavor  ( d_flavor   );
	//info.self    ( check_selfR(D)      );
	info.self    ( check_selfR2(D)      );
	info.m_org( D.mass() );

	kinematic_fit( D, 1 ); // vertex-constrined fit

	if( !masscut( D, masscut_D_L, masscut_D_H, chg ? PDG_DplusMass : PDG_D0Mass ) ){
	  if( fl_message ) std::cout << " -> mass cut rejection : " << D.momentum().mass() << std::endl;
	  continue;
	}
	HepLorentzVector d_4Vcm = D.p();
	d_4Vcm.boost( cmboost );
	info.Vcm( d_4Vcm );
	
	if( fl_message ) std::cout << " -> push_back : "
				   << "LUND="     << D.lund()        << ", "
				   << "rec_mode=" << info.rec_mode() << ", "
				   << "d_flavor=" << info.flavor()   << ", "
				   << "self="     << info.self()
				   << ")"         << std::endl;
	d_list.push_back( D );

      }
    }
    if( fl_message ) std::cout << "D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)D(2body)" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
