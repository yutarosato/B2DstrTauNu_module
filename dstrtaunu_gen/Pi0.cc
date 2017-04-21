#include "belle.h"
#include "DSTRTAUNU.h"

#include <vector>
#include <iostream>
#include <stdlib.h>

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Pi0_cand( std::vector<Particle>& pi0_list, const bool fl_message )
  {
    int cnt = 0;
    if( fl_message ) std::cout << "PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0" << std::endl;
    Mdst_pi0_Manager& pi0Mgr = Mdst_pi0_Manager::get_manager();
    for( Mdst_pi0_Manager::iterator p = pi0Mgr.begin(); p != pi0Mgr.end(); p++ ){
      cnt++;
      // make particle
      Particle pi0( (*p), true );
      double morg = (pi0.child(0).p() + pi0.child(1).p()).m();
      if( !(morg - PDG_PI0Mass < masscut_pi0_H && morg - PDG_PI0Mass >= masscut_pi0_L) ) continue;
      setUserInfo( pi0 );
      UserInfo& info = dynamic_cast<UserInfo&>(pi0.userInfo());
      info.cntid( cnt );
      
      
      Particle& gam0 = pi0.child(0);
      Particle& gam1 = pi0.child(1);
      
      if( gam0.e() < gam_E_cut ) continue;
      if( gam1.e() < gam_E_cut ) continue;

      setUserInfo( gam0 );
      setUserInfo( gam1 );
      UserInfo& info0 = dynamic_cast<UserInfo&>( gam0.userInfo() );
      UserInfo& info1 = dynamic_cast<UserInfo&>( gam1.userInfo() );
      info0.self( check_selfG(gam0) );
      info1.self( check_selfG(gam1) );

      // gen_Hepevt
      info.self( check_selfR(pi0) );
      Particle tmp_pi0( gam0.p() + gam1.p(), Ptype(PI0_LUND) );
      info.m_org( tmp_pi0.mass() );
      
      if( fl_message ) std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt
				 << " : Pi0 cand : "
				 << "gam0 ( "
				 << "gen-ID = " << std::setw(3) << std::right << gam0.genHepevt().get_ID() << ", "
				 << "self = "   << info0.self()              << "), "
				 << "gam1 ( "
				 << "gen-ID = " << std::setw(3) << std::right << gam1.genHepevt().get_ID() << ", "
				 << "self = "   << info1.self()              << ") "
				 << "-> "
				 << "pi0 ( "
				 << "self = " << info.self() << ") "
				 << std::endl;
      
      pi0_list.push_back( pi0 );
    }
    if( fl_message ) std::cout << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
			       << "                          # of pi0 candidates : " << pi0_list.size()                                                                        << std::endl
			       << "PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
