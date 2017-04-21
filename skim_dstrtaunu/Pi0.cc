#include "belle.h"
#include "DSTRTAUNU.h"

#include <vector>
#include <iostream>
#include <stdlib.h>

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Pi0_cand( std::vector<Particle>& pi0_list, double masscut_L, double masscut_H, double cos_cut, double gam_e_cut_L, double gam_e_cut_H, double pi0_e_cut, const bool fl_message )
  {
    int cnt = 0;
    if( fl_message ) std::cout << "PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0" << std::endl;
    Mdst_pi0_Manager& pi0Mgr = Mdst_pi0_Manager::get_manager();
    for( Mdst_pi0_Manager::iterator p = pi0Mgr.begin(); p != pi0Mgr.end(); p++ ){
      cnt++;
      // make particle
      Particle pi0( (*p), true );
      setGammaError( pi0.child(0), m_IP, m_IP_err );
      setGammaError( pi0.child(1), m_IP, m_IP_err );
      kmassvertexfitter kf;
      addTrack2fit( kf, pi0.child(0) );
      addTrack2fit( kf, pi0.child(1) );
      kf.invariantMass( PDG_PI0Mass );
      
      unsigned err = kf.fit();
      if( err ) continue;
      makeMother( kf, pi0 );

      if( pi0.e() < pi0_e_cut ) continue;

      double morg = (pi0.child(0).p() + pi0.child(1).p()).m();
      if( !(morg - PDG_PI0Mass < masscut_H && morg - PDG_PI0Mass >= masscut_L) ) continue;
      setUserInfo( pi0 );
      UserInfo& info = dynamic_cast<UserInfo&>(pi0.userInfo());
      info.cntid( cnt );
      
      Particle& gam0 = pi0.child(0);
      Particle& gam1 = pi0.child(1);

      if( gam0.e() < gam_e_cut_L || gam1.e() < gam_e_cut_L ) continue;
      if( gam0.e() < gam_e_cut_H && gam1.e() < gam_e_cut_H ) continue;
      
      double cos_2gam = cos2track( pi0.child(0), pi0.child(1) );
      if( cos_2gam < cos_cut ) continue;
      info.cos ( cos_2gam );

      setUserInfo( gam0 );
      setUserInfo( gam1 );
      UserInfo& info0 = dynamic_cast<UserInfo&>( gam0.userInfo() );
      UserInfo& info1 = dynamic_cast<UserInfo&>( gam1.userInfo() );
      info0.self( check_selfG(gam0) );
      info1.self( check_selfG(gam1) );

      // gen_Hepevt
      info.self( check_selfPi0(pi0) );
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
				 << "self = " << info.self() << ", "
				 << "cos = "  << info.cos ()
				 << ") "
				 << std::endl;
      
      pi0_list.push_back( pi0 );

      // for pi0 study
      /*
      pi0_dist->column( "pi0e",   pi0.e()        );
      pi0_dist->column( "pi0m",   tmp_pi0.mass() );
      pi0_dist->column( "ghe",    gam0.e() > gam1.e() ? gam0.e() : gam1.e() );
      pi0_dist->column( "gle",    gam0.e() > gam1.e() ? gam1.e() : gam0.e() );
      pi0_dist->column( "cos",    info.cos()     );
      pi0_dist->column( "self",   info.self()    );
      pi0_dist->column( "selfg1", info0.self()   );
      pi0_dist->column( "selfg2", info1.self()   );
      pi0_dist->dumpData();
      */
    }
    if( fl_message ) std::cout << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
			       << "                          # of pi0 candidates : " << pi0_list.size()                                                                        << std::endl
			       << "PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0PI0" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
