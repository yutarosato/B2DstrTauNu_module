#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Gamma_cand( std::vector<Particle>& gamma_list, const bool fl_message )
  {
    if( fl_message ) std::cout << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" << std::endl;
    
    Mdst_gamma_Manager& gamMgr = Mdst_gamma_Manager::get_manager();
    for( Mdst_gamma_Manager::iterator g = gamMgr.begin(); g!=gamMgr.end(); g++ ){
      // make particle
      Particle gamma( *g );
      gamma.relation().mdstEcl( g->ecl() );
      setUserInfo( gamma );
      UserInfo& info = dynamic_cast<UserInfo&>( gamma.userInfo() );

      // gen_Hepevt
      info.self( check_selfG(gamma) );
      
      HepLorentzVector gamma_4Vcm = gamma.p();
      gamma_4Vcm.boost( cmboost );
      info.Vcm( gamma_4Vcm );
      setGammaError(gamma, IpProfile::position(), IpProfile::position_err());
      
      if( fl_message ){
	std::cout << "                          Gamma cand : "
		  << "gam-ID = " << std::setw(3)  << std::right << gamma.mdstEcl().get_ID();
	if( gamma.genHepevt() ) std::cout << ", "
					  << " gen-ID = " << std::setw(3)  << std::right << gamma.genHepevt().get_ID() << ", "
					  << " E = "      << std::setw(10) << std::right << gamma.mdstEcl().energy()   << ", "
					  << "self = "    << info.self();
	std::cout << std::endl;
      }
      gamma_list.push_back( gamma );
    }
    if( fl_message ){
      std::cout << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
		<< "                          # of gamma candidates : " << gamma_list.size()                                                                    << std::endl
		<< "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" << std::endl;
      }
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
