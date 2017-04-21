#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Kaon_cand( std::vector<Particle>& kaon_list, const bool fl_message )
  {
    if( fl_message ) std::cout << "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK" << std::endl;
    Mdst_charged_Manager& chgMgr = Mdst_charged_Manager::get_manager();
    for( Mdst_charged_Manager::iterator c = chgMgr.begin();
	c != chgMgr.end();
	c++ )
      {
	// reduce beam B.G.
	double dr = correct_dr( (*c), m_IP, KAON_CODE );
	double dz = correct_dz( (*c), m_IP, KAON_CODE );
	if( fabs(dr) > Dr_cut ) continue;
	if( fabs(dz) > Dz_cut ) continue;

	// atckpi
	double atcKPI = selKPI.prob( *c );
	if( atcKPI <= selKPI_cut ) continue;
	
	// make particle
	Particle kaon( (*c), Ptype(c->charge() > 0 ? Kplus_LUND : Kminus_LUND) );
	setUserInfo( kaon );
	UserInfo& info = dynamic_cast<UserInfo&>( kaon.userInfo() );
	info.dr    ( fabs(dr) );
	info.dz    ( fabs(dz) );
	info.selKPI( atcKPI   );

	// gen_hepevt
	info.self    ( check_selfF(kaon)     );
	info.id      ( check_idF(kaon)       );
	info.motherid( check_motheridF(kaon) );
	
	if( fl_message ) std::cout << "                          Kaon cand : "
				   << "gen-ID = "         << std::setw( 3) << std::right << kaon.genHepevt().get_ID()   << ", "
				   << "mdstCharged-ID = " << std::setw( 3) << std::right << kaon.mdstCharged().get_ID() << ", "
				   << "gen-LUND = "       << std::setw( 6) << std::right << kaon.genHepevt().idhep()    << ", "
				   << "K-ID = "           << std::setw(15) << std::right << info.selKPI()               << " : "
				   << "self = "           << std::setw( 3) << std::right << info.self()
				   << std::endl;
	
	kaon_list.push_back( kaon );
      }
    if( fl_message ) std::cout << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
			       << "                          # of charged K candidates : " << kaon_list.size()                                                                 << std::endl
			       << "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
