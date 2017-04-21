#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Tracks_cand( std::vector<Particle>& trk_list, const bool fl_message )
  {
    if( fl_message ) std::cout << "TRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRK" << std::endl;
    Mdst_charged_Manager& chgMgr = Mdst_charged_Manager::get_manager();
    for( Mdst_charged_Manager::iterator c = chgMgr.begin();
	c != chgMgr.end();
	c++ )
      {
	// reduce beam B.G.
	double dr = correct_dr( (*c), m_IP, PION_CODE );
	double dz = correct_dz( (*c), m_IP, PION_CODE );
	
	// make particle
	Particle trk( (*c), Ptype(c->charge() > 0 ? PIplus_LUND : PIminus_LUND) );
	setUserInfo( trk );
	UserInfo& info = dynamic_cast<UserInfo&>( trk.userInfo() );
	info.dr   ( fabs(dr)          );
	info.dz   ( fabs(dz)          );

	// gen_hepevt
	info.self ( check_selfF(trk)  );

	HepLorentzVector trk_4Vcm = trk.p();
	trk_4Vcm.boost( cmboost );
	info.Vcm( trk_4Vcm );

	if( fl_message ) std::cout << "                         Track cand : "
				   << "gen-ID = "         << std::setw(3) << std::right << trk.genHepevt().get_ID()   << ", "
				   << "mdstCharged-ID = " << std::setw(3) << std::right << trk.mdstCharged().get_ID() << ", "
				   << "gen-LUND = "       << std::setw(6) << std::right << trk.genHepevt().idhep()    << ", "
				   << std::endl;
	trk_list.push_back( trk );
      }
    if( fl_message ) std::cout << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
			       << "                          # of tracks candidates : " << trk_list.size()                                                                     << std::endl
			       << "TRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRKTRK" << std::endl;
  }

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
