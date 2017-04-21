#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Mu_cand( std::vector<Particle>& mu_list, double pid_cut, double mom_cut, const bool fl_message )
  {
    if( fl_message ) std::cout << "MUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMU" << std::endl;
    Mdst_charged_Manager& chgMgr = Mdst_charged_Manager::get_manager();
    for( Mdst_charged_Manager::iterator c = chgMgr.begin();
	 c != chgMgr.end();
	 c++ )
      {
	// reduce beam B.G.
	double dr = correct_dr( (*c), m_IP, MUON_CODE );
	double dz = correct_dz( (*c), m_IP, MUON_CODE );
	if( fabs(dr) > Dr_cut ) continue;
	if( fabs(dz) > Dz_cut ) continue;

	// muid
	Muid_mdst muid( *c );
	if( muid.Chi_2()           <= 0       ) continue;
	if( muid.Muon_likelihood() <= pid_cut ) continue;

	// make particle
	Particle muon( (*c), Ptype(c->charge() > 0 ? MUplus_LUND : MUminus_LUND) );

	if( muon.ptot() < mom_cut ) continue;
		
	setUserInfo( muon );
	UserInfo& info = dynamic_cast<UserInfo&>(muon.userInfo());
	info.dr( fabs(dr) );
	info.dz( fabs(dz) );
	info.muonLikelihood( muid.Muon_likelihood() );

	// gen_hepevt
	info.self( check_selfF(muon)         );
	info.id      ( check_idF(muon)       );
	info.motherid( check_motheridF(muon) );

	// eid
	eid e( *c );
	double eidProb = e.prob( 3,-1,5 );
	info.eidProb( eidProb );

	// atckpi
	double atcKPI = selKPI.prob( *c );
	info.selKPI( atcKPI );

	HepLorentzVector muon_4Vcm = muon.p();
	muon_4Vcm.boost( cmboost );
	info.Vcm( muon_4Vcm );
	
	if( fl_message ) std::cout << "                          Muon cand : "
				   << "gen-ID = "         << std::setw( 3) << std::right << muon.genHepevt().get_ID()   << ", "
				   << "mdstCharged-ID = " << std::setw( 3) << std::right << muon.mdstCharged().get_ID() << ", "
				   << "gen-LUND = "       << std::setw( 6) << std::right << muon.genHepevt().idhep()    << ", "
				   << "mu-ID = "          << std::setw(15) << std::right << info.muonLikelihood()       << " : "
				   << "self = "           << std::setw( 3) << std::right << info.self()
				   << std::endl;
	
	mu_list.push_back( muon );
	
      }
    if( fl_message ) std::cout << "------------------------------------------------------------------------------" << std::endl
			       << "                          # of muon candidates : " << mu_list.size()            << std::endl
			       << "MUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMUMU" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
