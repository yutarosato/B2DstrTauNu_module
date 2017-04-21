#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Ks_cand( std::vector<Particle>& ks_list, int flag_kfitter, const bool fl_message )
  {
    int cnt = 0;
    if( fl_message ) std::cout << "KSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKS" << std::endl;
    Mdst_vee2_Manager &vee2_mgr = Mdst_vee2_Manager::get_manager();
    for( std::vector<Mdst_vee2>::iterator i = vee2_mgr.begin();
	 i != vee2_mgr.end(); i++ )
      {
	cnt++;
	Mdst_vee2 &vee2 = *i;
	if ( vee2.kind() != 1 ) continue; // kind: 1(Ks), 2(Lambda), 3(Lambda_bar), 4(gamma->ee)
	
	// FindKs
	FindKs findks;
	findks.candidates( vee2, m_IP );
	if( findks.goodKs() != 1 ) continue; // 1 = good Ks
	
	// make particle
	Particle ks( vee2 );
	setUserInfo( ks );
	UserInfo& info = dynamic_cast<UserInfo&>( ks.userInfo() );
	info.cntid( cnt );

	kvertexfitter kf;
	addTrack2fit( kf, ks.child(0) );
	addTrack2fit( kf, ks.child(1) );
	unsigned err = kf.fit();
	if( err ) continue;
	makeMother( kf, ks );


	HepLorentzVector ks_4Vcm = ks.p();
	ks_4Vcm.boost( cmboost );
	info.Vcm( ks_4Vcm );

	Particle& pip = ks.child( 0 );
	setUserInfo( pip );
	UserInfo& infop = dynamic_cast<UserInfo&>( pip.userInfo() );
	infop.self  ( check_selfF(pip) );
	infop.selKPI( selKPI.prob(pip.mdstCharged())  );
	infop.dr    ( fabs(correct_dr(pip.mdstCharged(), m_IP, PION_CODE)) );
	infop.dz    ( fabs(correct_dz(pip.mdstCharged(), m_IP, PION_CODE)) );

	Particle& pim = ks.child(1);
	setUserInfo( pim );
	UserInfo& infom = dynamic_cast<UserInfo&>( pim.userInfo() );
	infom.self  ( check_selfF(pim) );
	infom.selKPI( selKPI.prob(pim.mdstCharged())  );
	infom.dr    ( fabs(correct_dr(pim.mdstCharged(), m_IP, PION_CODE)) );
	infom.dz    ( fabs(correct_dz(pim.mdstCharged(), m_IP, PION_CODE)) );

	// gen_Hepevt
	info.self( check_selfR(ks) );

	if( fl_message ) std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt
				   << " : Ks cand : "
				   << "pi+ ( "
				   << "gen-ID = " << std::setw(3) << std::right << pip.genHepevt().get_ID() << ", "
				   << "self = "   << infop.self()              << "), "
				   << "pi- ( "
				   << "gen-ID = " << std::setw(3) << std::right << pim.genHepevt().get_ID() << ", "
				   << "self = "   << infom.self()              << ") "
				   << "-> "
				   << "Ks ( "
				   << "self = " << info.self() << ") "
				   << std::endl;
	
	ks_list.push_back( ks );
      }
    if( fl_message ) std::cout << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
			       << "                          # of Ks candidates : " << ks_list.size()                                                                          << std::endl
			       << "KSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKSKS" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
