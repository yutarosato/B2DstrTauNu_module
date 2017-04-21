#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Rec_CC( std::vector<Particle>& cc_list,
			  std::vector<Particle>& lep_list,
			  std::vector<Particle>& gamma_list,
			  const bool fl_message
			  )
  {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( lep_list.size()==0 ) return;
    int cnt = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message ) std::cout << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << std::endl;
    for( std::vector<Particle>::iterator p = lep_list.begin(); p!=lep_list.end(); p++ ){
      for(std::vector<Particle>::iterator m = lep_list.begin(); m!=lep_list.end(); m++ ){
	if( fl_message ) std::cout << "Lplus_LUND = "
				   << std::setw(4) << std::right << p->lund() << "("
				   << std::setw(2) << std::right << p->mdstCharged().get_ID() << ")"
				   << ", Lminus_LUND = "
				   << std::setw(4) << std::right << m->lund() << "("
				   << std::setw(2) << std::right << m->mdstCharged().get_ID() << ")";
	
	
	// charge check
	if( !(p->charge()>0 && m->charge()<0) ){
	  if( fl_message ) std::cout << " -> not allowed charge : " << p->charge() << " * " << m->charge() << std::endl;
	  continue;
	}
	// flavor check
	if( !((p->lund()*m->lund()==-121) || (p->lund()*m->lund()==-169)) ){
	  if( fl_message ) std::cout << " -> not allowed ee or mumu combination : " << p->lund() << " * " << m->lund() << std::endl;
	  continue;
	}
	
	// reconstruction from lepton pair
	Particle cc( p->p() + m->p(), Ptype(JPsi_LUND) );
	
	// mass cut
	//if( cc.mass() < pseudo_cc_M_cut ){
	//if( fl_message ) std::cout << " -> not allowed by masscut : " << pseudo_cc.mass() << std::endl;
	//continue;
	//}
	
	cc.relation().append( *p );
	cc.relation().append( *m );
	
	// userinfo
	setUserInfo( cc );
	UserInfo& info = dynamic_cast<UserInfo&>( cc.userInfo() );
	int rec_mode_ll = 0;
	if     ( abs(p->lund())==Electron_LUND ) rec_mode_ll = 10;
	else if( abs(p->lund())==MUminus_LUND  ) rec_mode_ll =  1;
	info.rec_mode( rec_mode_ll );
	
	HepLorentzVector cc_4Vcm = cc.p();
	cc_4Vcm.boost( cmboost );
	info.Vcm( cc_4Vcm );
	
	info.m_org( cc.mass() ); // original mass before correction of radiation gamma
	add_bremsstrahlung( cc, gamma_list );
	
	cc_list.push_back( cc );
	if( fl_message ) std::cout << std::endl;	    
      }
    }
    if( fl_message ) std::cout << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
