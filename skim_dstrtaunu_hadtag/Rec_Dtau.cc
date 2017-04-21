#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Rec_Dtau( std::vector<Particle>& dtau_list,
			    std::vector<Particle>& d_list, // D(*)
			    std::vector<Particle>& lep_list,
			    std::vector<Particle>& gamma_list,
			    const bool fl_message
			    )
  {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( d_list.size()  ==0 ) return;
    if( lep_list.size()==0 ) return;
    int cnt = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message ) std::cout << "D(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tau" << std::endl;
    for( std::vector<Particle>::iterator d = d_list.begin(); d != d_list.end(); d++ ){
      for( std::vector<Particle>::iterator lep = lep_list.begin(); lep != lep_list.end(); lep++ ){
	cnt++;
	UserInfo& info_D  = dynamic_cast<UserInfo&>( d->userInfo()  );
	if( fl_message ) std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt
				   << " : D(*)_LUND = "
				   << std::setw(4) << std::right << d->lund() << "("
				   << std::setw(2) << std::right << info_D.rec_mode() << ","
				   << std::setw(3) << std::right << info_D.cntid()    << ")"
				   << ", lep_LUND = "
				   << std::setw(4) << std::right << lep->lund() << "("
				   << std::setw(2) << std::right << lep->mdstCharged().get_ID() << ")";

	int chg  = int( d->charge() + lep->charge() ); // total charge
	
	// Charge check
	if( abs(chg)>1 ){
	  if( fl_message ) std::cout << " -> not allowed charge : "
				     << d->charge() << " * " << lep->charge()
				     << std::endl;
	  continue;
	}
	// flavor-charge check
	if( info_D.flavor() && info_D.flavor() == lep->charge() ){
	  if( fl_message ) std::cout << " -> not allowed flavor-charge : "
				     << info_D.flavor() << " * " << lep->charge()
				     << std::endl;
	  continue;
	}
	
	// Check duplication
	if( !check_dupli_daughter(d,lep) ){
	  if( fl_message ) std::cout << " -> duplication " << std::endl;
	  continue;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int b_lund = 0;
	int b_flavor = 0;

	if( info_D.flavor()==1 && chg!=0 ){ // B-
	  b_lund   = Bminus_LUND;
	  b_flavor = -1;
	}else if( info_D.flavor()==1 && chg==0 ){ // anti-B0
	  b_lund   = antiB0_LUND;
	  b_flavor = -1;
	}else if( info_D.flavor()==-1 && chg!=0){ // B+
	  b_lund   = Bplus_LUND;
	  b_flavor = 1;
	}else if( info_D.flavor()==-1 && chg==0){ // B0
	  b_lund   = B0_LUND;
	  b_flavor = 1;
	}else{ // unknown D(*) flavor
	  if( chg==1 ){ // B+
	    b_lund   = Bplus_LUND;
	    b_flavor = 1;
	  }else if( chg==-1 ){ // B-
	    b_lund   = Bminus_LUND;
	    b_flavor = -1;
	  }else{ // unknown case (B0)
	    b_lund   = B0_LUND;
	    b_flavor = 0; // flavor can be determined from lepton-charge, however zero is assigned in this analysis.
	  }
	}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int rec_mode_ll = 0;
	if     ( abs(lep->lund())==Electron_LUND ) rec_mode_ll = 10;
	else if( abs(lep->lund())==MUminus_LUND  ) rec_mode_ll =  1;
	else   std::cerr << "[ABORT] Wrong lepton flavor : " << lep->lund() << std::endl, abort();
	
	// reconstruction
	Particle Dtau( d->p() + lep->p(), Ptype(b_lund) ); // dummy LUND code

	Dtau.relation().append( *d   );
	Dtau.relation().append( *lep );
	setUserInfo( Dtau );
	UserInfo& info = dynamic_cast<UserInfo&>( Dtau.userInfo() );
	info.flavor  ( info_D.flavor() );
	info.rec_mode( rec_mode_ll     );
	info.cntid   ( cnt             );
	info.flavor  ( b_flavor        );
	add_bremsstrahlung( Dtau, gamma_list );
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	/* not improved @ 20140414
	kvertexfitter kvf;
	unsigned err = 1;
	if( info_D.kf_ndf()>=0 ){
	  for( int i=1; i<Dtau.nChildren(); i++ ) addTrack2fit( kvf, Dtau.child(i) );
	  kvf.initialVertex( d->momentum().decayVertex() );
	  kvf.knownVertex();
	  err = kvf.fit();
	  if( err==0 ){ // success.
	    HepLorentzVector Dtau_mom = Dtau.child(0).p();
	    for( int i=1; i<Dtau.nChildren(); i++ ) Dtau_mom += kvf.momentum(i-1);
	    Dtau.momentum().momentum( Dtau_mom );
	    info.kf_cl   ( kvf.cl()    );
	    info.kf_chisq( kvf.chisq() );
	    info.kf_ndf  ( kvf.dgf()   );
	  }
	}
	*/
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	HepLorentzVector Dtau_4Vcm = Dtau.p();
	Dtau_4Vcm.boost( cmboost );
	info.Vcm( Dtau_4Vcm );

	HepLorentzVector Dtau_ctrl_4Vcm = Dtau.p();
	if( (abs(d->lund())==Dstr0_LUND || abs(d->lund())==Dstrp_LUND) ){
	  HepLorentzVector slow_pi_4Vcm = d->p() - d->child(0).p();
	  slow_pi_4Vcm.boost( cmboost );
	  Dtau_ctrl_4Vcm -= slow_pi_4Vcm;
	}

	double Bmass = ( chg ? PDG_BplusMass : PDG_B0Mass );
	double cosBDl = ( 2.0*eb*Dtau_4Vcm.e() - Bmass*Bmass-Dtau.mass()*Dtau.mass() ) / ( 2.0*sqrt(eb*eb-Bmass*Bmass)*Dtau_4Vcm.vect().mag() );
	info.cos( cosBDl );

	double cosBDl_ctrl = ( 2.0*eb*Dtau_ctrl_4Vcm.e() - Bmass*Bmass-Dtau_ctrl_4Vcm.mag()*Dtau_ctrl_4Vcm.mag() ) / ( 2.0*sqrt(eb*eb-Bmass*Bmass)*Dtau_ctrl_4Vcm.vect().mag() );
	info.cos_ctrl( cosBDl_ctrl ); // cos(theta) calculated from D and lepton in D*l mode (missing slow pion)

	if( fl_message ) std::cout << " -> push_back !" << std::endl;
	dtau_list.push_back( Dtau );

      }
    }
    if( fl_message ) std::cout << "D(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tauD(*)tau" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
