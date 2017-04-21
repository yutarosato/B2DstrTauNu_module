#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  void DSTRTAUNU::Rec_Dstr( std::vector<Particle>& dstr_list,
			    std::vector<Particle>& d_list, // charged D,  neutral D
			    std::vector<Particle>& p_list, // charged pi, neutral pi, gamma
			    const bool fl_message
			    )
  {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( d_list.size()==0 ) return;
    if( p_list.size()==0 ) return;
    int cnt = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message ) std::cout << "D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*" << std::endl;
    for( std::vector<Particle>::iterator d = d_list.begin(); d != d_list.end(); d++ ){
      for( std::vector<Particle>::iterator p = p_list.begin(); p != p_list.end(); p++ ){
	cnt++;
	UserInfo& info_D  = dynamic_cast<UserInfo&>( d->userInfo()  );
	UserInfo& info_p  = dynamic_cast<UserInfo&>( p->userInfo()  );
	if( fl_message ){
	  std::cout << "          CNTID = "   << std::setw(3) << std::right << cnt
		    << " : D_LUND = "
		    << std::setw(4) << std::right << d->lund() << "("
		    << "mode="  << std::setw(4) << std::right << info_D.rec_mode() << ","
		    << "cntID=" << std::setw(3) << std::right << info_D.cntid()    << ","
		    << "self="  << std::setw(3) << std::right << info_D.self()     << "), ";
	  rec_message( p );
	}

	int chg = int( d->charge() + p->charge() ); // total charge

	// charge check
	if( abs(chg) > 1 ){
	  if( fl_message ) std::cout << " -> not allowed charge : " << chg << std::endl;
	  continue;
	}

	// mode check
	if( d->charge() && p->lund()!=PI0_LUND ){ // NOTE : D+ is combined only with pi0
	                                          // D*+ -> D+ gamma (Low Br), D*0 -> D+ pi- (Not allowed)
	  if( fl_message ) std::cout << " -> not allowed mode : "
				     << d->charge() << " * " << p->charge()
				     << std::endl;
	  continue;
	}
	
	// flavor-charge check
	if( info_D.flavor() && info_D.flavor()*p->charge()==-1 ){
	  if( fl_message ) std::cout << " -> not allowed flavor-charge : "
				     << info_D.flavor() << " * " << p->charge()
				     << std::endl;
	  continue;
	}
	
	// Check duplication
	if( !check_dupli_daughter(d,p) ){
	  if( fl_message ) std::cout << " -> duplication " << std::endl;
	  continue;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	int dstr_lund = 0;
	int dstr_flavor = 0;
	if( info_D.flavor()==1 && chg!=0 ){ // D*+
	  dstr_lund   = Dstrp_LUND;
	  dstr_flavor = info_D.flavor();
	}else if( info_D.flavor()==1 && chg==0 ){ // D*0
	  dstr_lund   = Dstr0_LUND;
	  dstr_flavor = info_D.flavor();
	}else if( info_D.flavor()==-1 && chg!=0){ // D*-
	  dstr_lund   = Dstrm_LUND;
	  dstr_flavor = info_D.flavor();
	}else if( info_D.flavor()==-1 && chg==0){ // anti-D*0
	  dstr_lund   = antiDstr0_LUND;
	  dstr_flavor = info_D.flavor();
	}else{ // unknown D flavor
	  if( chg==1 ){ // D*+
	    dstr_lund   = Dstrp_LUND;
	    dstr_flavor = 1;
	  }else if( chg==-1 ){ // D*-
	    dstr_lund   = Dstrm_LUND;
	    dstr_flavor = -1;
	  }else{ // unknown case (D*0)
	    dstr_lund = Dstr0_LUND;
	    dstr_flavor = 0;
	  }
	}
	int rec_mode_dstr = p->lund();
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
	
	// reconstruction
	Particle Dstr( d->p() + p->p(), Ptype(dstr_lund) );
	
	Dstr.relation().append( *d );
	Dstr.relation().append( *p );
	setUserInfo( Dstr );
	UserInfo& info = dynamic_cast<UserInfo&>( Dstr.userInfo() );
	info.rec_mode( rec_mode_dstr );
	info.cntid   ( cnt           );
	info.flavor  ( dstr_flavor   );
	info.m_org   ( Dstr.mass()   );
	double delta_m = Dstr.mass()-info_D.m_org();
	if( delta_m > delta_m_cut ){
	  if( fl_message ) std::cout << " -> delta-m cut : " << delta_m << std::endl;
	  continue;
	}
	info.dm      ( delta_m );
	//info.self    ( check_selfR(Dstr) );
	info.self    ( check_selfR2(Dstr) );

	HepLorentzVector dstr_4Vcm = Dstr.p();
	dstr_4Vcm.boost( cmboost );
	info.Vcm( dstr_4Vcm );

	if( fl_message ) std::cout << " -> push_back : "
				   << "LUND="     << Dstr.lund()     << ", "
				   << "rec_mode=" << info.rec_mode() << ", "
				   << "d_flavor=" << info.flavor()   << ", "
				   << "self=" << info.self()
				   << ")"  << std::endl;
	dstr_list.push_back( Dstr );

      }
    }
    if( fl_message ) std::cout << "D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*D*" << std::endl;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
