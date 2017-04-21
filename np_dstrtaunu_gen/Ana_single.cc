#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif


  int DSTRTAUNU::Ana_single( std::vector<Particle>& dtau_list, 
			     int gen_b_decay_info[][40], 
			     int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
			     std::vector<Particle>& gam_list, std::vector<Particle>& lep_list,
			     BelleTuple* dist, const bool fl_dump, const bool fl_message )
  {
    int cnt_dump = 0;
    if( fl_message ) std::cout << "ANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANa" << std::endl;
    for( std::vector<Particle>::iterator dtau1 = dtau_list.begin(); dtau1 != dtau_list.end(); dtau1++ ){
      Ana_single_sub( dtau1, gen_b_decay_info, gen_d_mode_info, gen_tau_mode_info, gam_list, lep_list, dist, fl_dump, fl_message );
      cnt_dump++;
    }
    if( fl_message ) std::cout << "ANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANaANa" << std::endl;
    return cnt_dump;
  }



  void DSTRTAUNU::Ana_single_sub( std::vector<Particle>::iterator dtau1,
				  int gen_b_decay_info[][40], 
				  int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
				  std::vector<Particle>& gam_list, std::vector<Particle>& lep_list,
				  BelleTuple* dist, const bool fl_dump, const bool fl_message )
  {
      // ++++++++++++++++++++++++++++++++++++++++ [ number ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    double exprunNo = expNo*10000+runNo;
    if( runNo==0 ) exprunNo = expNo*10000+runNoforsigMC; // for signal MC
    
    dist->column( "exprun",    exprunNo   );
    dist->column( "event",     evtNo      );
    dist->column( "evt_type",  event_type );
    
    // +++++++++++++++++++++++++++++++++ [ make Particl object ] ++++++++++++++++++++++++++++++++++++++++++++++++++
    // [ check (DD/D*D/D*D*) ]
    // (if D*D, recD1 is D* and recD2 is D)
    int fl_dstr1  = 0; // 1(D*), 0(D)
    //int fl_switch = 0; // 1(DD*), other(0)
    if     ( abs(dtau1->child(0).lund())==Dstrp_LUND || abs(dtau1->child(0).lund())==Dstr0_LUND ) fl_dstr1 = 1;
    else if( abs(dtau1->child(0).lund())==Dplus_LUND || abs(dtau1->child(0).lund())==D0_LUND    ) fl_dstr1 = 0;
    else std::cerr << "[ABORT] Wrong D(*) LUND" << std::endl, abort();
    //if( fl_dstr1==0 && fl_dstr2==1 ) fl_switch = 1;
    Particle& Dtau1 = *dtau1;
    Particle& Dstr1 = Dtau1.child(0); // if Dtaunu mode, "Dstr1" is D (not D*)
    UserInfo& info_Dtau1 = dynamic_cast<UserInfo&>( Dtau1.userInfo() );
    UserInfo& info_Dstr1 = dynamic_cast<UserInfo&>( Dstr1.userInfo() );
    
    Particle& D1    = ( info_Dstr1.dm()==-1 ? Dtau1.child(0) : Dstr1.child(0) );
    Particle& lep1  = Dtau1.child(1);
    UserInfo& info_D1    = dynamic_cast<UserInfo&>( D1.userInfo()    );
    UserInfo& info_lep1  = dynamic_cast<UserInfo&>( lep1.userInfo()  );
    
    if( fl_message ) std::cout << "Dtau_CNTID(1st) = " << std::setw(3) << std::right << info_Dtau1.cntid();
    
    // ++++++++++++++++++++++++++++++++++++ [ rec-mode(BB) ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    int rm_bb = (Dtau1.charge() ? 4 : 3); // 3(B0B0bar) 4(B+B-)
    dist->column( "rm_bb", rm_bb );
    
    // ++++++++++++++++++++++++++++++++++++ [ rec-mode(B) ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    dist->column( "rm_l1",       info_Dtau1.rec_mode() );
    dist->column( "rm_d1",       info_D1.rec_mode()    );
    dist->column( "rm_d1lund",   D1.lund()             );
    dist->column( "rm_dst1lund", info_Dstr1.dm()==-1 ? 0 : Dstr1.lund() );
    dist->column( "dfl1",        info_Dtau1.flavor() );
    dist->column( "dtchg1",      Dtau1.charge()      );
    
    if( info_Dstr1.dm()==-1 ) dist->column( "rm_dst1",   0                     ); // B -> D  lnu
    else                      dist->column( "rm_dst1",   info_Dstr1.rec_mode() ); // B -> D* lnu

    // +++++++++++++++++++++++++++++++++++++ [ momentum ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    dist->column( "l1pcm",   info_lep1.Vcm().vect().mag() );
    dist->column( "l1p",     lep1.p().vect().mag()        );
    dist->column( "l1pt",    lep1.p().perp()              );
    dist->column( "l1pc",    lep1.p().cosTheta()          ); // @20141202
    
    if( info_Dstr1.dm()==-1 ){ // B -> D  lnu
      dist->column( "dst1pcm", -1 );
      dist->column( "dst1p",   -1 );
      dist->column( "acc1p",   -1 );
      dist->column( "acc1pc",  -1 ); // @ 20141202
      dist->column( "acc1m",   -1 ); // @ 20141105
    }else{ // B -> D* lnu
      dist->column( "dst1pcm", info_Dstr1.Vcm().vect().mag()   );
      dist->column( "dst1p",   Dstr1.p().vect().mag()          );
      dist->column( "acc1p",   Dstr1.child(1).p().vect().mag() );
      dist->column( "acc1pc",  Dstr1.child(1).p().cosTheta()   ); // @20141202
      UserInfo& info_accpi0 = dynamic_cast<UserInfo&>( Dstr1.child(1).userInfo() );
      dist->column( "acc1m",   info_accpi0.m_org()             ); // @ 20141105
    }
    
    dist->column( "d1pcm",   info_D1.Vcm().vect().mag() );
    dist->column( "d1p",     Dstr1.p().vect().mag()     );

    double thrust_angle;
    double r2;
    qq_suppress( Dtau1, thrust_angle, r2 );
    dist->column( "thrust",    thrust_angle );
    dist->column( "r2",        r2           );

    // ++++++++++++++++++++++++++++++++++++++ [ D(*)tau ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    dist->column( "cosdl1", info_Dtau1.cos()      );
    
    //dist->column( "kfl1cl",   info_Dtau1.kf_cl()    );
    //dist->column( "kfl1chi2", info_Dtau1.kf_chisq() );
    //dist->column( "kfl1ndf",  info_Dtau1.kf_ndf()   );
    // ++++++++++++++++++++++++++++ [D* and accomapany particle] ++++++++++++++++++++++++++++++++++++++++++++++++++

    int dst1self = 0;
    int acc1self = 0; // accompany particle (e.g. slow pion)
    int dst1org  = -1;
    int acc1org  = -1;
    int acc1mo   = 0;
    int acc1gmo  = 0;
    
    if( info_Dstr1.dm()==-1 ){ // B -> D lnu
      dist->column( "dm1",       -1 );
      dist->column( "dst1_morg", -1 );
      dist->column( "dst1_m",    -1 );
      dist->column( "kfs1cl",    -1 );
      dist->column( "kfs1chi2",  -1 );
      dist->column( "kfs1ndf",   -1 );
      dst1self = -1;
      acc1self = -1;
    }else{ // B -> D* lnu
      dist->column( "dm1",        info_Dstr1.dm()       );
      dist->column( "dst1_morg",  info_Dstr1.m_org()    );
      dist->column( "dst1_m",     Dstr1.mass()          );
      dist->column( "kfs1cl",     info_Dstr1.kf_cl()    );
      dist->column( "kfs1chi2",   info_Dstr1.kf_chisq() );
      dist->column( "kfs1ndf",    info_Dstr1.kf_ndf()   );
      dst1self = info_Dstr1.self();
      dst1org = which_B( Dstr1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
      UserInfo& info_acc1 = dynamic_cast<UserInfo&>( Dstr1.child(1).userInfo() );
      acc1self = info_acc1.self();
      acc1org = which_B( Dstr1.child(1).genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
      if( Dstr1.child(1).genHepevt() && Dstr1.child(1).genHepevt().mother()                                                 ) acc1mo  = Dstr1.child(1).genHepevt().mother().idhep();
      if( Dstr1.child(1).genHepevt() && Dstr1.child(1).genHepevt().mother() && Dstr1.child(1).genHepevt().mother().mother() ) acc1gmo = Dstr1.child(1).genHepevt().mother().mother().idhep();
    }
    
    dist->column( "dst1self",  dst1self );
    dist->column( "dst1org",   dst1org  );
    dist->column( "acc1self",  acc1self );
    dist->column( "acc1org",   acc1org  );
    dist->column( "acc1mo",    acc1mo   );
    dist->column( "acc1gmo",   acc1gmo  );
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++

    if( info_Dstr1.rec_mode()==PI0_LUND ){
      UserInfo& info_pi0 = dynamic_cast<UserInfo&>( Dstr1.child(1).userInfo() );
	  dist->column( "dst1_pi0cos", info_pi0.cos() );
    }else{
      dist->column( "dst1_pi0cos", 2 );
    }
    
    // +++++++++++++++++++++++++++++++++++++++++ [D] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    dist->column( "d1_morg", info_D1.m_org() );
    dist->column( "d1_m",    D1.mass()       );
    
    dist->column( "kfd1cl",   info_D1.kf_cl()    );
    dist->column( "kfd1chi2", info_D1.kf_chisq() );
    dist->column( "kfd1ndf",  info_D1.kf_ndf()   );
    
    int d1self = info_D1.self();
    dist->column( "d1self", d1self );

    if( D1.genHepevt() ) dist->column( "d1moid", D1.genHepevt().mother().idhep() );
    else                 dist->column( "d1moid", 0                               );
    
    int d1org = which_B( D1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
    
    dist->column( "d1org", d1org );
    
    // ++++++++++++++++++++++++++++++++++++++ [lepton] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int l1org = which_B( lep1.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
    int l1self;
    
    if( lep1.lund()==info_lep1.id() ){ // true lepton
      if     ( (gen_b_decay_info[0][1]==1 || gen_b_decay_info[0][1]==2) && lep1.genHepevt().get_ID()==gen_b_decay_info[0][ 2] ) l1self = 2; // tag-side true lepton
      else if( (gen_b_decay_info[1][1]==1 || gen_b_decay_info[1][1]==2) && lep1.genHepevt().get_ID()==gen_b_decay_info[1][ 2] ) l1self = 2; // tag-side true lepton
      else if(                               gen_b_decay_info[0][1]==3  && lep1.genHepevt().get_ID()==gen_b_decay_info[0][22] ) l1self = 1; // sig-side true lepton
      else if(                               gen_b_decay_info[1][1]==3  && lep1.genHepevt().get_ID()==gen_b_decay_info[1][22] ) l1self = 1; // sig-side true lepton
      else if( info_D1.self()==1 ){ // true D
	if( d1org==l1org )                                                                                                      l1self = 3; // random   true lepton from same     D side
	else                                                                                                                    l1self = 5; // random   true lepton from opposite D side
      }else                                                                                                                     l1self = 7; // random   true lepton (fake D)
    }else{ // fake lepton
      if( info_D1.self()==1 ){ // true D
	if( d1org==l1org ) l1self = -3; // same D side
	else               l1self = -5; // opposite D side
      }else                l1self = -7; // fake D
	}
    
    dist->column( "l1self",   l1self               );
    dist->column( "l1selfid", info_lep1.id()       );
    dist->column( "l1moid",   info_lep1.motherid() );
    dist->column( "l1org",    l1org                );
    
    if     ( info_Dtau1.rec_mode() == 10 ) dist->column( "l1pid",    info_lep1.eidProb()        ); // 10(e)
    else if( info_Dtau1.rec_mode() ==  1 ) dist->column( "l1pid",    info_lep1.muonLikelihood() ); //  1(mu)

    // +++++++++++++++++++++++++[ remaining lepton @20141202 ] +++++++++++++++++++++++++;

    int nleps1 = 0; //   same   sign charged lepton with PID>0.1
    int nlepo1 = 0; // opposite sign charged lepton with PID>0.1
    int nleps9 = 0; //   same   sign charged lepton with PID>0.9
    int nlepo9 = 0; // opposite sign charged lepton with PID>0.9
    for( std::vector<Particle>::iterator remlep = lep_list.begin(); remlep != lep_list.end(); remlep++ ){
      UserInfo& info_remlep = dynamic_cast<UserInfo&>( remlep->userInfo()  );
      remlep->charge();
      double pid = 0;
      if     ( abs(remlep->lund())==Electron_LUND ) pid = info_remlep.eidProb();
      else if( abs(remlep->lund())==MUminus_LUND  ) pid = info_remlep.muonLikelihood();

      if( remlep->charge()==lep1.charge() ){
	if( pid > 0.1 ) nleps1++;
	if( pid > 0.9 ) nleps9++;
      }else{
	if( pid > 0.1 ) nlepo1++;
	if( pid > 0.9 ) nlepo9++;
      }
    }
    dist->column( "nls1", nleps1 );
    dist->column( "nls9", nleps9 );
    dist->column( "nlo1", nlepo1 );
    dist->column( "nlo9", nlepo9 );

    // +++++++++++++++++++++++++++++++ [D's daughter particles] +++++++++++++++++++++++++++++++++++++++++++++++++++
    dist->column( "d1nch", D1.nChildren() );

    for( int i=0; i<4; i++ ){
      int reclund =  -1;
      int genlund =  -1;
      int moid    =  -1;
      int gmoid   =  -1; // @20150123 
      int org     =  -1;
      double mom  =  -1; // @20141202
      double cos  = -10; // @20141202
      
      if( i < D1.nChildren() ){
	Particle& Ddaughter = D1.child(i);
	reclund = Ddaughter.lund();
	if( Ddaughter.genHepevt() ){
	  genlund = Ddaughter.genHepevt().idhep();
	  if( Ddaughter.genHepevt().mother() ) moid = Ddaughter.genHepevt().mother().idhep();
	  else                                 moid  = 0;
	  if( Ddaughter.genHepevt().mother() && Ddaughter.genHepevt().mother().mother() ) gmoid = Ddaughter.genHepevt().mother().mother().idhep();
	  else                                                                            gmoid = 0;
	  org = which_B( Ddaughter.genHepevt(), gen_b_decay_info[0][0], gen_b_decay_info[1][0] );
	}
	mom = Ddaughter.p().vect().mag();
	cos = Ddaughter.p().cosTheta();
      }
      if      ( i==0 ){ dist->column( "d1ch0recid", reclund ); dist->column( "d1ch0selfid", genlund ); dist->column( "d1ch0mo", moid ); dist->column( "d1ch0gmo", gmoid ); dist->column( "d1ch0org", org ); dist->column( "d1ch0p", mom ); dist->column( "d1ch0pc", cos );
      }else if( i==1 ){ dist->column( "d1ch1recid", reclund ); dist->column( "d1ch1selfid", genlund ); dist->column( "d1ch1mo", moid ); dist->column( "d1ch1gmo", gmoid ); dist->column( "d1ch1org", org ); dist->column( "d1ch1p", mom ); dist->column( "d1ch1pc", cos );
      }else if( i==2 ){ dist->column( "d1ch2recid", reclund ); dist->column( "d1ch2selfid", genlund ); dist->column( "d1ch2mo", moid ); dist->column( "d1ch2gmo", gmoid ); dist->column( "d1ch2org", org ); dist->column( "d1ch2p", mom ); dist->column( "d1ch2pc", cos );
      }else if( i==3 ){ dist->column( "d1ch3recid", reclund ); dist->column( "d1ch3selfid", genlund ); dist->column( "d1ch3mo", moid ); dist->column( "d1ch3gmo", gmoid ); dist->column( "d1ch3org", org ); dist->column( "d1ch3p", mom ); dist->column( "d1ch3pc", cos );
      }
    }
    // +++++++++++++++++++++++++++++++++++++++++++++++++++
    
    double k1pidmax = -1;
    double k1pidmin = 2;
    for( int i=0; i<D1.nChildren(); i++ ){
      if( abs(D1.child(i).lund())==Kplus_LUND ){
	UserInfo& info_k = dynamic_cast<UserInfo&>( D1.child(i).userInfo() );
	if( k1pidmax < info_k.selKPI() ) k1pidmax = info_k.selKPI();
	if( k1pidmin > info_k.selKPI() ) k1pidmin = info_k.selKPI();
      }
    }
    if( k1pidmax==-1 ) k1pidmax = 2;
    dist->column( "k1pidmin", k1pidmin );
    dist->column( "k1pidmax", k1pidmax );

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
    // additional pi0 information
    if( info_D1.rec_mode() > 999 ){
      for( int i=0; i<D1.nChildren(); i++ ){
	if( D1.child(i).lund()==PI0_LUND ){
	  UserInfo& info_pi0 = dynamic_cast<UserInfo&>( D1.child(i).userInfo() );
	  dist->column( "d1_pi0cos", info_pi0.cos()   );
	  dist->column( "pi01m",     info_pi0.m_org() );
	  dist->column( "pi01self",  info_pi0.self()  );
	}
      }
    }else{
      dist->column( "d1_pi0cos",  2 );
      dist->column( "pi01m",     -1 );
      dist->column( "pi01self",   1 );
    }

    //+++++++++++++++++++++++++++++++++++
    double d1_pi0_mass = -1;
    for( int i=0; i<D1.nChildren(); i++ ){
      if( D1.child(i).lund()==PI0_LUND ){
	UserInfo& info_pi0 = dynamic_cast<UserInfo&>( D1.child(i).userInfo() );

      }
    }


    // ++++++++++++++++++++++++++++++++ [bremsstrahlung study] ++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( mc==1 ){
      // <bremsstrahlung from generator inforamation>
      dist->column( "l1geg", radgam_energy(lep1.genHepevt()) );
      dist->column( "l1gng", radgam_cnt   (lep1.genHepevt()) );
      // +++++++++++++++++++++++++++++++	
      Hep3Vector lep1_mom_gen( lep1.genHepevt().PX(), lep1.genHepevt().PY(), lep1.genHepevt().PZ() );
      dist->column( "l1ag", lep1.p3().angle(lep1_mom_gen) ); // deviation angle between rec and gen
      // +++++++++++++++++++++++++++++++
      // <reconstructed bremsstrahlung gamma>
      dist->column( "l1gnr", info_lep1.radgam_n() );
      // +++++++++++++++++++++++++++++++
      // <truely reconstructed bremsstrahlung gamma>
      int    cnt_true_radgam1   = 0;
      double cnt_true_radgam1_e = 0;
      for( int i=0; i<info_lep1.radgam_n(); i++ ){
	if( Dtau1.child(2+i).genHepevt() && lep1.genHepevt().mother().get_ID() == Dtau1.child(2+i).genHepevt().mother().get_ID() ){
	  cnt_true_radgam1   += 1;
	  cnt_true_radgam1_e += Dtau1.child(2+i).mdstEcl().energy();
	}
      }
      dist->column( "l1gntr", cnt_true_radgam1   );
      dist->column( "l1getr", cnt_true_radgam1_e );
      // +++++++++++++++++++++++++++++++
      //<detected bremsstrahlung gamma>
      int    cnt_det_radgam1   = 0;
      double cnt_det_radgam1_e = 0;

      dist->column( "rgd1e", -10 );
      dist->column( "rgd1a", -10 );
      dist->column( "rgd1r",   0 );
      dist->column( "rgd2e", -10 );
      dist->column( "rgd2a", -10 );
      dist->column( "rgd2r",   0 );
      dist->column( "rgd3e", -10 );
      dist->column( "rgd3a", -10 );
      dist->column( "rgd3r",   0 );
      dist->column( "rgd4e", -10 );
      dist->column( "rgd4a", -10 );
      dist->column( "rgd4r",   0 );
      
      
      for( std::vector<Particle>::iterator gam = gam_list.begin(); gam != gam_list.end(); gam++ ){
	if( !gam->genHepevt() ) continue;
	// ++++++++++++++++++++++++++++++++++++
	double tmp_angle = -10;
	if( gam->genHepevt().mother().get_ID() == lep1.genHepevt().mother().get_ID() ){
	  cnt_det_radgam1   += 1;
	  cnt_det_radgam1_e += gam->mdstEcl().energy();
	  tmp_angle = lep1.p3().angle( gam->p3() );
	}
	// ++++++++++++++++++++++++++++++++++++
	int tmp_rec = 0;
	for( int i=0; i<info_lep1.radgam_n(); i++ ){
	  if( Dtau1.child(2+i).genHepevt() && Dtau1.child(2+i).genHepevt().get_ID() == gam->genHepevt().get_ID() ){
	    tmp_rec = 1;
	  }
	}
	// ++++++++++++++++++++++++++++++++++++
	if( gam->genHepevt().mother().get_ID() == lep1.genHepevt().mother().get_ID() ){
	  int cnt = cnt_det_radgam1;
	  
	  if( cnt==1 ){
	    dist->column( "rgd1e", gam->mdstEcl().energy() );
	    dist->column( "rgd1a", tmp_angle               );
	    dist->column( "rgd1r", tmp_rec                 );
	  }else if( cnt==2 ){
		dist->column( "rgd2e", gam->mdstEcl().energy() );
		dist->column( "rgd2a", tmp_angle               );
		dist->column( "rgd2r", tmp_rec                 );
	  }else if( cnt==3 ){
	    dist->column( "rgd3e", gam->mdstEcl().energy() );
	    dist->column( "rgd3a", tmp_angle               );
	    dist->column( "rgd3r", tmp_rec                 );
	  }else{
	    dist->column( "rgd4e", gam->mdstEcl().energy() );
	    dist->column( "rgd4a", tmp_angle               );
	    dist->column( "rgd4r", tmp_rec                 );
	  }
	}
      }

      dist->column( "l1ged", cnt_det_radgam1_e );
      dist->column( "l1gnd", cnt_det_radgam1   );
    }
    // +++++++++++++++++++++++++++++++
    // make more bremsstrahlung recovery
    dist->column( "cosdla1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.05, MUminus_LUND) );
    dist->column( "cosdlb1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.10, MUminus_LUND) );
    dist->column( "cosdlc1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.10              ) );
    dist->column( "cosdld1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.20, MUminus_LUND) );
    dist->column( "cosdle1",  cosBDl_with_more_bremsstrahlung       (Dtau1,        gam_list, 0.20              ) );

    // +++++++++++++++++++++++++++++++++ [ gen-rec matching ] +++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    int recfB_id = 0;
    int recsB_id = 0;
    int recgen   = 0;
    if( info_D1.self()==1 ){ // true D1
      if( d1org==1 ){
	recfB_id = 0;
	recsB_id = 1;
	recgen   = 1;
      }else if( d1org==2 ){
	recfB_id =  1;
	recsB_id =  0;
	recgen   = -1;
      }
    }else{ // fake D1
      recfB_id = 0;
      recsB_id = 1;
      recgen   = 0;
    }
    
    dist->column( "recgen",  recgen ); // 1(rec1=gen1,rec2=gen2), -1(rec1=gen2,rec2=gen1), 0(rec1=gen1,rec2=gen2)
    
    
    // +++++++++++++++++++++++++++++++++++ [ gen-mode(BB) ] +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int gm_dd   = make_flag_dd  ( gen_b_decay_info[0][13], gen_b_decay_info[1][13] );
    int gm_semi = make_flag_semi( gen_b_decay_info[0][ 1], gen_b_decay_info[1][ 1] );
    dist->column( "gm_dd",   gm_dd   );
    dist->column( "gm_semi", gm_semi );
    // +++++++++++++++++++++++++++++++++++ [ gen-mode(B) ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int nb1  = gen_b_decay_info[0][3];
    int nb2  = gen_b_decay_info[1][3];
    int nb1g = gen_b_decay_info[0][4];
    int nb2g = gen_b_decay_info[1][4];
    int semi1 = gen_b_decay_info[0][1];
    int semi2 = gen_b_decay_info[1][1];
    dist->column( "nb1",       nb1                     ); dist->column( "nb2",       nb2                     );
    dist->column( "nb1g",      nb1g                    ); dist->column( "nb2g",      nb2g                    );
    dist->column( "nb1pip",    gen_b_decay_info[0][29] ); dist->column( "nb2pip",    gen_b_decay_info[1][29] );
    dist->column( "nb1pi0",    gen_b_decay_info[0][30] ); dist->column( "nb2pi0",    gen_b_decay_info[1][30] );
    dist->column( "semi1",     semi1                   ); dist->column( "semi2",     semi2                   );
    dist->column( "nrootd1",   gen_b_decay_info[0][ 5] ); dist->column( "nrootd2",   gen_b_decay_info[1][ 5] );
    dist->column( "rootdf1",   gen_b_decay_info[0][ 6] ); dist->column( "rootdf2",   gen_b_decay_info[1][ 6] ); 
    dist->column( "rootds1",   gen_b_decay_info[0][ 8] ); dist->column( "rootds2",   gen_b_decay_info[1][ 8] ); // fixed @ 20150123
    dist->column( "nd1",       gen_b_decay_info[0][10] ); dist->column( "nd2",       gen_b_decay_info[1][10] );
    int gm_ddst1 = gen_b_decay_info[0][13];
    int gm_ddst2 = gen_b_decay_info[1][13];
    dist->column( "gm_ddst1",  gm_ddst1                ); dist->column( "gm_ddst2",  gm_ddst2                );
    dist->column( "dst1_acc",  gen_b_decay_info[0][14] ); dist->column( "dst2_acc",  gen_b_decay_info[1][14] );
    dist->column( "dst1_d",    gen_b_decay_info[0][27] ); dist->column( "dst2_d",    gen_b_decay_info[1][27] ); // 20141002
    dist->column( "cc1",       gen_b_decay_info[0][16] ); dist->column( "cc2",       gen_b_decay_info[1][16] );
    
    dist->column( "rootd1nc",  gen_b_decay_info[0][31] ); dist->column( "rootd2nc",  gen_b_decay_info[1][31] ); // 20141028
    dist->column( "fldstst1",  gen_b_decay_info[0][32] ); dist->column( "fldstst2",  gen_b_decay_info[1][32] ); // 20141028
    dist->column( "incl1gmc",  gen_b_decay_info[0][33] ); dist->column( "incl2gmc",  gen_b_decay_info[1][33] ); // 20141028
    dist->column( "nddst1pp",  gen_b_decay_info[0][34] ); dist->column( "nddst2pp",  gen_b_decay_info[1][34] ); // 20141028
    dist->column( "nddst1p0",  gen_b_decay_info[0][35] ); dist->column( "nddst2p0",  gen_b_decay_info[1][35] ); // 20141028
    dist->column( "flmcddst",  flag_DststMC );
    
    dist->column( "taulep1",   gen_b_decay_info[0][21] ); dist->column( "taulep2",   gen_b_decay_info[1][21] ); // 20140925
    dist->column( "ddst1_d",   gen_b_decay_info[0][23] ); dist->column( "ddst2_d",   gen_b_decay_info[1][23] ); // 20140703
    dist->column( "ddst1_acc", gen_b_decay_info[0][25] ); dist->column( "ddst2_acc", gen_b_decay_info[1][25] ); // 20140703
    
    // +++++++++++++++++++++++++++++++++++ [ gen-mode(D) ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    dist->column( "gm_df1hr",  gen_d_mode_info [0][0][0]  ); dist->column( "gm_df2hr",  gen_d_mode_info [1][0][0] );
    dist->column( "gm_ds1hr",  gen_d_mode_info [0][1][0]  ); dist->column( "gm_ds2hr",  gen_d_mode_info [1][1][0] );
    dist->column( "gm_df1hnr", gen_d_mode_info [0][0][1]  ); dist->column( "gm_df2hnr", gen_d_mode_info [1][0][1] );
    dist->column( "gm_ds1hnr", gen_d_mode_info [0][1][1]  ); dist->column( "gm_ds2hnr", gen_d_mode_info [1][1][1] );
    dist->column( "gm_df1lep", gen_d_mode_info [0][0][2]  ); dist->column( "gm_df2lep", gen_d_mode_info [1][0][2] );
    dist->column( "gm_ds1lep", gen_d_mode_info [0][1][2]  ); dist->column( "gm_ds2lep", gen_d_mode_info [1][1][2] );
    dist->column( "gm_df1nu",  gen_d_mode_info [0][0][3]  ); dist->column( "gm_df2nu",  gen_d_mode_info [1][0][3] );
    dist->column( "gm_ds1nu",  gen_d_mode_info [0][1][3]  ); dist->column( "gm_ds2nu",  gen_d_mode_info [1][1][3] );
    dist->column( "gm_df1gam", gen_d_mode_info [0][0][4]  ); dist->column( "gm_df2gam", gen_d_mode_info [1][0][4] );
    dist->column( "gm_ds1gam", gen_d_mode_info [0][1][4]  ); dist->column( "gm_ds2gam", gen_d_mode_info [1][1][4] );
    
    dist->column( "gm_t1hr",  gen_tau_mode_info [0][0]  ); dist->column( "gm_t2hr",  gen_tau_mode_info [1][0] );
    dist->column( "gm_t1hnr", gen_tau_mode_info [0][1]  ); dist->column( "gm_t2hnr", gen_tau_mode_info [1][1] );
    dist->column( "gm_t1lep", gen_tau_mode_info [0][2]  ); dist->column( "gm_t2lep", gen_tau_mode_info [1][2] );
    dist->column( "gm_t1nu",  gen_tau_mode_info [0][3]  ); dist->column( "gm_t2nu",  gen_tau_mode_info [1][3] );
    dist->column( "gm_t1gam", gen_tau_mode_info [0][4]  ); dist->column( "gm_t2gam", gen_tau_mode_info [1][4] );
    dist->column( "gm_t1pro", gen_tau_mode_info [0][5]  ); dist->column( "gm_t2pro", gen_tau_mode_info [1][5] );

    // ++++++++++++++++++++++++++++++++++++++++ [    D** mass(gen)     ] ++++++++++++++++++++++++++++++++++++++++++
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    if( gen_b_decay_info[0][13]==3 ) dist->column( "ddst1m", genMgr( Panther_ID(gen_b_decay_info[0][7]) ).M() );
    else                             dist->column( "ddst1m", -10                                              );
    if( gen_b_decay_info[1][13]==3 ) dist->column( "ddst2m", genMgr( Panther_ID(gen_b_decay_info[1][7]) ).M() );
    else                             dist->column( "ddst2m", -10                                              );
    
    // ++++++++++++++++++++++++++++++++++++++++ [kinematics for weight ] ++++++++++++++++++++++++++++++++++++++++++
    // for MC model correction, calculate q^2, p_l*, w, cos(theta_l) in B->D(**) l nu.
    for( int i=0; i<2; i++ ){
      double q2     = -10;
      double plep   = -10;
      double w      = -10;
      double coslep = -10;
      if( (gen_b_decay_info[i][1]==1 || gen_b_decay_info[i][1]==2) &&
	  (
	   ( gen_b_decay_info[i][13]==1 ) || // B -> D  l nu
	   ( gen_b_decay_info[i][13]==2 ) || // B -> D* l nu
	   ( gen_b_decay_info[i][13]==3 && gen_b_decay_info[i][32]==1 ) )
	  ){
	calKinematics( gen_b_decay_info[i][0], gen_b_decay_info[i][7], gen_b_decay_info[i][2], q2, plep, w, coslep );
      }
      if( fl_message ) std::cout << "Kinematics Calculation " << i << " : "
				 << "q2 = "     << q2   << ", "
				 << "plep = "   << plep << ", "
				 << "w = "      << w    << ", "
				 << "coslep = " << coslep
				 << std::endl;
      if( i==0 ){
	dist->column( "k1_q2",  q2     );
	dist->column( "k1_pl",  plep   );
	dist->column( "k1_w",   w      );
	dist->column( "k1_cos", coslep );
      }else{
	dist->column( "k2_q2",  q2     );
	dist->column( "k2_pl",  plep   );
	dist->column( "k2_w",   w      );
	dist->column( "k2_cos", coslep );
      }
    }

    // ++++++++++++++++++++++++++++++++++++++++ [ self ] ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    int old_self = 0;
    if     ( ( abs(event_type)==rm_bb ) && ( dst1self!=0&&d1self==1&&l1self>0 ) && ( l1self==1 ) ) old_self =  1; // sig
    else if( ( abs(event_type)==rm_bb ) && ( dst1self!=0&&d1self==1&&l1self>0 ) && ( l1self==2 ) ) old_self =  2; // tag
    else if( ( abs(event_type)==rm_bb )                                                          ) old_self =  0; // fake
    else if( ( abs(event_type)==rm_bb )                                                          ) old_self = -1; // B+B-
    else if( ( abs(event_type)!=rm_bb )                                                          ) old_self = -2; // qq
    dist->column( "old_self", old_self );

    int self = 0;
    if     ( ( abs(event_type)==rm_bb ) && ( (recgen==1 && nb1-nb1g==3 && semi1>0 && gm_ddst1==2) || (recgen==-1 && nb2-nb2g==3 && semi2>0 && gm_ddst2==2) ) && ( dst1self!=0&&d1self==1&&l1self>0 && l1self==1 ) ) self =  1; // sig
    else if( ( abs(event_type)==rm_bb ) && ( (recgen==1 && nb1-nb1g==3 && semi1>0 && gm_ddst1==2) || (recgen==-1 && nb2-nb2g==3 && semi2>0 && gm_ddst2==2) ) && ( dst1self!=0&&d1self==1&&l1self>0 && l1self==2 ) ) self =  2; // tag
    else if( ( abs(event_type)==rm_bb )                                                          ) self =  0; // fake
    else if( ( abs(event_type)!=rm_bb&&abs(event_type)>2 )                                       ) self = -1; // other BB(charged)
    else if( ( abs(event_type)!=rm_bb&&abs(event_type)<3 )                                       ) self = -2; // continuum
    dist->column( "self", self );

    if( fl_dump ) dist->dumpData();
    if( fl_message && fl_dump ) std::cout << " -> dump !" << std::endl;

    return;
  }

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
