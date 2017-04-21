#include "belle.h"
#include "DSTRTAUNU.h"
//#include "toolbox/FoxWolfr.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  double DSTRTAUNU::correct_dr( const Mdst_charged chg, HepPoint3D m_ip, int id )
  {
    double co_dr = 10000.;
    Mdst_trk& trk = chg.trk();
    Mdst_trk_fit& fit = trk.mhyp(id);

    if( trk && fit ){

      // Get Helix parameter
      HepVector a(5,0);
      a[0] = fit.helix(0);
      a[1] = fit.helix(1);
      a[2] = fit.helix(2);
      a[3] = fit.helix(3);
      a[4] = fit.helix(4);

      // Get the error matrix
      HepSymMatrix Ea(5,0);
      Ea[0][0] = fit.error(0);
      Ea[1][0] = fit.error(1); Ea[1][1] = fit.error(2);
      Ea[2][0] = fit.error(3); Ea[2][1] = fit.error(4);
      Ea[2][2] = fit.error(5);
      Ea[3][0] = fit.error(6); Ea[3][1] = fit.error(7);
      Ea[3][2] = fit.error(8);
      Ea[3][3] = fit.error(9);
      Ea[4][0] = fit.error(10); Ea[4][1] = fit.error(11);
      Ea[4][2] = fit.error(12);
      Ea[4][3] = fit.error(13); Ea[4][4] = fit.error(14);

      HepPoint3D pivot(fit.pivot(0), fit.pivot(1), fit.pivot(2));
      Helix crtrk(pivot, a, Ea);
      crtrk.pivot(m_ip);

      co_dr=crtrk.dr();
    }
    return(co_dr);
  }

  double DSTRTAUNU::correct_dz( const Mdst_charged chg, HepPoint3D m_ip, int id )
  {
    double co_dz = 10000.;
    Mdst_trk& trk = chg.trk();
    Mdst_trk_fit& fit = trk.mhyp(id);

    if( trk && fit ){

      // Get Helix parameter
      HepVector a(5,0);
      a[0] = fit.helix(0);
      a[1] = fit.helix(1);
      a[2] = fit.helix(2);
      a[3] = fit.helix(3);
      a[4] = fit.helix(4);

      // Get the error matrix
      HepSymMatrix Ea(5,0);
      Ea[0][0] = fit.error(0);
      Ea[1][0] = fit.error(1); Ea[1][1] = fit.error(2);
      Ea[2][0] = fit.error(3); Ea[2][1] = fit.error(4);
      Ea[2][2] = fit.error(5);
      Ea[3][0] = fit.error(6); Ea[3][1] = fit.error(7);
      Ea[3][2] = fit.error(8);
      Ea[3][3] = fit.error(9);
      Ea[4][0] = fit.error(10); Ea[4][1] = fit.error(11);
      Ea[4][2] = fit.error(12);
      Ea[4][3] = fit.error(13); Ea[4][4] = fit.error(14);

      HepPoint3D pivot(fit.pivot(0), fit.pivot(1), fit.pivot(2));
      Helix crtrk(pivot, a, Ea);
      crtrk.pivot(m_ip);

      co_dz=crtrk.dz();
    }
    return(co_dz);
  }

  int DSTRTAUNU::add_bremsstrahlung( Particle& particle, std::vector<Particle>& gamma_list, double angle_cut, bool fl_write, int fl_par1, int fl_par2, int fl_par3 )
  {
    int nchild   = particle.nChildren();
    int ngam_tot = 0;
    for( int i=0; i<nchild; i++ ){
      if( !(
	    abs(particle.child(i).pType().lund()) == Electron_LUND ||
	    abs(particle.child(i).pType().lund()) == fl_par1 ||
	    abs(particle.child(i).pType().lund()) == fl_par2 ||
	    abs(particle.child(i).pType().lund()) == fl_par3
	    )
	  ) continue;
      int ngam = 0; // the number of added gamma
      HepLorentzVector gamma_momentum;
      for(std::vector<Particle>::iterator gam = gamma_list.begin();
	  gam != gamma_list.end();
	  gam++)
	{
	  double angle = particle.child(i).p3().angle( gam->p3() ); // bug fixed @ 20140312
	  if( angle < angle_cut  ){
	    int fl_dupli = 0;
	    for( int j=0; j<nchild; j++ ){
	      Particle already_added_gamma = particle.child(j);
	      int gamID1 = already_added_gamma.mdstGamma().get_ID();
	      int gamID2 = gam->mdstGamma().get_ID();
	      if( gamID1 == gamID2 ){
		fl_dupli = 1;
		break;
	      }
	    }
	    if( fl_dupli )  continue;

	    gamma_momentum += gam->p();
	    particle.relation().append( *gam );
	    ngam++;
	  }
	}
      particle.momentum().momentum( particle.p() + gamma_momentum );
      if( fl_write ) dynamic_cast<UserInfo&>(particle.child(i).userInfo()).radgam_n( ngam );
      ngam_tot += ngam;
    }
    dynamic_cast<UserInfo&>(particle.userInfo()).radgam_n( ngam_tot );
    return ngam_tot;
  }

  
  // ---------------------------------------------------
  //   Check duplications of the final-state particles
  //   1 -> not duplication, 0 -> duplication
  // ---------------------------------------------------
  int DSTRTAUNU::check_dupli_daughter( std::vector<Particle>::iterator i, std::vector<Particle>::iterator j )
  {
    Particle& pi = *i;
    Particle& pj = *j;
    return check_dupli_daughter( pi, pj );
  }
  
  
  int DSTRTAUNU::check_dupli_daughter( Particle i, Particle j )
  {
      int b_nonduplication_flag = 1;
    const int chk1_nptcle = i.relation().nFinalStateParticles();
    const int chk2_nptcle = j.relation().nFinalStateParticles();
    for (int  k = 0; k < chk1_nptcle; k++) {
      for (int  l = 0; l < chk2_nptcle; l++) {
	Particle chk1_child = i.relation().finalStateParticle(k);
	Particle chk2_child = j.relation().finalStateParticle(l);
	if ( !chk1_child.mdstCharged() || !chk2_child.mdstCharged() ) continue;
	int chk1_chgID = chk1_child.mdstCharged().get_ID();
	int chk2_chgID = chk2_child.mdstCharged().get_ID();
	if( chk1_chgID * chk2_chgID != 0  && (chk1_chgID == chk2_chgID) ) b_nonduplication_flag = 0;
	if( b_nonduplication_flag == 0 ) break;
      }
      if( b_nonduplication_flag == 0 ) break;
    }

    for (int  k = 0; k < chk1_nptcle; k++) {
      for (int  l = 0; l < chk2_nptcle; l++) {
	Particle chk1_child = i.relation().finalStateParticle(k);
	Particle chk2_child = j.relation().finalStateParticle(l);
	
	if ( !chk1_child.mdstGamma() || !chk2_child.mdstGamma() ) continue;
	int chk1_gamID = chk1_child.mdstGamma().get_ID();
	int chk2_gamID = chk2_child.mdstGamma().get_ID();
	if( chk1_gamID * chk2_gamID != 0  && (chk1_gamID == chk2_gamID) ) b_nonduplication_flag = 0;
	if( b_nonduplication_flag == 0 ) break;
      }
      if( b_nonduplication_flag == 0 ) break;
    }
    
    for (int  m = 0; m < chk1_nptcle; m++) {
      for (int  n = m+1; n < chk1_nptcle; n++) {
	Particle d_child1 = i.relation().finalStateParticle(m);
	Particle d_child2 = i.relation().finalStateParticle(n); 
	if ( !d_child1.mdstCharged() || !d_child2.mdstCharged() ) continue;
	int d_ch1_chgID = d_child1.mdstCharged().get_ID();
	int d_ch2_chgID = d_child2.mdstCharged().get_ID();
	if(d_ch1_chgID * d_ch2_chgID != 0  && (d_ch1_chgID == d_ch2_chgID))
	  b_nonduplication_flag = 0;
	if( b_nonduplication_flag == 0 ) break;
      }
      if( b_nonduplication_flag == 0 ) break;
    }

    for (int  m = 0; m < chk1_nptcle; m++) {
      for (int  n = m+1; n < chk1_nptcle; n++) {
	Particle d_child1 = i.relation().finalStateParticle(m);
	Particle d_child2 = i.relation().finalStateParticle(n); 
	if ( !d_child1.mdstGamma() || !d_child2.mdstGamma() ) continue;
	int d_ch1_gamID = d_child1.mdstGamma().get_ID();
	int d_ch2_gamID = d_child2.mdstGamma().get_ID();
	if(d_ch1_gamID * d_ch2_gamID != 0  && (d_ch1_gamID == d_ch2_gamID))
	  b_nonduplication_flag = 0;
	if( b_nonduplication_flag == 0 ) break;
      }
      if( b_nonduplication_flag == 0 ) break;
    }
    
    for (int  m = 0; m < chk2_nptcle; m++) {
      for (int  n = m+1; n < chk2_nptcle; n++) {
	Particle d_child1 = j.relation().finalStateParticle(m);
	Particle d_child2 = j.relation().finalStateParticle(n); 
	if ( !d_child1.mdstCharged() || !d_child2.mdstCharged() ) continue;
	int d_ch1_chgID = d_child1.mdstCharged().get_ID();
	int d_ch2_chgID = d_child2.mdstCharged().get_ID();
	if(d_ch1_chgID * d_ch2_chgID != 0  && (d_ch1_chgID == d_ch2_chgID))
	  b_nonduplication_flag = 0;
	if( b_nonduplication_flag == 0 ) break;
      }
      if( b_nonduplication_flag == 0 ) break;
    }
    
    for (int  m = 0; m < chk2_nptcle; m++) {
      for (int  n = m+1; n < chk2_nptcle; n++) {
	Particle d_child1 = j.relation().finalStateParticle(m);
	Particle d_child2 = j.relation().finalStateParticle(n); 
	if ( !d_child1.mdstGamma() || !d_child2.mdstGamma() ) continue;
	int d_ch1_gamID = d_child1.mdstGamma().get_ID();
	int d_ch2_gamID = d_child2.mdstGamma().get_ID();
	if(d_ch1_gamID * d_ch2_gamID != 0  && (d_ch1_gamID == d_ch2_gamID))
	  b_nonduplication_flag = 0;
	if( b_nonduplication_flag == 0 ) break;
      }
      if( b_nonduplication_flag == 0 ) break;
    }
    
    return b_nonduplication_flag;
  }

  int DSTRTAUNU::kinematic_fit( Particle& particle, int flag_kfitter ) // for D
  {
    
    setUserInfo( particle );
    UserInfo& info = dynamic_cast<UserInfo&>( particle.userInfo() );
    info.kf( flag_kfitter );
    /* removed @20141107
    // pre-fit to know vertex
    int cnt_chg = 0;
    kvertexfitter pre_kvf;
    for( int i=0; i<particle.nChildren(); i++ ){
      if( abs(particle.child(i).pType().lund()) == Kplus_LUND    ||
	  abs(particle.child(i).pType().lund()) == PIplus_LUND   ||
	  abs(particle.child(i).pType().lund()) == Electron_LUND ||
	  abs(particle.child(i).pType().lund()) == MUminus_LUND  ||
	  abs(particle.child(i).pType().lund()) == Proton_LUND ){
	addTrack2fit( pre_kvf, particle.child(i) );
	cnt_chg++;
      }
    }
    addTube2fit(pre_kvf); // for IP tube constraint
    unsigned pre_err = pre_kvf.fit();
    */
    //unsigned pre_err = 1; // tmpppp
    //kvertexfitter pre_kvf; // tmpppp

    //std::cout << "mode : " << info.rec_mode() << ", "
    //<< "chg : "  << cnt_chg
    //<< std::endl;

    if( flag_kfitter==1 ){ // kvertexfitter *********************************************************
      kvertexfitter kf;
      for( int i=0; i<particle.nChildren(); i++ ) addTrack2fit( kf, particle.child(i) );
      /* removed @20141107
      if( pre_err==0 ){
	kf.initialVertex( pre_kvf.vertex() );
	kf.knownVertex();
      }
      */
      addTube2fit(kf); // for IP tube constraint // added @20141107
      unsigned err = kf.fit();
      if( err==0 ){ // success.
	info.kf_cl   ( kf.cl()      );
	info.kf_chisq( kf.chisq()   );
	info.kf_ndf  ( kf.dgf()     );
	makeMother   ( kf,particle  );
      }else{ // fail
	//std::cout << "false kvertexfitter : error code " << err << std::endl;
	return 0;
      }
    }else if( flag_kfitter==2 ){ // kmassvertexfitter ****************************************************************************
      kmassvertexfitter kf;
      for( int i=0; i<particle.nChildren(); i++ ) addTrack2fit( kf, particle.child(i) );
      kf.invariantMass( particle.pType().mass() );
      unsigned err = kf.fit();
      if( err==0 ){ // success.
	info.kf_cl   ( kf.cl()     );
	info.kf_chisq( kf.chisq()  );
	info.kf_ndf  ( kf.dgf()    );
	makeMother   ( kf,particle );
      }else{ // fail
	//std::cout << "false kmassvertexfitter : error code " << err << std::endl;
	return 0;
      }
    }else if( flag_kfitter==3 ){ // kmassfitter
      kmassfitter kf;
      for( int i=0; i<particle.nChildren(); i++ ) addTrack2fit( kf, particle.child(i) );
      kf.invariantMass(particle.pType().mass());
      unsigned err = kf.fit();
      if( err==0 ){ // success.
	info.kf_cl   ( kf.cl()     );
	info.kf_chisq( kf.chisq()  );
	info.kf_ndf  ( kf.dgf()    );
	makeMother   ( kf,particle );
      }else{ // fail
	//std::cout << "false kmassvertexfitter : error code " << err << std::endl;
	return 0;
      }
    }
    
    return 1;
  }


  int DSTRTAUNU::masscut(Particle particle, const double low, const double high, const double mass){
    // mass(PDG) + low < mass(REC) < mass(PDG) + high
    // 1 -> masscut through, 0 -> masscut rejection
    if( particle.momentum().mass() - mass < high && particle.momentum().mass() - mass > low ) return 1;
    else return 0;
  }
  
  void DSTRTAUNU::cal_Mmiss_Evis( BelleTuple* dist, std::vector<Particle>& trk_list, std::vector<Particle>& gamma_list,
				  const char* mmiss_name, const char* evis_name,
				  Particle* veto_pi0 ){
    
    // [ Combinatorial B.G. Suppression ]
    const double dr_cut1 = 1;
    const double dz_cut1 = 5;
    
    double Evis_trk=0;
    double Evis_gam=0;
    
    Hep3Vector Ptot_trk(0);
    Hep3Vector Ptot_gam(0);
    
    // <tracks>
    for( std::vector<Particle>::iterator p = trk_list.begin();
	 p != trk_list.end();
	 p++ )
      {
	UserInfo& info_trk  = dynamic_cast<UserInfo&>( p->userInfo()  );
	if( info_trk.dr() < dr_cut1 && info_trk.dz() < dz_cut1 ){
	  Evis_trk += info_trk.Vcm().e();
	  Ptot_trk += info_trk.Vcm().vect();
	}
      }	

    // <gamma>
    for( std::vector<Particle>::iterator g = gamma_list.begin();
	 g != gamma_list.end();
	 g++ )
      {
	UserInfo& info_gam = dynamic_cast<UserInfo&>( g->userInfo() );
	if( veto_pi0!=NULL && (
			       veto_pi0->child(0).mdstGamma().get_ID()==g->mdstGamma().get_ID() || 
			       veto_pi0->child(1).mdstGamma().get_ID()==g->mdstGamma().get_ID() )
	    ) continue;
	Evis_gam += info_gam.Vcm().e();
	Ptot_gam += info_gam.Vcm().vect();
      }
    
    double     mmiss2=0;
    double     Evis=0;
    Hep3Vector Ptot(0);
    
    Evis = Evis_trk + Evis_gam;
    Ptot = Ptot_trk + Ptot_gam;

    mmiss2 = (2*eb-Evis)*(2*eb-Evis) - Ptot.mag2();
    dist->column( evis_name,  Evis                                        );
    dist->column( mmiss_name, mmiss2 > 0 ? sqrt(mmiss2) : -sqrt(-mmiss2)  );
    
    return;
  }

  void DSTRTAUNU::cal_Mmiss_Evis( BelleTuple* dist, Particle par1, Particle par2,
				  const char* mmiss_name, const char* evis_name,
				  Particle* veto_pi0 ){
    double     mmiss2 = 0;
    double     Evis   = 0;
    Hep3Vector Ptot(0);

    UserInfo& info_par1  = dynamic_cast<UserInfo&>( par1.userInfo()  );
    Evis += info_par1.Vcm().e();
    Ptot += info_par1.Vcm().vect();

    UserInfo& info_par2  = dynamic_cast<UserInfo&>( par2.userInfo()  );
    Evis += info_par2.Vcm().e();
    Ptot += info_par2.Vcm().vect();

    if( veto_pi0!=NULL ){
      HepLorentzVector slow_pion_4Vcm = veto_pi0->p();
      slow_pion_4Vcm.boost( cmboost );
      Evis -= slow_pion_4Vcm.e();
      Ptot -= slow_pion_4Vcm.vect();
    }
    
    mmiss2 = (2*eb-Evis)*(2*eb-Evis) - Ptot.mag2();

    dist->column( evis_name,  Evis                                        );
    dist->column( mmiss_name, mmiss2 > 0 ? sqrt(mmiss2) : -sqrt(-mmiss2)  );

    return;
  }
  
  double DSTRTAUNU::cos2track( Particle par1, Particle par2 ){

    Hep3Vector par1_mom(par1.px(), par1.py(), par1.pz());
    Hep3Vector par2_mom(par2.px(), par2.py(), par2.pz());
    double cos = (par1_mom * par2_mom)/(par1_mom.mag()*par2_mom.mag());
    return cos;
}

  void DSTRTAUNU::qq_suppress(Particle b, double& thrust_angle, double& R2){
    // signal (rec-B) side
    std::vector<Hep3Vector> alltrk;
    std::vector<Hep3Vector> sigtrk;
    std::vector<Hep3Vector> othtrk;
    std::set<int> sig_id;
    std::set<int> gam_id;
    //std::cout << "**************** < thrust_angle start > *******************"<< std::endl;    
    for(int i=0; i<b.relation().nFinalStateParticles(); i++)
      {
	HepLorentzVector b_ch(b.relation().finalStateParticle(i).p());
	b_ch.boost( cmboost );
	sigtrk.push_back( b_ch.vect() );
	alltrk.push_back( b_ch.vect() );
	if( b.relation().finalStateParticle(i).mdstCharged().get_ID() ){
	  sig_id.insert( b.relation().finalStateParticle(i).mdstCharged().get_ID() );
	}else{
	  gam_id.insert( b.relation().finalStateParticle(i).mdstGamma().get_ID() );
	}
      }

    // remains
    //std::cout << " remains " << std::endl;
    Mdst_charged_Manager& chgMgr = Mdst_charged_Manager::get_manager();
    for(Mdst_charged_Manager::iterator c = chgMgr.begin();
	c != chgMgr.end();
	c++)
      {
	//std::cout << "mdst_charged : " << c->get_ID() << " -> ";
	if( sig_id.find( c->get_ID() ) != sig_id.end() ) {/*std::cout << "continue" << std::endl;*/ continue;}
	Particle tmp( (*c), Ptype(-211) );
	HepLorentzVector tmp_vec( tmp.p());
	tmp_vec.boost( cmboost );
	othtrk.push_back( tmp_vec.vect() );
	alltrk.push_back( tmp_vec.vect() );
	//std::cout << "push_back" << std::endl;
      }
    
    Mdst_gamma_Manager& gamMgr = Mdst_gamma_Manager::get_manager();
    for( Mdst_gamma_Manager::iterator gam = gamMgr.begin(); gam!=gamMgr.end(); gam++){
      //std::cout << "mdst_gamma : " << gam->get_ID() << " -> ";
      if( gam_id.find( gam->get_ID() ) != gam_id.end() ){/*std::cout << "continue(id)" << std::endl;*/ continue;}
      Particle tmp_b_gam( *gam );
      HepLorentzVector b_gam( tmp_b_gam.p() );
      b_gam.boost( cmboost );
      othtrk.push_back( b_gam.vect() );
      alltrk.push_back( b_gam.vect() );
      //std::cout << "push_back" << std::endl;
    }

    Vector3 thrust_signal = thrust( sigtrk.begin(), sigtrk.end(), SelfFunc(Vector3()) );
    Vector3 thrust_other  = thrust( othtrk.begin(), othtrk.end(), SelfFunc(Vector3()) );
    double  costhr        = thrust_signal.unit().dot( thrust_other.unit() );
    thrust_angle          = fabs(costhr);
    FoxWolfram fw( foxwolfram(alltrk.begin(), alltrk.end(), SelfFunc(Vector3())) );
    R2 = fw.R(2);
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
