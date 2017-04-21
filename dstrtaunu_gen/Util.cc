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

  int DSTRTAUNU::add_bremsstrahlung(Particle& particle, std::vector<Particle>& gamma_list)
  {
    int nchild   = particle.nChildren();
    int ngam_tot = 0;
    for( int i=0; i<nchild; i++ ){
      if( abs(particle.child(i).pType().lund()) != Electron_LUND ) continue;
      int ngam = 0; // the number of added gamma
      HepLorentzVector gamma_momentum;
      for(std::vector<Particle>::iterator gam = gamma_list.begin();
	  gam != gamma_list.end();
	  gam++)
	{
	  double angle = particle.child(i).p3().angle( gam->p3() );
	  //std::cout << "angle = " << angle << std::endl;
	  if( angle < Rad_angle_cut ){
	    gamma_momentum += gam->p();
	    particle.relation().append( *gam );
	    ngam++;
	  }
	}
      particle.momentum().momentum( particle.p() + gamma_momentum );
      dynamic_cast<UserInfo&>(particle.child(i).userInfo()).radgam_n( ngam );
      ngam_tot += ngam;
    }
    dynamic_cast<UserInfo&>(particle.userInfo()).radgam_n( ngam_tot );
    return ngam_tot;
  }

  
  // ---------------------------------------------------
  //   Check duplications of the final-state particles
  //   1 -> not duplication, 0 -> duplication
  // ---------------------------------------------------
  int DSTRTAUNU::check_dupli_daughter( std::vector<Particle>::iterator i,
				       std::vector<Particle>::iterator j )
  {
    int b_nonduplication_flag = 1;
    const int chk1_nptcle = i->relation().nFinalStateParticles();
    const int chk2_nptcle = j->relation().nFinalStateParticles();
    for (int  k = 0; k < chk1_nptcle; k++) {
      for (int  l = 0; l < chk2_nptcle; l++) {
	Particle chk1_child = i->relation().finalStateParticle(k);
	Particle chk2_child = j->relation().finalStateParticle(l);

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
	Particle chk1_child = i->relation().finalStateParticle(k);
	Particle chk2_child = j->relation().finalStateParticle(l);
	
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
	Particle d_child1 = i->relation().finalStateParticle(m);
	Particle d_child2 = i->relation().finalStateParticle(n); 
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
	Particle d_child1 = i->relation().finalStateParticle(m);
	Particle d_child2 = i->relation().finalStateParticle(n); 
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
	Particle d_child1 = j->relation().finalStateParticle(m);
	Particle d_child2 = j->relation().finalStateParticle(n); 
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
	Particle d_child1 = j->relation().finalStateParticle(m);
	Particle d_child2 = j->relation().finalStateParticle(n); 
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

  int DSTRTAUNU::kinematic_fit( Particle& particle, int flag_kfitter,
			  int ntrk)
  {
    if(!ntrk) ntrk = particle.nChildren();
    //std::cout << "kinematic fitting " << flag_kfitter
    //<< " : Particle( "      << particle.pType().name() << ",  "
    //<< ntrk                 << " tracks, "
    //<< particle.nChildren() << "children )"
    //<< std::endl;
    
    setUserInfo( particle );
    UserInfo& info = dynamic_cast<UserInfo&>( particle.userInfo() );
    info.kf( flag_kfitter );
    if( flag_kfitter==1 || flag_kfitter==2 ){ // kvertexfitter (+ IP tube) *********************************************************
      kvertexfitter kf;
      for( int i=0; i<particle.nChildren(); i++ ) addTrack2fit( kf, particle.child(i) );      
      if( flag_kfitter==2 ) addTube2fit(kf); // for IP tube constraint
      unsigned err = kf.fit();
      
      if( err==0 ){ // success.
	info.kf_cl   ( kf.cl()      );
	info.kf_chisq( kf.chisq()   );
	info.kf_ndf  ( kf.dgf()     );
	makeMother   ( kf, particle );
      }else{ // false.
	//std::cout << "false kvertexfitter : error code " << err << std::endl;
	return 0;
      }
    }else if( flag_kfitter==3 ){ // kmassvertexfitter ****************************************************************************
      if( particle.nChildren()-ntrk ){ // set gamma error(find vertex)
	kvertexfitter tmp_kf;
	for( int i=0; i<ntrk; i++ ) addTrack2fit( tmp_kf, particle.child(i) );
	unsigned tmp_err = tmp_kf.fit();
	if(tmp_err != 0) return 0;
	for(int i=ntrk; i<particle.nChildren(); i++){
	  setGammaError( particle.child(i), tmp_kf.vertex(), tmp_kf.errVertex() );
	}
      }
      
      kmassvertexfitter kf;
      for(int i=0; i<particle.nChildren(); i++) addTrack2fit( kf, particle.child(i) );
      kf.invariantMass(particle.pType().mass());
      unsigned err = kf.fit();
      if( err==0 ){ // success.
	info.kf_cl( kf.cl() );
	info.kf_chisq( kf.chisq() );
	info.kf_ndf( kf.dgf() );
	makeMother(kf, particle);
      }else{ // false.
	//std::cout << "false kmassvertexfitter : error code " << err << std::endl;
	return 0;
      }	
    
    }else if(flag_kfitter==4){ // kmassfitter ***********************************************************************************
      kmassfitter kf;
      for( int i=0; i<particle.nChildren(); i++ ){
	if( abs(particle.child(i).pType().lund()) == Kplus_LUND    ||
	    abs(particle.child(i).pType().lund()) == PIplus_LUND   ||
	    abs(particle.child(i).pType().lund()) == Electron_LUND ||
	    abs(particle.child(i).pType().lund()) == MUminus_LUND  ||
	    abs(particle.child(i).pType().lund()) == Proton_LUND )
	addTrack2fit( kf, particle.child(i) );
      }
      kf.invariantMass( particle.pType().mass() );
      unsigned err = kf.fit();
      if( err==0 ){ // success.
	info.kf_cl   ( kf.cl()      );
	info.kf_chisq( kf.chisq()   );
	info.kf_ndf  ( kf.dgf()     );
	makeMother   ( kf, particle );
      }else{ // false.
	//std::cout << "false kmassfitter : error code " << err << std::endl;
	return 0;
      }
    }else{
      info.kf_cl   ( -1 );
      info.kf_chisq( -1 );
      info.kf_ndf  ( -1 );
    }
    return 1;
  }

  int DSTRTAUNU::masscut(Particle particle, const double low, const double high, const double mass){
    // mass(PDG) + low < mass(REC) < mass(PDG) + high
    // 1 -> masscut through, 0 -> masscut rejection
    if( particle.momentum().mass() - mass < high && particle.momentum().mass() - mass > low ) return 1;
    else return 0;

  }
  
  void DSTRTAUNU::cal_Mmiss_Evis( std::vector<Particle>& trk_list, std::vector<Particle>& gamma_list ){
    
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
	Evis_gam += info_gam.Vcm().e();
	Ptot_gam += info_gam.Vcm().vect();
      }
    
    double     mmiss2=0;
    double     Evis=0;
    Hep3Vector Ptot(0);
    
    Evis = Evis_trk + Evis_gam;
    Ptot = Ptot_trk + Ptot_gam;
    mmiss2 = (2*eb-Evis)*(2*eb-Evis) - Ptot.mag2();
    Rec_dist->column( "evis",  Evis                                        );
    Rec_dist->column( "mmiss", mmiss2 > 0 ? sqrt(mmiss2) : -sqrt(-mmiss2)  );
    
    return;
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
