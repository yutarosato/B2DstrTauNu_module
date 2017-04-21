#include "DSTRTAUNU.h"

using namespace std;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



  void DSTRTAUNU::ckEclEnergyWithMatch2GoodGamma( Particle& B1, Particle& B2, const int fl_self, const int neecl, const bool fl_message,
						  const double width, const double e9e25,
						  const double th_fw, const double th_br, const double th_bw
						  ){
    Mdst_ecl_Manager&     eclmgr = Mdst_ecl_Manager::get_manager();
    Mdst_sim_ecl_Manager& simmgr = Mdst_sim_ecl_Manager::get_manager();
    
    Mdst_sim_ecl_Manager::iterator isim =simmgr.begin();
    for( Mdst_ecl_Manager::iterator iecl = eclmgr.begin();
	 iecl!=eclmgr.end(); iecl++ ){
      if( iecl!=eclmgr.begin() ) isim++;
      // check matching with charged tracks
      if( iecl->match()==1 || (iecl->match()==2&&!good_gamma(*iecl, 0.02, e9e25, width)) ) continue;
      
      // check if cluster is used for the reconstruction
      if( myCheckSame( B1, *iecl ) ) continue;
      if( myCheckSame( B2, *iecl ) ) continue;
      
      // check if cluster energy exceed the threshold
      int fl_region    = checkEECLCluster( iecl->energy(), iecl->theta(),   0.0,   0.0,   0.0 ); // 1(forward), 2(barrel), 3(backward) without threshold
      int fl_threshold = checkEECLCluster( iecl->energy(), iecl->theta(), th_fw, th_br, th_bw ); // 1(forward), 2(barrel), 3(backward)   with  threshold
      int fl_origin    = ( mc==0 ? 0 : checkEECLClusterOrigin( isim->hepevt(), fl_message ) );
      
      if( mc!=0 && isim->ecl().energy() != iecl->energy() ) std::cerr << "Wrong iterator between Mdst_sim_ecl and Mdst_ecl" << std::endl, abort();


      double exprunNo = expNo*10000+runNo;
      if( runNo==0 ) exprunNo = expNo*10000+runNoforsigMC; // for signal MC
      eecl_dist->column( "exprun",    exprunNo       );
      eecl_dist->column( "event",     evtNo          );
      eecl_dist->column( "evt_type",  event_type     );

      eecl_dist->column( "eecl",      iecl->energy() );
      eecl_dist->column( "region",    fl_region      );
      eecl_dist->column( "origin",    fl_origin      );
      eecl_dist->column( "threshold", fl_threshold   );
      eecl_dist->column( "self",      fl_self        );
      eecl_dist->column( "neecl",     neecl          );
      eecl_dist->column( "dtau1",     B1.mass()      );
      eecl_dist->column( "dtau2",     B2.mass()      );
      eecl_dist->dumpData();
    }

  }
    

  double DSTRTAUNU::calEclEnergyWithMatch2GoodGamma( const Particle& B1, const Particle& B2, const bool fl_message,
						     double eecl[39], int neecl[39],
						     const Particle* veto_pi0,
						     const int id1, const int id2,
						     const double width, const double e9e25,
						     const double th_fw, const double th_br, const double th_bw
						     )
  {
    Mdst_ecl_Manager&     eclmgr = Mdst_ecl_Manager::get_manager();
    Mdst_sim_ecl_Manager& simmgr = Mdst_sim_ecl_Manager::get_manager();
    double eclenergy(0);
    //if( fl_message ) std::cout << std::endl
    //<< "[ECL calculation : region(1[forward], 2[barrel], 3[backward]), side(1[B1], 2[B2], 3[other])"
    //<< std::endl;

    Mdst_sim_ecl_Manager::iterator isim =simmgr.begin();
    for( Mdst_ecl_Manager::iterator iecl = eclmgr.begin();
	 iecl!=eclmgr.end(); iecl++ ){
      if( iecl!=eclmgr.begin() ) isim++;
      // check matching with charged tracks
      if( iecl->match()==1 || (iecl->match()==2&&!good_gamma(*iecl, 0.02, e9e25, width)) ) continue;
      
      // check if cluster energy exceed the threshold
      int fl_region = checkEECLCluster( iecl->energy(), iecl->theta(), th_fw, th_br, th_bw ); // 1(forward), 2(barrel), 3(backward)
      if( !fl_region ) continue;

       // check if cluster is used for the reconstruction
      if( veto_pi0==NULL || (veto_pi0!=NULL && !myCheckSame(*veto_pi0, *iecl)) ){ // control samples study for D**
	if( myCheckSame( B1, *iecl ) ) continue;
	if( myCheckSame( B2, *iecl ) ) continue;
      }else{ // modified for extra pi0 samples @20150728
	continue;
      }
      
      if( mc!=0 && isim->ecl().energy() != iecl->energy() ) std::cerr << "Wrong iterator between Mdst_sim_ecl and Mdst_ecl" << std::endl, abort();

      int fl_side = ( mc==0 ? -1 : which_B( gen_level(isim->hepevt()), id1, id2 ) ); // 1(B1), 2(B2), -1(other)
      if( fl_side==-1 ) fl_side = 3;                                                 //             -> 3(other)
      /*
      if( fl_message ){
	std::cout << "[ECLCluster] "
		  << "match = "  << iecl->match()  << ", "
		  << "energy = " << iecl->energy() << ", "
		  << "r = "      << iecl->r()      << ", "
		  << "phi = "    << iecl->phi()    << ", "
		  << "theta = "  << iecl->theta()  << ", "
		  << "region = " << fl_region      << ", "
		  << "side = "   << fl_side        << ", ";
	//int nchild = roothep.daLast()==0 ? 0 : roothep.daLast() - roothep.daFirst() + 1;
	//std::cout << "nchild = " << nchild << std::endl;
      }
      */
      
      int fl_origin = ( mc==0 ? 0 : checkEECLClusterOrigin( isim->hepevt(), fl_message ) );
      
      eclenergy += iecl->energy();
      
      if( eecl!=NULL ){
	eecl[fl_side  -1] += iecl->energy(); // [0-2] 0(B1), 1(B2), 2(other->beam B.G.)
	eecl[fl_region+2] += iecl->energy(); // [3-5] 3(forward), 4(barrel), 6(backward)
	eecl[fl_origin+6] += iecl->energy(); // [6-34]
	if( fl_origin==10 ) eecl[36+fl_region] += iecl->energy(); // [35-38] for D** study @ 20141002
      }
      if( neecl!=NULL ){
	neecl[fl_side  -1]++; // [0-2] 0(B1), 1(B2), 2(other->beam B.G.)
	neecl[fl_region+2]++; // [3-5] 3(forward), 4(barrel), 6(backward)
	neecl[fl_origin+6]++; // [6-34]
	if( fl_origin==10 ) neecl[36+fl_region]++; // [35-38] for D** study @ 20141002 
      }
    }
    return eclenergy;
  }


  double DSTRTAUNU::calEclEnergyWithMatch2GoodGamma_mask( const Particle& B1, const Particle& B2, const bool fl_message,
							  double eecl[39], int neecl[39],
							  const Particle* l_lep,
							  const Particle* veto_pi0,
							  const int id1, const int id2,
							  const double width, const double e9e25,
							  const double th_fw, const double th_br, const double th_bw
							  )
  {
    Mdst_ecl_Manager&     eclmgr = Mdst_ecl_Manager::get_manager();
    Mdst_sim_ecl_Manager& simmgr = Mdst_sim_ecl_Manager::get_manager();
    double eclenergy(0);
    //if( fl_message ) std::cout << std::endl
    //<< "[ECL calculation : region(1[forward], 2[barrel], 3[backward]), side(1[B1], 2[B2], 3[other])"
    //<< std::endl;

    Mdst_sim_ecl_Manager::iterator isim =simmgr.begin();
    for( Mdst_ecl_Manager::iterator iecl = eclmgr.begin();
	 iecl!=eclmgr.end(); iecl++ ){
      if( iecl!=eclmgr.begin() ) isim++;
      // check matching with charged tracks
      //if( ( mc==0 ? 0 : checkEECLClusterOrigin( isim->hepevt(), fl_message ) )==16 ){
      Hep3Vector p_lep( (HepDouble)l_lep->px(), (HepDouble)l_lep->py(), (HepDouble)l_lep->pz() );
      Hep3Vector p_ecl(1,1,1);
      p_ecl.setTheta(iecl->theta());
      p_ecl.setPhi  (iecl->phi()  );
      p_ecl.setMag  (iecl->r()    );
      /*
	std::cout << std::endl;
	std::cout << "origin = " << checkEECLClusterOrigin( isim->hepevt(), fl_message ) << std::endl;
	std::cout << "phi1 = "   << p_lep.phi()   << ", "
	<< "theta1 = " << p_lep.theta() << std::endl;
	std::cout << "phi2 = "   << p_ecl.phi()   << ", "
	<< "theta2 = " << p_ecl.theta() << std::endl;
	std::cout << "angle = "  << p_lep.angle( p_ecl ) << std::endl;
	std::cout << "match = " << iecl->match() << std::endl;
	std::cout << "energy = " << iecl->energy() << std::endl;
      */
      double l_lep_angle = p_lep.angle( p_ecl );
      Mdst_ecl_aux_Manager& aux = Mdst_ecl_aux_Manager::get_manager();
      Mdst_ecl_aux & shower_aux(aux(Panther_ID(iecl->get_ID())));
      /*
      std::cout << "e9/e25 = " << shower_aux.e9oe25() << std::endl;
      */
      //}
      if( iecl->match()==1 ) continue;
      ///*
      if( iecl->match()==2 ){
	if( l_lep_angle > 0.150 ){ // for hz1 and hz2
	  //if( l_lep_angle > 0.400 ){ // hz3
	  if( !good_gamma(*iecl, 0.02, e9e25, width) ) continue;
	}else{ // special region
	  // hz1 and hz3 (no cut)
	  //if( !good_gamma(*iecl, 0.02, 0.0,   width) ) continue; // hz2
	}
      }
      //*/
      
      
      // check if cluster energy exceed the threshold
      int fl_region = checkEECLCluster( iecl->energy(), iecl->theta(), th_fw, th_br, th_bw ); // 1(forward), 2(barrel), 3(backward)
      if( !fl_region ) continue;

       // check if cluster is used for the reconstruction
      if( veto_pi0==NULL || (veto_pi0!=NULL && !myCheckSame(*veto_pi0, *iecl)) ){ // control samples study for D**
	if( myCheckSame( B1, *iecl ) ) continue;
	if( myCheckSame( B2, *iecl ) ) continue;
      }else{ // modified for extra pi0 samples @20150728
	continue;
      }
      
      if( mc!=0 && isim->ecl().energy() != iecl->energy() ) std::cerr << "Wrong iterator between Mdst_sim_ecl and Mdst_ecl" << std::endl, abort();

      int fl_side = ( mc==0 ? -1 : which_B( gen_level(isim->hepevt()), id1, id2 ) ); // 1(B1), 2(B2), -1(other)
      if( fl_side==-1 ) fl_side = 3;              //                3(other)
      /*
      if( fl_message ){
	std::cout << "[ECLCluster] "
		  << "match = "  << iecl->match()  << ", "
		  << "energy = " << iecl->energy() << ", "
		  << "r = "      << iecl->r()      << ", "
		  << "phi = "    << iecl->phi()    << ", "
		  << "theta = "  << iecl->theta()  << ", "
		  << "region = " << fl_region      << ", "
		  << "side = "   << fl_side        << ", ";
	//int nchild = roothep.daLast()==0 ? 0 : roothep.daLast() - roothep.daFirst() + 1;
	//std::cout << "nchild = " << nchild << std::endl;
      }
      */
      
      int fl_origin = ( mc==0 ? 0 : checkEECLClusterOrigin( isim->hepevt(), fl_message ) );

      eclenergy += iecl->energy();
      
      if( eecl!=NULL ){
	eecl[fl_side  -1] += iecl->energy(); // [0-2] 0(B1), 1(B2), 2(other->beam B.G.)
	eecl[fl_region+2] += iecl->energy(); // [3-5] 3(forward), 4(barrel), 6(backward)
	eecl[fl_origin+6] += iecl->energy(); // [6-34]
	if( fl_origin==10 ) eecl[36+fl_region] += iecl->energy(); // [35-38] for D** study @ 20141002
      }
      if( neecl!=NULL ){
	neecl[fl_side  -1]++; // [0-2] 0(B1), 1(B2), 2(other->beam B.G.)
	neecl[fl_region+2]++; // [3-5] 3(forward), 4(barrel), 6(backward)
	neecl[fl_origin+6]++; // [6-34]
	if( fl_origin==10 ) neecl[36+fl_region]++; // [35-38] for D** study @ 20141002 
      }
    }
    return eclenergy;
  }


int DSTRTAUNU::checkEECLCluster( const double i_energy, const double i_theta, const double th_fw, const double th_br, const double th_bw )
{
  const double theta(i_theta*180./M_PI);
  if( theta<31.4 ){ // forward-endcap region
    if( i_energy>th_fw ) return 1;
    else                 return 0;
  }
  else if( theta>130.7 ){ // backward-endcap region
    if( i_energy>th_bw ) return 3;
    else                 return 0;
  }
  else{ // barrel region
    if( i_energy>th_br ) return 2;
    else                 return 0;
  }
  return 0;
}

  int DSTRTAUNU::checkEECLClusterOrigin( const Gen_hepevt& hep, bool fl_message ){
    // 0   (other?)
    // 1   (beam B.G.)
    // 2-7 (neutral particle : KS, Lambda, KL, neutron,baryon)
    // 8-13(gamma from pi0 decay)
    // 14  (gamma from eta decay)
    // 15  (gamma from D* decay)
    // 16  (brems from semileptonic B decay)
    // 17  (gamma from other decay)
    // 18-29(charged particles : K, pi, e, mu, proton, baryon)

    if( mc==0 ) return 0;
    
    const Gen_hepevt& roothep = gen_level( hep );
    if( fl_message ){
      std::cout << "     [ECLClusterOrigin] "
		<< hep.idhep()     << " (" << hep.isthep()     << ", " << hep.get_ID()     << ")" <<      " <- "
		<< roothep.idhep() << " (" << roothep.isthep() << ", " << roothep.get_ID() << ")";
      if( roothep.mother()                                                                                                                   ) std::cout << " <- " << roothep.mother().idhep()                            << " (" << roothep.mother().isthep()                            << ")";
      if( roothep.mother() && roothep.mother().mother()                                                                                      ) std::cout << " <- " << roothep.mother().mother().idhep()                   << " (" << roothep.mother().mother().isthep()                   << ")";
      if( roothep.mother() && roothep.mother().mother() && roothep.mother().mother().mother()                                                ) std::cout << " <- " << roothep.mother().mother().mother().idhep()          << " (" << roothep.mother().mother().mother().isthep()          << ")";
      if( roothep.mother() && roothep.mother().mother() && roothep.mother().mother().mother() && roothep.mother().mother().mother().mother() ) std::cout << " <- " << roothep.mother().mother().mother().mother().idhep() << " (" << roothep.mother().mother().mother().mother().isthep() << ")";
    }

    int fl_origin = 0;

    if     ( roothep.idhep()==911 ) fl_origin = 1; // 1(Beam B.G.)
    else if( checkEECLClusterOrigin_KsLambda( roothep ) ){ // [Lambda or Ks]
      int fl_kslambda = checkEECLClusterOrigin_KsLambda( roothep );
      if     ( fl_kslambda==1 ) fl_origin = 2; // 2(Lambda)
      else if( fl_kslambda==2 ) fl_origin = 3; // 3(Ks->pi+pi-)
      else if( fl_kslambda==3 ) fl_origin = 4; // 4(Ks->pi0pi0)
      else{
	//std::cout << "Strange Ks decay" << std::endl;
	//fl_origin = 0; // 0(for debug)
	fl_origin = 4; // rare case(Ks hadronization) is included as 4(Ks->pi0pi0)
      }
    }else if( roothep.idhep() == 130 ) fl_origin = 5; // 5(KL)
    else if( abs(roothep.idhep())==2112 ) fl_origin = 6; // 6(neutron)
    else if( abs(roothep.idhep())==3322 ) fl_origin = 7; // 7(neutral baryon) : xsi0
    else if( roothep.idhep()==22 && roothep.mother() && roothep.isthep()>0 ){ // [gamma from real-dacay]
      if ( roothep.mother().idhep()==111 && roothep.mother().mother() ){ // 8-12(gamma from pi0 decay)
	if     ( abs(roothep.mother().mother().idhep())==   411 || abs(roothep.mother().mother().idhep())==   421 ) fl_origin = 8; // 8(gamma from pi0 decay from D )
	else if( abs(roothep.mother().mother().idhep())==   413 || abs(roothep.mother().mother().idhep())==   423 ) fl_origin = 9; // 9(gamma from pi0 decay from D*)
	else if( abs(roothep.mother().mother().idhep())== 10411 || abs(roothep.mother().mother().idhep())== 10421 || // D**(L=1, broad)
		 abs(roothep.mother().mother().idhep())== 20413 || abs(roothep.mother().mother().idhep())== 20423 || // D**(L=1, broad)
		 abs(roothep.mother().mother().idhep())== 10413 || abs(roothep.mother().mother().idhep())== 10423 || // D**(L=1, narrow)
		 abs(roothep.mother().mother().idhep())==   415 || abs(roothep.mother().mother().idhep())==   425 || // D**(L=1, narrow)
		 abs(roothep.mother().mother().idhep())==100411 || abs(roothep.mother().mother().idhep())==100421 || // D**(radial excitation)
		 abs(roothep.mother().mother().idhep())==100413 || abs(roothep.mother().mother().idhep())==100423    // D**(radial excitation)
		 ) fl_origin = 10; // 10(gamma from pi0 decay from D**);
	else if( abs(roothep.mother().mother().idhep())==   433 || abs(roothep.mother().mother().idhep())== 20433 || abs(roothep.mother().mother().idhep())==431 ) fl_origin = 11; // 11(gamma from pi0 from Ds(*))
	else if( abs(roothep.mother().mother().idhep())==   511 || abs(roothep.mother().mother().idhep())== 521                                                  ) fl_origin = 12; // 12(gamma from pi0 from B)
	else fl_origin = 13; // 13(gamma from pi0 decay from others)
      }else if( roothep.mother().idhep()==221 ) fl_origin = 14; // 14(gamma from eta decay)
      else  if( abs(roothep.mother().idhep())==413 || abs(roothep.mother().idhep())==  423 ||
		abs(roothep.mother().idhep())==433 || abs(roothep.mother().idhep())==20433
		) fl_origin = 15;  // 15(gamma from D* meson decay) : D*0/D*+/Ds*/Ds1+
      else if( abs(roothep.mother().idhep())==511 || abs(roothep.mother().idhep())==521 ) fl_origin = 16; // 16(gamma from B ; Bremsstrahlung gammas of semileptonic B decay)
      else fl_origin = 17; // 17(gamma from others) : D0, rho(770)+, w(782), sigma0, xi_c'+, xi_c'0,...
    }
    else if( abs(roothep.idhep())== 321 && hep.isthep()>0 ) fl_origin = 18; // 18(K track)
    else if( abs(roothep.idhep())== 321 && hep.isthep()<0 ) fl_origin = 19; // 19(K shower)
    else if( abs(roothep.idhep())== 211 && hep.isthep()>0 ) fl_origin = 20; // 20(pi track)
    else if( abs(roothep.idhep())== 211 && hep.isthep()<0 ) fl_origin = 21; // 21(pi shower)
    else if( abs(roothep.idhep())==  11 && hep.isthep()>0 ) fl_origin = 22; // 22(electron track)
    else if( abs(roothep.idhep())==  11 && hep.isthep()<0 ) fl_origin = 23; // 23(electron shower)
    else if( abs(roothep.idhep())==  13 && hep.isthep()>0 ) fl_origin = 24; // 24(muon track)
    else if( abs(roothep.idhep())==  13 && hep.isthep()<0 ) fl_origin = 25; // 25(muon shower)
    else if( abs(roothep.idhep())==2212 && hep.isthep()>0 ) fl_origin = 26; // 26(proton track)
    else if( abs(roothep.idhep())==2212 && hep.isthep()<0 ) fl_origin = 27; // 27(proton shower)
    else if( abs(roothep.idhep())==3112 || abs(roothep.idhep())==3222 || abs(roothep.idhep())==3312 ){
      if( hep.isthep()>0 ) fl_origin = 28; // 28(charged baryon track ) : sigma, xsi
      else                 fl_origin = 29; // 29(charged baryon shower)
    }

    if( fl_message ) std::cout << " : -> " << fl_origin << ";" << std::endl;

    return fl_origin;
  }

  int DSTRTAUNU::checkEECLClusterOrigin_KsLambda( const Gen_hepevt& hep ){ // 1(Lambda), 2(Ks->pi+pi-), 3(Ks->pi0pi0), 4(other Ks decay?), 5(direct Ks?)
    if( mc==0 ) return 0;

    int fl_kslambda = 0;
    const Gen_hepevt* gen = &hep;
    while(1){
      if( gen->idhep()==3122 ){ // Lambda
	fl_kslambda = 1;
	break;
      }else if( gen->idhep()==310 ){ // Ks
	int nchild = gen->daLast()==0 ? 0 : gen->daLast() - gen->daFirst() + 1;
	if( nchild ){
	  Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
	  int ks_daughter = genMgr(Panther_ID(gen->daFirst())).idhep();
	  if     ( abs(ks_daughter)==211 ) fl_kslambda = 2; // Ks->pi+pi-
	  else if(     ks_daughter ==111 ) fl_kslambda = 3; // Ks->pi0pi0
	  else                             fl_kslambda = 4; // other Ks decay? 
	}else{ // direct Ks?
	  fl_kslambda = 5;
	}
	break;
      }else if( gen->mother() ){ // iteration
	gen = &gen->mother();
      }else break;
    }
    return fl_kslambda;
  }

  int DSTRTAUNU::cnt_remain_trk( const Particle& B1, const Particle& B2,
				 double dr_cut, double dz_cut
				 ){
    int cnt = 0;
    Mdst_charged_Manager& chgMgr = Mdst_charged_Manager::get_manager();
    for( Mdst_charged_Manager::iterator c = chgMgr.begin();
	 c != chgMgr.end();
	 c++ )
      {
	double dr = correct_dr( (*c), m_IP, PION_CODE );
	double dz = correct_dz( (*c), m_IP, PION_CODE );
	if( fabs(dr) > dr_cut ) continue;
	if( fabs(dz) > dz_cut ) continue;
	if( !myCheckSame( B1, *c ) && !myCheckSame( B2, *c ) ) cnt++;
      }
    return cnt;
  }

  int DSTRTAUNU::cnt_remain_pi0( const Particle& B1, const Particle& B2, std::vector<Particle>& pi0_list ){
    int cnt = 0;
    for( std::vector<Particle>::iterator p = pi0_list.begin(); p != pi0_list.end(); p++ ){
      if( !myCheckSame( B1, p->mdstPi0() ) && !myCheckSame( B2, p->mdstPi0() ) ) cnt++;
    }
    return cnt;
  }

  Particle* DSTRTAUNU::search_remain_pi0( const Particle& B1, const Particle& B2, std::vector<Particle>& pi0_list ){
    for( std::vector<Particle>::iterator p = pi0_list.begin(); p != pi0_list.end(); p++ ){
      if( !myCheckSame( B1, p->mdstPi0() ) && !myCheckSame( B2, p->mdstPi0() ) ){
	Particle& extra_pi0 = *p;
	Particle* par = &extra_pi0;
	return par;
      }
    }
    return NULL;
  }

  int DSTRTAUNU::cnt_remain_ks( const Particle& B1, const Particle& B2, std::vector<Particle>& ks_list ){
    int cnt = 0;
    for( std::vector<Particle>::iterator k = ks_list.begin(); k != ks_list.end(); k++ ){
      if( !myCheckSame( B1, k->child(0).mdstCharged() ) && !myCheckSame( B2, k->child(0).mdstCharged() ) &&
	  !myCheckSame( B1, k->child(1).mdstCharged() ) && !myCheckSame( B2, k->child(1).mdstCharged() )
	  ) cnt++;
    }
    return cnt;
  }

  int DSTRTAUNU::cnt_remain_gamma( const Particle& B1, const Particle& B2, std::vector<Particle>& gam_list ){
    int cnt = 0;
    for( std::vector<Particle>::iterator g = gam_list.begin(); g != gam_list.end(); g++ ){
      if( !myCheckSame( B1, g->mdstGamma() ) && !myCheckSame( B2, g->mdstGamma() ) ) cnt++;
    }
    return cnt;
  }

  int DSTRTAUNU::cnt_remain_gamma_energy( const Particle& B1, const Particle& B2, std::vector<Particle>& gam_list,
					  const double th_fw, const double th_br, const double th_bw ){
    int cnt = 0;
    for( std::vector<Particle>::iterator g = gam_list.begin(); g != gam_list.end(); g++ ){
      const Mdst_ecl& ecl  = g->mdstEcl();
      int fl_threshold = checkEECLCluster( ecl.energy(), ecl.theta(), th_fw, th_br, th_bw ); // 1(forward), 2(barrel), 3(backward)   with  threshold
      if( fl_threshold && !myCheckSame( B1, g->mdstGamma() ) && !myCheckSame( B2, g->mdstGamma() ) ) cnt++;
    }
    return cnt;
  }

  double DSTRTAUNU::cal_nearestgammadirection( const Particle par1, const Particle par2, const Particle& B1, const Particle& B2, std::vector<Particle>& gam_list ){
    double nearest_angle = 10;
    for( std::vector<Particle>::iterator g = gam_list.begin(); g != gam_list.end(); g++ ){
      if( !myCheckSame( B1, g->mdstGamma() ) && !myCheckSame( B2, g->mdstGamma() ) ){
	double angle1 = par1.p3().angle( g->p3() );
	double angle2 = par2.p3().angle( g->p3() );
	if( angle1 < nearest_angle ) nearest_angle = angle1;
	if( angle2 < nearest_angle ) nearest_angle = angle2;
      }
    }
    return nearest_angle;
  }
  

  double DSTRTAUNU::calBdirection( Particle& Dtau1, Particle& Dtau2 ){

    UserInfo& info_Dtau1 = dynamic_cast<UserInfo&>( Dtau1.userInfo() );
    UserInfo& info_Dtau2 = dynamic_cast<UserInfo&>( Dtau2.userInfo() );

    //calculate B vector
    Hep3Vector x( Dtau1.p().vect().unit() ); // p_A
    Hep3Vector y( Dtau2.p().vect().unit() ); // p_B
    Hep3Vector z( x.cross( y ).unit()     ); // (p_A x p_B) / |p_A x p_B|
    double xy = x.dot( y );
    double cosBDl1 = info_Dtau1.cos();
    double cosBDl2 = info_Dtau2.cos();

    double u  = ( cosBDl1 + xy * cosBDl2 );
    u /= ( 1 - xy * xy );
    double v  = -1 * ( cosBDl2 + xy * cosBDl1 );
    v /= ( 1 - xy * xy );
    double ww = 1 - u * u - v * v - 2 * u * v * xy;

    return ww; // if ww is negative, the calculation of b-direction fails;
    /*
    if( ww<0. ) return 0;
    double w  = sqrt( ww );
    
    
    B_plus = u * x + v * y + w * z;
    double bmass = ( Dtau1.charge() ? PDG_BplusMass : PDG_B0Mass );
    B_plus *= sqrt(eb*eb-bmass*bmass);
    
    B_minus = u * x + v * y - w * z;
    B_minus *= sqrt(eb*eb-bmass*bmass); 

    return 1;
    */
  }

  double DSTRTAUNU::calBdirection_with_more_bremsstrahlung( Particle Dtau1, Particle Dtau2, std::vector<Particle>& gamma_list,
							    double angle_cut,
							    int fl_par1, int fl_par2, int fl_par3 ){
    add_bremsstrahlung( Dtau1, gamma_list, angle_cut, false, fl_par1, fl_par2, fl_par3 );
    add_bremsstrahlung( Dtau2, gamma_list, angle_cut, false, fl_par1, fl_par2, fl_par3 );

    double Bmass = ( Dtau1.charge() ? PDG_B0Mass : PDG_BplusMass );

    UserInfo& info_Dtau1 = dynamic_cast<UserInfo&>( Dtau1.userInfo() );
    HepLorentzVector Dtau1_4Vcm = Dtau1.p();
    Dtau1_4Vcm.boost( cmboost );
    info_Dtau1.Vcm( Dtau1_4Vcm );
    double cosBDl1 = ( 2.0*eb*Dtau1_4Vcm.e() - Bmass*Bmass-Dtau1.mass()*Dtau1.mass() ) / ( 2.0*sqrt(eb*eb-Bmass*Bmass)*Dtau1_4Vcm.vect().mag() );
    info_Dtau1.cos( cosBDl1 );

    UserInfo& info_Dtau2 = dynamic_cast<UserInfo&>( Dtau2.userInfo() );
    HepLorentzVector Dtau2_4Vcm = Dtau2.p();
    Dtau2_4Vcm.boost( cmboost );
    info_Dtau2.Vcm( Dtau2_4Vcm );
    double cosBDl2 = ( 2.0*eb*Dtau2_4Vcm.e() - Bmass*Bmass-Dtau2.mass()*Dtau2.mass() ) / ( 2.0*sqrt(eb*eb-Bmass*Bmass)*Dtau2_4Vcm.vect().mag() );
    info_Dtau2.cos( cosBDl2 );
    
    return calBdirection( Dtau1, Dtau2 );
  }

  double DSTRTAUNU::cosBDl_with_more_bremsstrahlung( Particle Dtau, std::vector<Particle>& gamma_list,
						     double angle_cut,
						     int fl_par1, int fl_par2, int fl_par3 ){
    
    add_bremsstrahlung( Dtau, gamma_list, angle_cut, false, fl_par1, fl_par2, fl_par3 );
    double Bmass = ( Dtau.charge() ? PDG_B0Mass : PDG_BplusMass );
    UserInfo& info_Dtau = dynamic_cast<UserInfo&>( Dtau.userInfo() );
    HepLorentzVector Dtau_4Vcm = Dtau.p();
    Dtau_4Vcm.boost( cmboost );
    info_Dtau.Vcm( Dtau_4Vcm );
    double cosBDl = ( 2.0*eb*Dtau_4Vcm.e() - Bmass*Bmass-Dtau.mass()*Dtau.mass() ) / ( 2.0*sqrt(eb*eb-Bmass*Bmass)*Dtau_4Vcm.vect().mag() );
    info_Dtau.cos( cosBDl );
    return cosBDl;
  }

  void DSTRTAUNU::extra_pi0_study( Particle Dtau_l, Particle Dtau_h, Particle D_l, Particle D_h,
				   Particle* extra_pi0, char* name, BelleTuple* dist
				   ){
    char epi_pie[16]; sprintf( epi_pie,"%s_pie", name );
    char epi_gle[16]; sprintf( epi_gle,"%s_gle", name );
    char epi_ghe[16]; sprintf( epi_ghe,"%s_ghe", name );
    char epi_cl [16]; sprintf( epi_cl, "%s_cl",  name );
    char epi_ml [16]; sprintf( epi_ml, "%s_ml",  name );
    char epi_ch [16]; sprintf( epi_ch, "%s_ch",  name );
    char epi_mh [16]; sprintf( epi_mh, "%s_mh",  name );

    
    if( extra_pi0!=NULL ){
      dist->column( epi_pie,   extra_pi0->e() );
      dist->column( epi_gle, ((extra_pi0->child(0).e() > extra_pi0->child(1).e()) ? extra_pi0->child(1).e() : extra_pi0->child(0).e()) );
      dist->column( epi_ghe, ((extra_pi0->child(0).e() < extra_pi0->child(1).e()) ? extra_pi0->child(1).e() : extra_pi0->child(0).e()) );
      
      Particle Dtaupi0_l( Dtau_l.p() + extra_pi0->p(), Dtau_l.pType() );
      Dtaupi0_l.relation().append( Dtau_l     );
      Dtaupi0_l.relation().append( *extra_pi0 );
      double Bmass_l = ( Dtaupi0_l.charge() ? PDG_BplusMass : PDG_B0Mass );
      HepLorentzVector Dtaupi0_l_4Vcm = Dtaupi0_l.p();
      Dtaupi0_l_4Vcm.boost( cmboost );
      double cosBDl_extra = ( 2.0*eb*Dtaupi0_l_4Vcm.e() - Bmass_l*Bmass_l-Dtaupi0_l.mass()*Dtaupi0_l.mass() ) / ( 2.0*sqrt(eb*eb-Bmass_l*Bmass_l)*Dtaupi0_l_4Vcm.vect().mag() );
      dist->column( epi_cl, cosBDl_extra     );
      Particle Dpi0_l( D_l.p() + extra_pi0->p(), D_l.pType() );	  
      dist->column( epi_ml, Dpi0_l.mass() );
      
      Particle Dtaupi0_h( Dtau_h.p() + extra_pi0->p(), Dtau_h.pType() );
      Dtaupi0_h.relation().append( Dtau_h     );
      Dtaupi0_h.relation().append( *extra_pi0 );
      double Bmass_h = ( Dtaupi0_h.charge() ? PDG_BplusMass : PDG_B0Mass );
      HepLorentzVector Dtaupi0_h_4Vcm = Dtaupi0_h.p();
      Dtaupi0_h_4Vcm.boost( cmboost );
      double cosBDh_extra = ( 2.0*eb*Dtaupi0_h_4Vcm.e() - Bmass_h*Bmass_h-Dtaupi0_h.mass()*Dtaupi0_h.mass() ) / ( 2.0*sqrt(eb*eb-Bmass_h*Bmass_h)*Dtaupi0_h_4Vcm.vect().mag() );
      dist->column( epi_ch, cosBDh_extra     );
      Particle Dpi0_h( D_h.p() + extra_pi0->p(), D_h.pType() );	  
      dist->column( epi_mh, Dpi0_h.mass() );
    }else{
      dist->column( epi_pie, -10 );
      dist->column( epi_gle, -10 );
      dist->column( epi_ghe, -10 );
      dist->column( epi_cl,  100 );
      dist->column( epi_ml,  -10 );
      dist->column( epi_ch,  100 );
      dist->column( epi_mh,  -10 );
    }
    
    
    
    return;
  }
    

#if defined(BELLE_NAMESPACE)
}
#endif
  
