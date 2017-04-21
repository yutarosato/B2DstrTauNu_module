#include "DSTRTAUNU.h"

using namespace std;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  double DSTRTAUNU::calEclEnergyWithMatch2GoodGamma(const Particle& B1, const Particle& B2, const double width, const double e9e25)
  {
    //  cout<<"good_gamma: "<<width<<", "<<e9e25<<endl;
    Mdst_ecl_Manager& eclmgr = Mdst_ecl_Manager::get_manager();
    double eclenergy(0);
    for (Mdst_ecl_Manager::iterator i = eclmgr.begin(); i != eclmgr.end(); ++i){
      if( i->match()==1 || (i->match()==2&&!good_gamma(*i, 0.02, e9e25, width)) )continue; // check matching with charged tracks
      if( myCheckSame(B1,*i) )continue; // check if cluster is used for the reconstruction of B1
      if( myCheckSame(B2,*i) )continue; // check if cluster is used for the reconstruction of B2
      if(!checkEECLCluster(i->energy(), i->theta()))continue; // check if cluster energy exceed the threshold
      eclenergy+=i->energy();
    }
    return eclenergy;
  }


int DSTRTAUNU::checkEECLCluster(const double i_energy, const double i_theta, const double th_fw, const double th_br, const double th_bw)
{
  const double theta(i_theta*180./M_PI);
  if(theta<31.4){ // forward-endcap region
    if(i_energy>th_fw)return 1;
    else return 0;
  }
  else if(theta>130.7){ // backward-endcap region
    if(i_energy>th_bw)return 1;
    else return 0;
  }
  else{ // barrel region
    if(i_energy>th_br)return 1;
    else return 0;
  }
  return 0;
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

  double DSTRTAUNU::calBdirection( Particle& Dtau1, Particle& Dtau2, Hep3Vector& B_plus, Hep3Vector& B_minus ){

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
    

#if defined(BELLE_NAMESPACE)
}
#endif
  
