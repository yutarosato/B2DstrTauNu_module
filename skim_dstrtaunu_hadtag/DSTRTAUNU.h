#ifndef DSTRTAUNU_H
#define DSTRTAUNU_H

#include "belle.h"
#include "eid/eid.h"
#include "kid/atc_pid.h"
#include "mdst/Muid_mdst.h"

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"

#include "basf/module.h"
#include "basf/module_descr.h"
#include "basf/basfout.h"

#include "panther/panther.h"

#include "particle/Particle.h"
#include "particle/utility.h"

#include "toolbox/FuncPtr.h"
#include "FoxWolfr.h"
#include "ksfwmoments.h"

#include "toolbox/Thrust.h"

#include BRECON_H
#include "tables/ekpfullrecon_panther.h"
#include "fullrecon/frec_util.h"
#include "hamlet/AnaBrecon.h"

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include MDST_H
#include BELLETDF_H
#include EVTCLS_H

#include "benergy/BeamEnergy.h"
#include "ip/IpProfile.h"

#include "mdst/mdst.h"
#include "mdst/findKs.h"

#include "exkfitter/ExKFitter.h"
#include "tables/hepevt.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "UserInfo.h"

class BelleEvent;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  
  // define parameter 
  // CODE
  const int EE_CODE     = 0;
  const int MUON_CODE   = 1;
  const int PION_CODE   = 2;
  const int KAON_CODE   = 3;
  const int PROTON_CODE = 4;
  // Mass
  // B
  const double PDG_B0Mass     =  5.27953;
  const double PDG_BplusMass  =  5.27915;
  const double PDG_BminusMass =  PDG_BplusMass;
  // D*
  const double PDG_D0strMass =  2.00696;
  const double PDG_DpstrMass =  2.01025;
  const double PDG_DmstrMass =  PDG_DpstrMass;
  // D
  const double PDG_D0Mass     =  1.86484;
  const double PDG_DplusMass  =  1.86962;
  const double PDG_DminusMass =  PDG_DplusMass;
  // other
  const double PDG_JPsiMass   =  3.096916;
  const double PDG_KsMass     =  0.497614;
  const double PDG_Kstr0Mass  =  0.89594;
  const double PDG_KstrpMass  =  0.89166;
  const double PDG_KstrmMass  =  PDG_KstrpMass;
  const double PDG_PIplusMass =  0.13957018;
  const double PDG_PI0Mass    =  0.1349766;
  const double PDG_MUMass     =  0.105658367;
  // Lund
  const int Upsilon5S_LUND = 9000553;
  const int Upsilon4S_LUND = 300553;
  // B
  const int B0_LUND        =  511;
  const int B0str_LUND     =  513;
  const int antiB0_LUND    = -511;
  const int antiB0str_LUND = -513;
  const int Bplus_LUND     =  521;
  const int Bpstr_LUND     =  523;
  const int Bminus_LUND    = -521;
  const int Bmstr_LUND     = -523;
  // Xs
  const int Xsd_LUND       =  30343;
  const int antiXsd_LUND   = -30343;
  const int Xsu_LUND       =  30353;
  const int antiXsu_LUND   = -30353;
  const int Xss_LUND       =  30363;
  const int antiXss_LUND   = -30363;
  // Ds
  const int Dsplus_LUND  =  431;
  const int Dsminus_LUND = -431;
  const int Dsstrp_LUND  =  433;
  const int Dsstrm_LUND  = -433;
  // D*
  const int Dstr0_LUND     =  423;
  const int antiDstr0_LUND = -423;
  const int Dstrp_LUND     =  413;
  const int Dstrm_LUND     = -413;
  // D
  const int D0_LUND     =  421;
  const int antiD0_LUND = -421;
  const int Dplus_LUND  =  411;
  const int Dminus_LUND = -411;
  // higher D
  const int NhigherD_LUND = 34;
  const int higherD_LUND[NhigherD_LUND] = { 10411, 10421,  100411,  100421,  10413,  10423,  20413,  20423,  100413,  100423,  415,  425,  10431,  10433,  20433,  435,  9000433,
					   -10411,-10421, -100411, -100421, -10413, -10423, -20413, -20423, -100413, -100423, -415, -425, -10431, -10433, -20433, -435, -9000433,
  }; // (four states with L=1[D0*,D1,D1',D2*]) + (two radial excited states[D(2S),D*(2S)]) + (five Ds states)

  // charmonium
  const int JPsi_LUND     =     443;
  const int Psi2S_LUND    =  100443;
  const int Psi3770_LUND  =   30443;
  const int Psi4040_LUND  = 9000443; // [bug LUND code] psi(4040)=-7776773
  const int Psi4160_LUND  = 9010443; // [bug LUND code] psi(4160)=-7766773
  const int Psi4415_LUND  = 9020443;
  // Lepton + gamma
  const int Electron_LUND   =   11;
  const int Positron_LUND   =  -11; 
  const int Nu_E_LUND       =   12;
  const int antiNu_E_LUND   =  -12; 
  const int MUplus_LUND     =  -13;
  const int MUminus_LUND    =   13;
  const int Nu_MU_LUND      =   14;
  const int antiNu_MU_LUND  =  -14; 
  const int TAUplus_LUND    =  -15;
  const int TAUminus_LUND   =   15;
  const int Nu_TAU_LUND     =   16;
  const int antiNu_TAU_LUND =  -16; 
  const int Gamma_LUND      =   22;
  // Kaon
  const int Kstr0_LUND     =  313;
  const int antiKstr0_LUND = -313;
  const int Kstrp_LUND     =  323;
  const int Kstrm_LUND     = -323;
  const int K0_LUND        =  311;
  const int antiK0_LUND    = -311;
  const int Ks_LUND        =  310;
  const int Kl_LUND        =  130;
  const int Kplus_LUND     =  321;
  const int Kminus_LUND    = -321;
  // Pion
  const int PI0_LUND     =  111;
  const int PIplus_LUND  =  211;
  const int PIminus_LUND = -211;
  // Light Mesons
  const int Eta_LUND     =  221;
  const int antiEta_LUND = -221;
  // Baryons
  const int Proton_LUND      =  2212;
  const int antiProton_LUND  = -2212;
  const int Neutron_LUND     =  2112;
  const int antiNeutron_LUND = -2112;
  const int Lambda_LUND      =  3122;
  const int antiLambda_LUND  = -3122;

  class DSTRTAUNU : public Module {
  public:
    DSTRTAUNU( void );
    ~DSTRTAUNU( void );
    void init( int* );
    void begin_run( BelleEvent*, int* );
    void disp_stat( const char* ){};
    void hist_def( void );
    void event( BelleEvent*, int* );
    void end_run( BelleEvent*, int* );
    void term( void );

    bool mc;
    int  expNo;
    int  runNo;
    int  evtNo;
    int  runNoforsigMC;

    int numberOfEvent;
    int numberOfSkim;

    std::set<int> recon_set;    // reconstruction particles(Ks, K+, K-, pi0, pi+, pi-, e+, e-, mu+, mu-, KL, gamma, neutron, proton )
    std::set<int> nonrecon_set; // not reconstruction particles but final particles (neutrinos)
    std::set<int> D_set;        // D particles

    HepPoint3D   m_IP;
    HepSymMatrix m_IP_err;

    // module parameters
    int  flag_DststMC;
    int  flag_single;
    int  flag_hadtag;
    int  flag_SkimFile;
    char SkimFileName[256];

    // pid
    atc_pid selKPI;
    
    // cut
    double Dr_cut;
    double Dz_cut;
    double masscut_D_L;
    double masscut_D_H;
    double masscut_pi0_L;
    double masscut_pi0_H;

    // gen
    int event_type;

  private:
    double E_LER, E_HER, X_ANGLE;
    HepLorentzVector lab; // 4-momentum in lab-frame
    Hep3Vector cmboost; //
    double eb; // the half of CM energy = B energy

    // ecl cluster lists
    std::vector<Hep3Vector> m_ecls;
    std::vector<int> m_ecl_ids;

    // histogram
    BelleTuple* Rec_dist;
    BelleTuple* Rec_single_dist;
    BelleTuple* Rec_hadtag_dist;
    BelleTuple* Gen_dist;
    //BelleTuple* pi0_dist;
    BelleTuple* eecl_dist;
    BelleTuple* hadtag_dist;

    BelleTuple* Skim_dist;
    BasfOutput* SkimFile;
        
  protected:
    int  event_start();                                          // DSTRTAUNU_event.cc
    void Skim_ana(); // DSTRTAUNU_event.cc
    
    void Gamma_cand ( std::vector<Particle>& gamma_list,                                                                      const bool fl_message=0 ); // Gamma.cc
    void EE_cand    ( std::vector<Particle>& lep_list,  double pid_cut=-10.0,  double mom_cut=-10.0,                          const bool fl_message=0 ); // EE.cc
    void Mu_cand    ( std::vector<Particle>& lep_list,  double pid_cut=-10.0,  double mom_cut=-10.0,                          const bool fl_message=0 ); // Mu.cc
    void Kaon_cand  ( std::vector<Particle>& kaon_list, double pid_cut=-10.0,                                                 const bool fl_message=0 ); // Kaon.cc    
    void Pion_cand  ( std::vector<Particle>& pion_list,                                                                       const bool fl_message=0 ); // Pion.cc
    void Tracks_cand( std::vector<Particle>& pion_list,                                                                       const bool fl_message=0 ); // Tracks.cc    
    void Ks_cand    ( std::vector<Particle>& ks_list,   int   flag_kfitter=0,                                                 const bool fl_message=0 ); // Ks.cc
    void Pi0_cand   ( std::vector<Particle>& pi0_list,  double masscut_L=-0.010, double masscut_H=0.010, double cos_cut=-10.0,
		      double gam_e_cut_L=0.050, double gam_e_cut_H=0.050,double pi0_e_cut=-10.0, const bool fl_message=0 ); // Pi0.cc

    void Rec_D_2body( std::vector<Particle>& d_list,    std::vector<Particle>& k_list, std::vector<Particle>& pi_list,  const bool fl_message    =0 ); // Rec_D_2body.cc
    void Rec_D_2body( std::vector<Particle>& d_list,    std::vector<Particle>& p_list,                                  const bool fl_message    =0 ); // Rec_D_2body.cc

    void Rec_D_3body( std::vector<Particle>& d_list,    std::vector<Particle>& k_list, std::vector<Particle>& p_list, std::vector<Particle>& pi0_list,  const bool fl_message    =0 ); // Rec_D_3body.cc
    void Rec_D_3body( std::vector<Particle>& d_list,    std::vector<Particle>& k_list, std::vector<Particle>& p_list,                                   const bool fl_message    =0 ); // Rec_D_3body.cc

    void Rec_D_4body( std::vector<Particle>& d_list,    std::vector<Particle>& k_list, std::vector<Particle>& p_list, std::vector<Particle>& pi0_list,  const bool fl_message    =0 ); // Rec_D_4body.cc
    void Rec_D_4body( std::vector<Particle>& d_list,    std::vector<Particle>& k_list, std::vector<Particle>& p_list,                                   const bool fl_message    =0 ); // Rec_D_4body.cc


    void Rec_Dstr ( std::vector<Particle>& dstr_list, std::vector<Particle>& d_list, std::vector<Particle>& pi_list,  double delta_m_cut=0.010, const bool fl_message    =0 ); // Rec_Dstr.cc
    void Rec_Dtau ( std::vector<Particle>& dtau_list, std::vector<Particle>& d_list, std::vector<Particle>& lep_list, std::vector<Particle>& gamma_list, const bool fl_message    =0 ); // Rec_Dtau.cc
    int FullReconN( const int lund=B0_LUND, const double nb_th=0.005, const bool fl_nb_cont=0 );
    void FullRecon( bool fl_message=false, const bool fl_dump=false ); // FullRecon.cc
    int FullReconMode( const int decay, const int decay1, const int decay2, const int decay3, const int decay4,
		       int& digit_B, int digit_D[4], int rootlund[4], int nroot[5]); // FullRecon.cc
    int FullReconMode_B   ( const int decay, int nroot[5] );
    int FullReconMode_Dstr( const int decay, int& lund );
    int FullReconMode_D   ( const int decay, int& lund );

    int  Ana      ( std::vector<Particle>& dtau_list, int gen_b_decay_info[][40], int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6],
		    std::vector<Particle>& trk_list, std::vector<Particle>& gam_list, std::vector<Particle>& pi0_list, std::vector<Particle>& ks_list,
		    std::vector<Particle>& pi0_list_test0, std::vector<Particle>& pi0_list_test1, std::vector<Particle>& pi0_list_test2, std::vector<Particle>& pi0_list_test3,
		    std::vector<Particle>& pi0_list_test4, std::vector<Particle>& pi0_list_test5, std::vector<Particle>& pi0_list_test6, std::vector<Particle>& pi0_list_test7, std::vector<Particle>& pi0_list_test8,
		    BelleTuple* dist, const bool fl_message_ana=0 ); // Ana.cc

    int  Ana_hadtag( std::vector<Particle>& dtau_list, int gen_b_decay_info[][40], int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6],
		     std::vector<Particle>& trk_list, std::vector<Particle>& gam_list, std::vector<Particle>& lep_list,
		     std::vector<Particle>& pi0_list, std::vector<Particle>& ks_list,
		     std::vector<Particle>& pi0_list_test0, std::vector<Particle>& pi0_list_test1, std::vector<Particle>& pi0_list_test2, std::vector<Particle>& pi0_list_test3,
		     std::vector<Particle>& pi0_list_test4, std::vector<Particle>& pi0_list_test5, std::vector<Particle>& pi0_list_test6, std::vector<Particle>& pi0_list_test7, std::vector<Particle>& pi0_list_test8,
		     BelleTuple* dist, const bool fl_message_ana=0 ); // Ana_hadtag.cc

    int Ana_single( std::vector<Particle>& dtau_list, int gen_b_decay_info[][40], int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], std::vector<Particle>& gam_list, std::vector<Particle>& lep_list, 
		    BelleTuple* dist, const bool fl_dump=true, const bool fl_message_ana=0 ); // Ana_single.cc

    void Ana_single_sub( std::vector<Particle>::iterator dtau1, int gen_b_decay_info[][40], int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], std::vector<Particle>& gam_list, std::vector<Particle>& lep_list,
			 BelleTuple* dist, const bool fl_dump=true, const bool fl_message_ana=0 ); // Ana_single.cc
    
    void Rec_CC ( std::vector<Particle>& cc_list, std::vector<Particle>& lep_list, std::vector<Particle>& gamma_list, const bool fl_message=0 ); // Rec_CC.cc

    double correct_dr( const Mdst_charged chg, HepPoint3D m_ip, int id );  // Util.cc
    double correct_dz( const Mdst_charged chg, HepPoint3D m_ip, int id );  // Util.cc
    
    int add_bremsstrahlung( Particle &particle, std::vector<Particle>& gamma_list,
			    double angle_cut=0.05, bool fl_write=true,
			    int fl_par1=0, int fl_par2=0, int fl_par3=0 );     // Util.cc
    
    int check_dupli_daughter( std::vector<Particle>::iterator i, std::vector<Particle>::iterator j ); // Util.cc
    int check_dupli_daughter(             Particle            i,             Particle            j ); // Util.cc
    int kinematic_fit( Particle& particle, int flag_kfitter );                 // Util.cc
    void qq_suppress( Particle b, double& thrust_angle, double& R2 );          // Util.cc
    void cal_Mmiss_Evis( BelleTuple* dist, std::vector<Particle>& trk_list, std::vector<Particle>& gamma_list, const char* mmiss_name="mmiss",  const char* evis_name="evis",  Particle* veto_pi0=NULL ); // Util.cc
    void cal_Mmiss_Evis( BelleTuple* dist, Particle par1,                   Particle par2                    , const char* mmiss_name="mmiss2", const char* evis_name="evis2", Particle* veto_pi0=NULL ); // Util.cc

    // [ Combinatorial B.G. Suppression ]
    
    void setUserInfo( Particle& particle ){
      if(  &(particle.userInfo())==NULL ) particle.userInfo( UserInfo() );
    }

    int masscut( Particle particle, double low, double high, double mass ); // Util.cc
    double cos2track( Particle par1, Particle par2 );
    
    int check_selfF( Particle& particle );                         // Util_Gen.cc
    int check_selfF( Particle& particle, int self_LUND,
		     int mother_LUND=0, int mothermother_LUND=0 ); // Util_Gen.cc
    int check_idF( Particle& particle );                           // Util_Gen.cc
    int check_motheridF( Particle& particle );                     // Util_Gen.cc
    int check_selfG( Particle& particle );                         // Util_Gen.cc
    int check_selfR ( Particle& particle, int specified_nchild=0 ); // Util_Gen.cc
    int check_selfR2( Particle& particle                         ); // Util_Gen.cc

    void display_particle( Particle particle );                              // Util_Gen.cc
    void display_hepevt  ( Gen_hepevt gen );                                 // Util_Gen.cc
    void display_list    ( std::set<int> set, std::multiset<int> multiset ); // Util_Gen.cc
    void rec_message     ( std::vector<Particle>::iterator p              ); // Util_Gen.cc
    void display_rec_particle( Particle particle, const bool fl_message = false, int indent=1 ); // Util_Gen.cc
    
    void find_fin_child( Gen_hepevt gen,
			 std::map<int, int>& child_id_map,
			 std::multiset<int>& n_particle_set,
			 bool fl_message = false,
			 int indent = 1 );                                  // Util_Gen.cc
    
    void Gen_info( std::map<int, int> child_id_map_d  [][2], std::multiset<int> n_particle_set_d  [][2],
		   std::map<int, int> child_id_map_tau[][2], std::multiset<int> n_particle_set_tau[][2],
		   int gen_b_decay_info[][40],
		   int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
		   const int fl_message=0, const bool fl_dump=false ); // Gen_info.cc
    
    int mapping_delete( Particle particle, std::map<int, int>& child_id_map,
			 int& n_recon, int gen_mode_bg,
			 char* name, const bool fl_message=false );  // Util_Gen.cc
    
    int search_daughter_D( Gen_hepevt gen, int& id1, int& id2, bool fl_message=false, double threshold=1.8 );

    int search_daughter_D_sub( Gen_hepevt gen, std::vector<int>& d_id_list,
			       bool fl_message, int indent=1, double threshold=1.8 );

    int search_daughter_Dstr( Gen_hepevt gen, int& id1, int& id2, bool fl_message=false, double threshold=1.8 );

    int search_daughter_Dstr_sub( Gen_hepevt gen, std::vector<int>& dstr_id_list,
				  bool fl_message, int indent=1, double threshold=1.8 );

    int search_daughter_charm( Gen_hepevt gen, int& id1, int& id2, bool fl_message=false, double threshold=1.8 );

    int search_daughter_charm_sub( Gen_hepevt gen, std::vector<int>& d_id_list,
				   bool fl_message, int indent=1, double threshold=1.8 );

    int search_origin( Gen_hepevt gen, int lund );

    int semilept_dec( Gen_hepevt gen, int& id, bool fl_message=false );

    int search_prompt_lepton( Gen_hepevt gen, int& lund, int& id, bool fl_message=false );

    int search_daughter_CC( Gen_hepevt gen, bool fl_message=false, int indent=1,
			    double threshold=3.0 ); // Util_Gen.cc

    int make_mode_digit( std::multiset<int>& n_particle_set, int digit[] ); // Util_Gen.cc
    int make_flag_semi( int fl1, int fl2 ); // Util_Gen.cc
    int make_flag_dd  ( int fl1, int fl2 ); // Util_Gen.cc

    int which_B( Gen_hepevt gen, int id1, int id2=-1 ); // Util_Gen.cc

    void    ckEclEnergyWithMatch2GoodGamma( Particle& B1, Particle& B2, const int fl_self, const int neecl, const bool fl_message=false,
					    const double width=6.0,     const double e9e25=0.94,
					    const double th_fw=0.10,    const double th_br=0.05, const double th_bw=0.15
					    ); // Util_dstrtaunu.cc
    double calEclEnergyWithMatch2GoodGamma( const Particle& B1,     const Particle& B2, const bool fl_message=false,
					    double eecl[39]=NULL, int neecl[39]=NULL,
					    const Particle* veto_pi0=NULL,
					    const int id1=0, const int id2=0,
					    const double width=6.0, const double e9e25=0.94,
					    const double th_fw=0.10,    const double th_br=0.05, const double th_bw=0.15
					    ); // Util_dstrtaunu.cc
    int checkEECLCluster( const double i_energy, const double i_theta, const double th_fw=0.10, const double th_br=0.05, const double th_bw=0.15); // Util_dstrtaunu.cc
    int checkEECLClusterOrigin( const Gen_hepevt& hep, bool fl_message=false ); // Util_dstrtaunu.cc
    int checkEECLClusterOrigin_KsLambda( const Gen_hepevt& hep ); // Util_dstrtaunu.cc
    int cnt_remain_trk  ( const Particle& B1, const Particle& B2, double dr=1.0, double dz=5.0      ); // Util_dstrtaunu.cc
    int cnt_remain_pi0  ( const Particle& B1, const Particle& B2, std::vector<Particle>& pi0_list   ); // Util_dstrtaunu.cc
    int cnt_remain_ks   ( const Particle& B1, const Particle& B2, std::vector<Particle>& ks_list    ); // Util_dstrtaunu.cc
    int cnt_remain_gamma( const Particle& B1, const Particle& B2, std::vector<Particle>& gamma_list ); // Util_dstrtaunu.cc

    double cal_nearestgammadirection( const Particle par1, const Particle par2, const Particle& B1, const Particle& B2, std::vector<Particle>& gam_list );

    double calBdirection( Particle& Dtau1, Particle& Dtau2 ); // Util_dstrtaunu.cc
    double calBdirection_with_more_bremsstrahlung( Particle Dtau1, Particle Dtau2, std::vector<Particle>& gamma_list,
						   double angle_cut,
						   int fl_par1=0, int fl_par2=0, int fl_par3=0 ); // Util_dstrtaunu.cc
    double cosBDl_with_more_bremsstrahlung( Particle Dtau1, std::vector<Particle>& gamma_list,
					    double angle_cut,
					    int fl_par1=0, int fl_par2=0, int fl_par3=0 ); // Util_dstrtaunu.cc
    
    bool myCheckSame(const Particle&, const Mdst_charged&);
    bool myCheckSame(const Particle&, const Mdst_gamma&);
    bool myCheckSame(const Particle&, const Mdst_pi0&);
    bool myCheckSame(const Particle&, const Mdst_ecl&);
    bool myCheckSame(const std::vector<Particle>& plist,const Mdst_charged&);

    double radgam_energy( Gen_hepevt gen );
    int    radgam_cnt   ( Gen_hepevt gen );
    void calKinematics  ( int id_b, int id_d, int id_lep,
			  double& q2, double& plep, double& w, double& coslep );
  };
#if defined(BELLE_NAMESPACE)
}
#endif

#endif
