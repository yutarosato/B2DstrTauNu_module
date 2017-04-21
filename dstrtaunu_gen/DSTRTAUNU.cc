#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  //my analysis code
  // Constructor
  DSTRTAUNU::DSTRTAUNU() :
    mc(true),
    expNo(0),
    runNo(0),
    evtNo(0),
    event_type(-1),
    selKPI         ( atc_pid(3, 1, 5, 3, 2) ),
    selKPI_cut     ( 0.6 ),
    electron_P_cut ( 0.000 ),
    muon_P_cut     ( 0.000 ),
    Dr_cut         ( 2.0 ), // cm
    Dz_cut         ( 5.0 ), //
    eidprob_cut    ( 0.80 ), 
    muLH_cut       ( 0.80 ),
    masscut_pi0_L  ( -0.015 ),
    masscut_pi0_H  (  0.015 ),
    masscut_D_L    ( -0.050 ),
    masscut_D_H    (  0.050 ),
    delta_m_cut    (  0.180 ),
    gam_E_cut      ( 0.05 ),
    Rad_angle_cut  ( 0.05 )
  {
    std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    std::cout << "[ DSTRTAUNU's constructor ]" << std::endl;
    
    std::cout << "          Dr_cut = " << Dr_cut           << std::endl
	      << "          Dz_cut = " << Dz_cut           << std::endl
	      << "      selKPI_cut = " << selKPI_cut       << std::endl
	      << "  electron_P_cut = " << electron_P_cut   << std::endl
	      << "      muon_P_cut = " << muon_P_cut       << std::endl
	      << "     eidprob_cut = " << eidprob_cut      << std::endl
	      << "        muLH_cut = " << muLH_cut         << std::endl
	      << "     masscut_D_L = " << masscut_D_L      << std::endl
	      << "     masscut_D_H = " << masscut_D_H      << std::endl
	      << "   masscut_pi0_L = " << masscut_pi0_L    << std::endl
	      << "   masscut_pi0_H = " << masscut_pi0_H    << std::endl
	      << "       gam_E_cut = " << gam_E_cut        << std::endl
 	      << "   Rad_angle_cut = " << Rad_angle_cut    << std::endl
      ;

    recon_set.insert( Ks_LUND      );
    recon_set.insert( Kplus_LUND   );
    recon_set.insert( Kminus_LUND  );
    recon_set.insert( PI0_LUND     );
    recon_set.insert( PIplus_LUND  );
    recon_set.insert( PIminus_LUND );
    recon_set.insert( Electron_LUND    );
    recon_set.insert( Positron_LUND    );
    recon_set.insert( MUplus_LUND      );
    recon_set.insert( MUminus_LUND     );
    recon_set.insert( Kl_LUND          );
    recon_set.insert( Gamma_LUND       );
    recon_set.insert( Proton_LUND      );
    recon_set.insert( antiProton_LUND  );
    recon_set.insert( Neutron_LUND     );
    recon_set.insert( antiNeutron_LUND );
    recon_set.insert( Lambda_LUND      );
    recon_set.insert( antiLambda_LUND  );

    std::cout << "   recon_set : ";
    std::set<int>::iterator it_recon = recon_set.begin();
    while( it_recon != recon_set.end() ){
      std::cout << *it_recon << ", ";
      ++it_recon;
    }
    std::cout << std::endl;

    nonrecon_set.insert(     Nu_E_LUND    );
    nonrecon_set.insert( antiNu_E_LUND    );
    nonrecon_set.insert(     Nu_MU_LUND   );
    nonrecon_set.insert( antiNu_MU_LUND   );
    nonrecon_set.insert(     Nu_TAU_LUND  );
    nonrecon_set.insert( antiNu_TAU_LUND  );

    std::cout << "   nonrecon_set : ";
    std::set<int>::iterator it_nonrecon = nonrecon_set.begin();
    while( it_nonrecon != nonrecon_set.end() ){
      std::cout << *it_nonrecon << ", ";
      ++it_nonrecon;
    }
    std::cout << std::endl;


    D_set.insert( Dsplus_LUND    );
    D_set.insert( Dsminus_LUND   );
    D_set.insert( Dsstrp_LUND    );
    D_set.insert( Dsstrm_LUND    );
    D_set.insert( Dstr0_LUND     );
    D_set.insert( antiDstr0_LUND );
    D_set.insert( Dstrp_LUND     );
    D_set.insert( Dstrm_LUND     );
    D_set.insert( D0_LUND        );
    D_set.insert( antiD0_LUND    );
    D_set.insert( Dplus_LUND     );
    D_set.insert( Dminus_LUND    );
    for( int i=0; i<NhigherD_LUND; i++ ) D_set.insert( higherD_LUND[i] );
    
    std::cout << "   D_set : ";
    std::set<int>::iterator it_D = D_set.begin();
    while( it_D != D_set.end() ){
      std::cout << *it_D << ", ";
      ++it_D;
    }
    std::cout << std::endl;
    
    std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    
    return;
  }

  extern "C" Module_descr *mdcl_DSTRTAUNU()
  {
    DSTRTAUNU* module = new DSTRTAUNU;
    Module_descr* dscr = new Module_descr( "DSTRTAUNU", module );

    IpProfile::define_global( dscr ); // for IP tube constraint
    return dscr;
  }
  
  // Destructor
  DSTRTAUNU::~DSTRTAUNU( void )
  {
    return;
  }

  // initilization
  void DSTRTAUNU::init( int *status )
  {
    std::cout << "[ DSTRTAUNU's initilization ]" << std::endl;
    extern BasfOutputManager* BASF_Output;
    return;
  }

  // begin_run function
  void DSTRTAUNU::begin_run( BelleEvent* evptr, int* status )
  {
    std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    std::cout << "[ DSTRTAUNU's begin_run function ]" << std::endl;

    // IP information
    IpProfile::begin_run();
    m_IP     = IpProfile::position();
    m_IP_err = IpProfile::position_err();
    std::cout << "    IpProfile::position() = " << m_IP
	      << " +- " << m_IP_err
	      << std::endl;

    // Beam Energy information
    BeamEnergy::begin_run();    
    E_HER = BeamEnergy::E_HER();
    E_LER = BeamEnergy::E_LER();
    std::cout << "    BeamEnergy::E_HER()   = " << E_HER << std::endl
	      << "    BeamEnergy::E_LER()   = " << E_LER << std::endl;
    X_ANGLE = 0.022; // rad
    lab = HepLorentzVector( E_HER * sin(X_ANGLE), 0.0,
			   E_HER * cos(X_ANGLE) - E_LER, E_HER + E_LER );
    cmboost = -lab.boostVector();
    HepLorentzVector cm = lab;

    cm.boost(cmboost);
    eb = cm.e()/2.0;

    std::cout << "       4-momentum(lab)    = " << lab     << std::endl
	      << "       cmboost vector     = " << cmboost << std::endl
	      << "    half of CM energy(eb) = " << eb      << std::endl;

    // for eid
    eid::init_data();

    // for Ptype
    Ptype dummy( "VPHO" );

    std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    return;
  }
  
  // hist_def function
  void DSTRTAUNU::hist_def()
  {
    std::cout << "[ DSTRTAUNU's hist_def function ]" << std::endl;
    extern BelleTupleManager* BASF_Histogram;
    BelleTupleManager* tm = BASF_Histogram;

    std::stringstream sTmp_gen;
    sTmp_gen << " exprun event " // [ EVENT ]
	     << " evt_type gm_dd gm_semi "
	     << " gm_df1hr gm_df1hnr gm_df1lep gm_df1nu gm_df1gam " // [ gen-MODE ]
	     << " gm_df2hr gm_df2hnr gm_df2lep gm_df2nu gm_df2gam "
	     << " gm_ds1hr gm_ds1hnr gm_ds1lep gm_ds1nu gm_ds1gam "
	     << " gm_ds2hr gm_ds2hnr gm_ds2lep gm_ds2nu gm_ds2gam "
	     << " gm_t1hr  gm_t1hnr  gm_t1lep  gm_t1nu  gm_t1gam gm_t1pro "
	     << " gm_t2hr  gm_t2hnr  gm_t2lep  gm_t2nu  gm_t2gam gm_t2pro "
	     << " nb1 semi1 "                   // [B (gen) ]
	     << " nrootd1 rootdf1 rootds1 nd1 " // [D (gen) ]
	     << " gm_ddst1 dst1_acc acc1_E "    // [D*(gen) ]
	     << " acc1_ghE acc1_glE acc1_cos "  // [D*(gen) ]
	     << " cc1 "                         // [cc(gen) ]
	     << " nb2 semi2 "                   // [B (gen) ]
	     << " nrootd2 rootdf2 rootds2 nd2 " // [D (gen) ]
	     << " gm_ddst2 dst2_acc acc2_E "    // [D*(gen) ]
	     << " acc2_ghE acc2_glE acc2_cos "  // [D*(gen) ]
	     << " cc2 "                         // [cc(gen) ]
	     << " pi0_E pi0_cos pi0_ghE pi0_glE " // [pi0(gen) ]
      ;

    Gen_dist = tm->ntuple( "gen",  sTmp_gen.str().c_str(), 16 );

    std::stringstream sTmp_rec;
    sTmp_rec << " exprun event evt_type gm_dd gm_semi "   // [ EVENT ]
	     << " eecl mmiss evis "              // [var]
	     << " remtrk rempi0 remks remgam "   // [remain]
	     << " l1pcm l2pcm d1pcm d2pcm "      // [momentum]
	     << " l1p   l2p   d1p   d2p "        //
	     << " l1pt  l2pt "                   //
	     << " l1pid l2pid "                  // [pid]
	     << " rm_bb rm_dd rm_bfl "           // [rec-MODE]
	     << " rm_l1 rm_d1 rm_dst1 "          //
	     << " cosdl1 "                       // [Dtau]
	     << " dm1 dst1_m "                   // [D*]
	     << " d1_m "                         // [D]
	     << " rm_l2 rm_d2 rm_dst2 "          // [rec-MODE]
	     << " cosdl2 "                       // [Dtau]
	     << " dm2 dst2_m "                   // [D*]
	     << " d2_m "                         // [D]
	     << " d1self dst1self " // [true/fake]
	     << " d2self dst2self "
	     << " rmdmo1 rmdmo2 "
	     << " l1self l2self "
	     << " recgen "
	     << " nb1 semi1 "                   // [B (gen) ]
	     << " nrootd1 rootdf1 rootds1 nd1 " // [D (gen) ]
	     << " gm_ddst1 dst1_acc "           // [D*(gen) ]
	     << " cc1 "                         // [cc(gen) ]
	     << " nb2 semi2 "                   // [B(gen)  ]
	     << " nrootd2 rootdf2 rootds2 nd2 " // [D(gen)  ]
	     << " gm_ddst2 dst2_acc "           // [D*(gen) ]
	     << " cc2 "                         // [cc(gen) ]
	     << " rm_bdir " // [B direction]
	     << " gm_df1hr gm_df1hnr gm_df1lep gm_df1nu gm_df1gam " // [ gen-MODE ]
	     << " gm_df2hr gm_df2hnr gm_df2lep gm_df2nu gm_df2gam "
	     << " gm_ds1hr gm_ds1hnr gm_ds1lep gm_ds1nu gm_ds1gam "
	     << " gm_ds2hr gm_ds2hnr gm_ds2lep gm_ds2nu gm_ds2gam "
	     << " gm_t1hr  gm_t1hnr  gm_t1lep  gm_t1nu  gm_t1gam gm_t1pro "
	     << " gm_t2hr  gm_t2hnr  gm_t2lep  gm_t2nu  gm_t2gam gm_t2pro "
      ;
    Rec_dist = tm->ntuple( "rec",  sTmp_rec.str().c_str(), 15 );


    return;
  }
  
  // event function
  void DSTRTAUNU::event( BelleEvent* evptr, int* status )
  {
    //std::cout << "[ DSTRTAUNU's event function ]" << std::endl;
    *status = 0;
    
    // Belle event info
    expNo = 0;
    runNo = 0;
    evtNo = 0;
    
    getEventInfo(expNo, runNo, evtNo, mc);

    int expmc = ((belle_event*)BsGetEnt(BELLE_EVENT,1,BBS_No_Index))->m_ExpMC; // for IP tube constraint
    if( expmc!=2 ){
      if( !IpProfile::usable() ) return;
      if( !IpProfile::b_life_smeared_usable() ) return;
    }
    *status = event_start();
    return;
  }
  
  // end_run function
  void DSTRTAUNU::end_run( BelleEvent* evtptr, int* status )
  {
    std::cout << "[ DSTRTAUNU's end_run function ]" << std::endl;
    return;
  }

  // term function
  void DSTRTAUNU::term()
  {
    std::cout << "[ DSTRTAUNU's term function ]" << std::endl;
    return;
  }
      
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
