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
    runNoforsigMC(0),
    numberOfEvent(0),
    numberOfSkim(0),
    flag_DststMC(0),
    flag_single(0),
    flag_hadtag(0),
    flag_SkimFile(0),
    event_type(-1),
    selKPI         ( atc_pid(3, 1, 5, 3, 2) ),
    Dr_cut         ( 2.0 ), // cm
    Dz_cut         ( 5.0 ), //
    masscut_pi0_L  ( -0.030 ), // -0.030 ?
    masscut_pi0_H  (  0.030 ), // +0.030 ?
    masscut_D_L    ( -0.100 ), // -0.100, 0.045
    masscut_D_H    (  0.060 )  //  0.060, 0.045
  {
    std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    std::cout << "[ DSTRTAUNU's constructor ]" << std::endl;
    strcpy( SkimFileName, "Skim_DSTRTAUNU_test.index" );
    
    std::cout << "          Dr_cut = " << Dr_cut           << std::endl
	      << "          Dz_cut = " << Dz_cut           << std::endl
	      << "     masscut_D_L = " << masscut_D_L      << std::endl
	      << "     masscut_D_H = " << masscut_D_H      << std::endl
	      << "   masscut_pi0_L = " << masscut_pi0_L    << std::endl
	      << "   masscut_pi0_H = " << masscut_pi0_H    << std::endl
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
    dscr->define_param( "flag_DststMC",  "flag_DststMC",  &module->flag_DststMC     );
    dscr->define_param( "flag_single",   "flag_single",   &module->flag_single      );
    dscr->define_param( "flag_hadtag",   "flag_hadtag",   &module->flag_hadtag      );
    dscr->define_param( "flag_SkimFile", "flag_SkimFile", &module->flag_SkimFile    );
    dscr->define_param( "SkimFileName",  "SkimFileName",  256, module->SkimFileName );

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
    if( flag_SkimFile ) SkimFile = BASF_Output->open( SkimFileName );
    return;
  }

  // begin_run function
  void DSTRTAUNU::begin_run( BelleEvent* evptr, int* status )
  {
    //std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    //std::cout << "[ DSTRTAUNU's begin_run function ]" << std::endl;

    // IP information
    IpProfile::begin_run();
    m_IP     = IpProfile::position();
    m_IP_err = IpProfile::position_err();
    //std::cout << "    IpProfile::position() = " << m_IP
    //<< " +- " << m_IP_err
    //<< std::endl;

    // Beam Energy information
    BeamEnergy::begin_run();    
    E_HER = BeamEnergy::E_HER();
    E_LER = BeamEnergy::E_LER();
    //std::cout << "    BeamEnergy::E_HER()   = " << E_HER << std::endl
    //<< "    BeamEnergy::E_LER()   = " << E_LER << std::endl;
    X_ANGLE = 0.022; // rad
    lab = HepLorentzVector( E_HER * sin(X_ANGLE), 0.0,
			   E_HER * cos(X_ANGLE) - E_LER, E_HER + E_LER );
    cmboost = -lab.boostVector();
    HepLorentzVector cm = lab;

    cm.boost(cmboost);
    eb = cm.e()/2.0;

    //std::cout << "       4-momentum(lab)    = " << lab     << std::endl
    //<< "       cmboost vector     = " << cmboost << std::endl
    //<< "    half of CM energy(eb) = " << eb      << std::endl;

    // for eid
    eid::init_data();

    // for Ptype
    Ptype dummy( "VPHO" );

    if( flag_SkimFile ) SkimFile->write();

    runNoforsigMC++;

    //std::cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    return;
  }
  
  // hist_def function
  void DSTRTAUNU::hist_def()
  {

    std::cout << "[ DSTRTAUNU's hist_def function ]" << std::endl;
    extern BelleTupleManager* BASF_Histogram;
    BelleTupleManager* tm = BASF_Histogram;
    Skim_dist = tm->ntuple( "Skim_DSTRTAUNU", "exp run evt skim event", 1 );

    std::stringstream sTmp_gen;
    sTmp_gen << " exprun event " // [ EVENT ]
	     << " evt_type gm_dd gm_semi "
	     << " gm_df1hr gm_df1hnr gm_df1lep gm_df1nu gm_df1gam " // [ gen-MODE ]
	     << " gm_df2hr gm_df2hnr gm_df2lep gm_df2nu gm_df2gam "
	     << " gm_ds1hr gm_ds1hnr gm_ds1lep gm_ds1nu gm_ds1gam "
	     << " gm_ds2hr gm_ds2hnr gm_ds2lep gm_ds2nu gm_ds2gam "
	     << " gm_t1hr  gm_t1hnr  gm_t1lep  gm_t1nu  gm_t1gam gm_t1pro "
	     << " gm_t2hr  gm_t2hnr  gm_t2lep  gm_t2nu  gm_t2gam gm_t2pro "
	     << " nb1 nb1g nb1pip nb1pi0 semi1 "  // [B (gen) ]
	     << " nrootd1 rootdf1 rootds1 nd1 "   // [D (gen) ]
      	     << " nb2 nb1g nb2pip nb2pi0 semi2 "  // [B (gen) ]
	     << " nrootd2 rootdf2 rootds2 nd2 "   // [D (gen) ]
	     << " gm_ddst1 dst1_acc "             // [D*(gen) ]
	     << " gm_ddst2 dst2_acc "             // [D*(gen) ]
	     << " ddst1_d ddst1_acc "             // [D**(gen)]
	     << " ddst2_d ddst2_acc "             // [D**(gen)]
	     << " ddst1m ddst2m "                 // [D**(gen)]
	     << " rootd1nc rootd2nc "             // [D**(gen)]
	     << " fldstst1 fldstst2 "             // [D**(gen)]
	     << " incl1gmc incl2gmc "             // [D**(gen)]
	     << " nddst1pp nddst2pp "             // [D**(gen)]
	     << " nddst1p0 nddst2p0 "             // [D**(gen)]
	     << " flmcddst "                      // [D**(gen)]
	     << " ac1E  ac1cos  ac1ghE  ac1glE "  // [D*(gen) ]
	     << " ac2E  ac2cos  ac2ghE  ac2glE "  // [D*(gen) ]
	     << " ac1adp0  ac1adgh  ac1adgl "     // [D*(gen) ]
	     << " ac2adp0  ac2adgh  ac2adgl "     // [D*(gen) ]
	     << " ac11E ac11cos ac11ghE ac11glE " // [D**(gen)]
	     << " ac22E ac22cos ac22ghE ac22glE " // [D**(gen)]
	     << " ac11adp0  ac11adgh  ac11adgl "  // [D*(gen) ]
	     << " ac22adp0  ac22adgh  ac22adgl "  // [D*(gen) ]
	     << " pi0E pi0cos pi0ghE pi0glE "     // [pi0(gen) ]
	     << " taulep1 taulep2 "               // [tau(gen)] @20140925
	     << " cc1 cc2 "                       // [cc(gen) ]
	     << " k1_q2 k1_pl k1_w k1_cos "       // [kinematics] @20141113
	     << " k2_q2 k2_pl k2_w k2_cos "       // [kinematics] @20141113
	     << " ddst1m ddst2m "                 // [D**(gen)]
      ;

    Gen_dist = tm->ntuple( "gen",  sTmp_gen.str().c_str(), 16 );

    std::stringstream sTmp_rec;
    sTmp_rec << " exprun event evt_type "                     // [ number ]
	     << " rm_bb rm_dd rm_bfl "                        // [rec-mode BB]
	     << " rm_l1 rm_d1 rm_dst1 rm_d1lund rm_dst1lund " // [rec-mode B]
	     << " rm_l2 rm_d2 rm_dst2 rm_d2lund rm_dst2lund " //
	     << " dfl1 dfl2 dtchg1 dtchg2 "                   //
	     << " eecl neecl mmiss evis mmiss2 evis2 "          // [event varibales]
	     << " eeclth1 eeclth2 eeclth3 "                     // [ECL threshold study]
	     << " eecl_b1    eecl_b2    eecl_bg "               //
	     << " eecl_ef    eecl_b     eecl_eb "               //
	     << " eecl_neu   eecl_pi0   eecl_trk eecl_had "     //
	     << " eecl_0  eecl_1  eecl_2  eecl_3 "              //
	     << " eecl_4  eecl_5  eecl_6  eecl_7 "              //
	     << " eecl_8  eecl_9  eecl_10 eecl_11 "             //
	     << " eecl_12 eecl_13 eecl_14 eecl_15 "             //
	     << " eecl_16 eecl_17 eecl_18 eecl_19 "             //
	     << " eecl_20 eecl_21 eecl_22 eecl_23 "             //
	     << " eecl_24 eecl_25 eecl_26 eecl_27 "             //
	     << " eecl_28 eecl_29 eecl_30 eecl_31 "             //
	     << " eecl_32 "                                     //
	     << " neecl_b1    neecl_b2    neecl_bg "            //
	     << " neecl_ef    neecl_b     neecl_eb "            //
	     << " neecl_neu   neecl_pi0   neecl_trk neecl_had " //
	     << " neecl_0  neecl_1  neecl_2  neecl_3 "          //
	     << " neecl_4  neecl_5  neecl_6  neecl_7 "          //
	     << " neecl_8  neecl_9  neecl_10 neecl_11 "         //
	     << " neecl_12 neecl_13 neecl_14 neecl_15 "         //
	     << " neecl_16 neecl_17 neecl_18 neecl_19 "         //
	     << " neecl_20 neecl_21 neecl_22 neecl_23 "         //
	     << " neecl_24 neecl_25 neecl_26 neecl_27 "         //
	     << " neecl_28 neecl_29 neecl_30 neecl_31 "         //
	     << " neecl_32 "                                    //
	     << " remtrk rempi0 remks remgam "                  // [remaining particles]
	     << " rempi0_0 rempi0_1 rempi0_2 rempi0_3 "         // added @20140919
	     << " rempi0_4 rempi0_5 rempi0_6 rempi0_7 rempi0_8 "// 
	     << " l1pcm l1p l1pt l1pc "                         // [momentum]
	     << " l2pcm l2p l2pt l2pc "                         //
	     << " dst1pcm dst1p acc1p acc1pc d1pcm d1p "        //
	     << " dst2pcm dst2p acc2p acc2pc d2pcm d2p "        //
	     << " cosdl1 cosdl2 cosdlh cosdll rm_bdir "         // [Dtau]
	     << " dm1 dst1_morg dst1_m dst1self dst1org "       // [D*]
	     << " dm2 dst2_morg dst2_m dst2self dst2org "       //
	     << " dst1_pi0cos acc1self acc1org "                // [accompany particles]
	     << " dst2_pi0cos acc2self acc2org "                //
	     << " d1_morg d1_m d1self d1moid d1org "            // [D]
	     << " d2_morg d2_m d2self d2moid d2org "            //
	     << " l1self l1selfid l1moid l1pid l1org "          // [lepton]
	     << " l2self l2selfid l2moid l2pid l2org "          //
	     << " d1ch0recid d1ch0selfid d1ch0mo d1ch0org "     // [D's daughter]
	     << " d1ch1recid d1ch1selfid d1ch1mo d1ch1org "     //
	     << " d1ch2recid d1ch2selfid d1ch2mo d1ch2org "     //
	     << " d1ch3recid d1ch3selfid d1ch3mo d1ch3org "     //
	     << " d2ch0recid d2ch0selfid d2ch0mo d2ch0org "     //
	     << " d2ch1recid d2ch1selfid d2ch1mo d2ch1org "     //
	     << " d2ch2recid d2ch2selfid d2ch2mo d2ch2org "     //
	     << " d2ch3recid d2ch3selfid d2ch3mo d2ch3org "     // 
	     << " d1ch0p d1ch0pc "                              // @20141202
	     << " d1ch1p d1ch1pc "                              //
	     << " d1ch2p d1ch2pc "                              //
	     << " d1ch3p d1ch3pc "                              //
      	     << " d2ch0p d2ch0pc "                              // @20141202
	     << " d2ch1p d2ch1pc "                              //
	     << " d2ch2p d2ch2pc "                              //
	     << " d2ch3p d2ch3pc "                              //
	     << " k1pidmin k1pidmax "                           // [D's daughter  ]
	     << " k2pidmin k2pidmax "                           //
	     << " d1_pi0cos d2_pi0cos "                         // [D's daughter  ]
	     << " kfd1cl kfd1chi2 kfd1ndf "                      // [kinematic fit D]
	     << " kfd2cl kfd2chi2 kfd2ndf "                      //
	     << " kfs1cl kfs1chi2 kfs1ndf "                      // [kinematic fit D*]
	     << " kfs2cl kfs2chi2 kfs2ndf "                      //
      //<< " kfl1cl kfl1chi2 kfl1ndf "                    // [kinematic fit lepton]
      //<< " kfl2cl kfl2chi2 kfl2ndf "                    //
	     << " l1geg l1gng l1gnr l1gntr l1getr l1ged l1gnd "//[bremsstrahlung]
	     << " l2geg l2gng l2gnr l2gntr l2getr l2ged l2gnd "//
	     << " l1ag l2ag "                                 //
	     << " rgd1e rgd1a rgd1r "                         //
	     << " rgd2e rgd2a rgd2r "                         //
	     << " rgd3e rgd3a rgd3r "                         //
	     << " rgd4e rgd4a rgd4r "                         //
	     << " cosdla1 cosdla2 rm_bdira "                  //
	     << " cosdlb1 cosdlb2 rm_bdirb "                  //
	     << " cosdlc1 cosdlc2 rm_bdirc "                  //
	     << " cosdld1 cosdld2 rm_bdird "                  //
	     << " cosdle1 cosdle2 rm_bdire "                  //
	     << " rm_ddst_angle "                             //
	     << " mmiss_ctrl  evis_ctrl "                     // D** study
	     << " mmiss2_ctrl evis2_ctrl "                    // D** study
	     << " eecl_ctrl "                                 // D** study
	     << " cosdl1_ctrl cosdl2_ctrl"                    // D** study
	     << " cosdll_ctrl cosdlh_ctrl"                    // D** study
	     << " recgen "                                    // [gen-rec matching]
	     << " gm_dd gm_semi "                             // [ gen-mode (BB)  ]
	     << " nb1 nb1g nb1pip nb1pi0 semi1 "              // [ gen-mode (B )  ]
	     << " nrootd1 rootdf1 rootds1 nd1 "               //
	     << " gm_ddst1 dst1_d dst1_acc "                  //
	     << " ddst1_d ddst1_acc "                         // [ gen-mode (D**) ]
	     << " cc1 "                                       //
	     << " rootd1nc rootd2nc "                         // [D**(gen)]
	     << " fldstst1 fldstst2 "                         // [D**(gen)]
	     << " incl1gmc incl2gmc "                         // [D**(gen)]
	     << " nddst1pp nddst2pp "                         // [D**(gen)]
	     << " nddst1p0 nddst2p0 "                         // [D**(gen)]
	     << " flmcddst "                                  // [D**(gen)]
	     << " nb2 nb2g nb2pip nb2pi0 semi2 "              // [ gen-mode (B )  ]
	     << " nrootd2 rootdf2 rootds2 nd2 "               //
	     << " gm_ddst2 dst2_d dst2_acc "                  //
	     << " ddst2_d ddst2_acc "                         // [ gen-mode (D**) ]
	     << " cc2 "                                       //
	     << " acc1m acc2m pi01m pi02m "                   // [pi0 mass] 20141105
	     << " taulep1 taulep2 "                           // [tau(gen)] 20140925
	     << " gm_df1hr gm_df1hnr gm_df1lep gm_df1nu gm_df1gam "         // [ gen-mode (D)]
	     << " gm_df2hr gm_df2hnr gm_df2lep gm_df2nu gm_df2gam "         //
	     << " gm_ds1hr gm_ds1hnr gm_ds1lep gm_ds1nu gm_ds1gam "         //
	     << " gm_ds2hr gm_ds2hnr gm_ds2lep gm_ds2nu gm_ds2gam "         //
	     << " gm_t1hr  gm_t1hnr  gm_t1lep  gm_t1nu  gm_t1gam gm_t1pro " //
	     << " gm_t2hr  gm_t2hnr  gm_t2lep  gm_t2nu  gm_t2gam gm_t2pro " //
	     << " self "                                                    // [self] added @ 20140122, self is updated @ 20140730
	     << " k1_q2 k1_pl k1_w k1_cos "       // [kinematics] @20141113
	     << " k2_q2 k2_pl k2_w k2_cos "       // [kinematics] @20141113
      ;
    Rec_dist = tm->ntuple( "rec",  sTmp_rec.str().c_str(), 15 );
    
   
    std::stringstream sTmp_rec_single;
    sTmp_rec_single << " exprun event evt_type "                    // [ number ]
		    << " rm_bb "                                     // [rec-mode BB]
		    << " rm_l1 rm_d1 rm_dst1 rm_d1lund rm_dst1lund " // [rec-mode B]
		    << " dfl1 dtchg1 "                               //
		    << " l1pcm l1p l1pt l1pc "                       // [momentum]
		    << " dst1pcm dst1p acc1p acc1pc d1pcm d1p "      //
		    << " cosdl1 "                                    // [Dtau]
		    << " dm1 dst1_morg dst1_m dst1self dst1org "     // [D*]
		    << " dst1_pi0cos acc1self acc1org "              // [accompany particles]
		    << " d1_morg d1_m d1self d1moid d1org "          // [D]
		    << " l1self l1selfid l1moid l1pid l1org "        // [lepton]
		    << " nls1 nls9 nlo1 nlo9 "                       // [remaining lepton] @20141202
		    << " d1ch0recid d1ch0selfid d1ch0mo d1ch0org "   // [D's daughter]
		    << " d1ch1recid d1ch1selfid d1ch1mo d1ch1org "   //
		    << " d1ch2recid d1ch2selfid d1ch2mo d1ch2org "   //
		    << " d1ch3recid d1ch3selfid d1ch3mo d1ch3org "   //
		    << " d1ch0p d1ch0pc "                            // @20141202
		    << " d1ch1p d1ch1pc "                            //
		    << " d1ch2p d1ch2pc "                            //
		    << " d1ch3p d1ch3pc "                            //
		    << " k1pidmin k1pidmax "                         // [D's daughter  ]
		    << " d1_pi0cos "                                 // [D's daughter  ]
		    << " kfd1cl kfd1chi2 kfd1ndf "                    // [kinematic fit D]
		    << " kfs1cl kfs1chi2 kfs1ndf "                    // [kinematic fit D*]
      //<< " kfl1cl kfl1chi2 kfl1ndf "                    // [kinematic fit lepton]
		    << " l1geg l1gng l1gnr l1gntr l1getr l1ged l1gnd "//[bremsstrahlung]
		    << " l1ag "                                      //
		    << " rgd1e rgd1a rgd1r "                         //
		    << " rgd2e rgd2a rgd2r "                         //
		    << " rgd3e rgd3a rgd3r "                         //
		    << " rgd4e rgd4a rgd4r "                         //
		    << " cosdla1 cosdlb1 cosdlc1 cosdld1 cosdle1 "   //
		    << " recgen "                                    // [gen-rec matching]
		    << " gm_dd gm_semi "                             // [ gen-mode (BB)  ]
		    << " nb1 nb1g nb1pip nb1pi0 semi1 "              // [ gen-mode (B )  ]
		    << " nrootd1 rootdf1 rootds1 nd1 "               //
		    << " gm_ddst1 dst1_d dst1_acc "                  //
		    << " rootd1nc rootd2nc "                         // [D**(gen)]
		    << " fldstst1 fldstst2 "                         // [D**(gen)]
		    << " incl1gmc incl2gmc "                         // [D**(gen)]
		    << " nddst1pp nddst2pp "                         // [D**(gen)]
		    << " nddst1p0 nddst2p0 "                         // [D**(gen)]
		    << " flmcddst "                                  // [D**(gen)]
		    << " ddst1_d ddst1_acc "                         // [ gen-mode (D**) ]
		    << " cc1 "                                       //
		    << " nb2 nb2g semi2 "                            //
		    << " nrootd2 rootdf2 rootds2 nd2 "               //
		    << " gm_ddst2 dst2_acc "                         //
		    << " ddst2_d ddst2_acc "                         // [ gen-mode (D**) ]
		    << " cc2 "                                       //
		    << " acc1m pi01m "                               // [pi0 mass] 20141105
		    << " taulep1 taulep2 "                           // [tau(gen)] 20140925
		    << " gm_df1hr gm_df1hnr gm_df1lep gm_df1nu gm_df1gam "         // [ gen-mode (D)]
		    << " gm_df2hr gm_df2hnr gm_df2lep gm_df2nu gm_df2gam "         //
		    << " gm_ds1hr gm_ds1hnr gm_ds1lep gm_ds1nu gm_ds1gam "         //
		    << " gm_ds2hr gm_ds2hnr gm_ds2lep gm_ds2nu gm_ds2gam "         //
		    << " gm_t1hr  gm_t1hnr  gm_t1lep  gm_t1nu  gm_t1gam gm_t1pro " //
		    << " gm_t2hr  gm_t2hnr  gm_t2lep  gm_t2nu  gm_t2gam gm_t2pro " //
		    << " self old_self"                                            // [self] added @ 20140122
      ;
    Rec_single_dist = tm->ntuple( "rec_single",  sTmp_rec_single.str().c_str(), 17 );

    std::stringstream sTmp_rec_hadtag;
    sTmp_rec_hadtag << sTmp_rec_single.str().c_str()
		    << " tmbc tde "
		    << " tmcinfo tself "
		    << " tblund tbdecay tdecay1 tdecay2 tdecay3 tdecay4 tnfs " 
		    << " tnboutdef tbestdef tnboutcont tbestcont "
		    << " tdgb tdgdst1 tdgd1 tdgdst2 tdgd2 tdg "
		    << " tndst tnd tndsst tnds tnjpsi "
		    << " tdst1lund td1lund tdst2lund td2lund "
		    << " eecl neecl mmiss evis mmiss2 evis2 "
		    << " remtrk rempi0 remks remgam "
		    << " rempi0_0 rempi0_1 rempi0_2 rempi0_3 "         // added @20140919
		    << " rempi0_4 rempi0_5 rempi0_6 rempi0_7 rempi0_8 "// 
      ;
    Rec_hadtag_dist = tm->ntuple( "rec_hadtag",  sTmp_rec_hadtag.str().c_str(), 18 );

    //pi0_dist = tm->ntuple( "pi0", "pi0e pi0m ghe gle cos self", 111 );
    eecl_dist = tm->ntuple( "eecl", "exprun event evt_type eecl region origin threshold self neecl dtau1 dtau2", 22 );

    std::stringstream sTmp_hadtag;
    sTmp_hadtag << " exp run evt nbtag "
		<< " tmbc tde "
		<< " tmcinfo tself "
		<< " tblund tbdecay tdecay1 tdecay2 tdecay3 tdecay4 tnfs " 
		<< " tnboutdef tbestdef tnboutcont tbestcont "
		<< " tdgb tdgdst1 tdgd1 tdgdst2 tdgd2 tdg "
		<< " tndst tnd tndsst tnds tnjpsi "
		<< " tdst1lund td1lund tdst2lund td2lund "
      ;
    hadtag_dist = tm->ntuple( "hadtag", sTmp_hadtag.str().c_str(), 511 );
    
    
    return;
  }
  
  // event function
  void DSTRTAUNU::event( BelleEvent* evptr, int* status )
  {
    //std::cout << "[ DSTRTAUNU's event function ]" << std::endl;
    *status = 0;
    numberOfEvent++;
    
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
    if( flag_SkimFile ) delete SkimFile;
    return;
  }
      
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
