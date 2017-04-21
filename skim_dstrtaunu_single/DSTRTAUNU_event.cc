#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  int DSTRTAUNU::event_start()
  {
    std::vector<Particle> trk_list;
    std::vector<Particle> kaon_list;
    std::vector<Particle> ks_list;
    std::vector<Particle> lep_list;
    std::vector<Particle> gamma_list;
    std::vector<Particle> pion_list;
    std::vector<Particle> pi0_list;
    std::vector<Particle> pi0_for_dstr_list;

    std::vector<Particle> recd0_list;
    std::vector<Particle> recdp_list;
    std::vector<Particle> recdstr0_list;
    std::vector<Particle> recdstrp_list;
    std::vector<Particle> recdtau_list;
    
    std::map<int,int>  child_id_map_D    [2][2]; // < hep.get_ID(), hep.idhep() >
    std::multiset<int> n_particle_set_D  [2][2]; // < LUND >
    std::map<int,int>  child_id_map_tau  [2][2]; // < hep.get_ID(), hep.idhep() >
    std::multiset<int> n_particle_set_tau[2][2]; // < LUND >

    int gen_b_decay_info[2][40]={{0},{0}};
    // ( [0][] is signal-B, [1][] is normalization-B, if there is signal events. )
    // [][ 0] B->Get_ID()
    // [][ 1] B semi-decay : 0(no lepton), 1(e-nu), 2(mu-nu), 3(tau-nu), -1(other two leptons), -2(others)
    // [][ 2] lepton->Get_ID(), lepton from semileptonic-B decay
    // [][ 3] # of B's (direct) children
    // [][ 4] # of gamma directly from B
    // [][29] # of pi+-  directly from B
    // [][30] # of pi0   directly from B
    // [][ 5] # of rootD 'directly' from B
    // [][ 6] 1st-rootD->idhep()
    // [][ 7] 1st-rootD->get_ID()
    // [][ 8] 2nd-rootD->idhep()
    // [][ 9] 2nd-rootD->get_ID()
    // [][10] # of D in B decay, where D is [D0/D+/Ds]
    // [][11] 1st-D->Get_ID()
    // [][12] 2nd-D->Get_ID()
    // [][13] flag of D/D*/D**/DD, 1(D), 2(D*), 3(D**), 4(two D), 5(no D)
    // [][14] accompany particle in D* decay->idhep()
    // [][15] accompany particle in D* decay->get_ID()
    // [][27] daughter  particle in D* decay->idhep()  @ 20141002
    // [][28] daughter  particle in D* decay->get_ID() @ 20141002

    // [][23] daughter  particle in rootD(D**) decay->idhep()  for D** study @ 20140703
    // [][24] daughter  particle in rootD(D**) decay->get_ID() for D** study @ 20140703
    // [][25] accompany particle in rootD(D**) decay->idhep()  for D** study @ 20140703
    // [][26] accompany particle in rootD(D**) decay->get_ID() for D** study @ 20140703

    // [][31] # of rootD's(D**) children @20141028
    // [][32] flag of D** : 0(not D**), 1(L=1), 2(radial excitation), 3(other D**?) @20141028
    // [][33] flag of inclusion in gMC : 1(included [D**lnu MC]), 0(not included [D**lnu MC]), 0([gMC]) @20141028
    // [][34] # of pi+- in rootD(D**) decay @20141028
    // [][35] # of pi0  in rootD(D**) decay @20141028
    
    
    // [][16] LUND code of charmonium from B

    // [][17] 1stD semi-decay
    // [][18] lepton->Get_ID(), lepton from 1stD
    // [][19] 2ndD semi-decay
    // [][20] lepton->Get_ID(), lepton from 2ndD

    // [][21] prompt lepton from tau->idhep()
    // [][22] prompt lepton from tau->Get_ID()
    
    
    int gen_d_mode_info  [2][2][6] = {{{0}},{{0}}};
    int gen_tau_mode_info[2][6]    = {{0},{0}};

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const int  fl_message_gen    = 0; // 0(no-message),1(gen),2(gen+table)
    const bool fl_message_funda  = !true;
    const bool fl_message_rec    = !true;
    const bool fl_message_ana    = !true;
    const bool fl_message_hadtag = !true;

    const bool fl_gen_dump      = !true;
    const bool fl_hadtag_dump   = !true;
    //flag_single = true;
    //flag_hadtag = true;
    //if( !(runNo==1507 && evtNo<1000) ) return 0;
    //if( evtNo!=3264) return 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message_gen || fl_message_hadtag || fl_message_funda || fl_message_rec || fl_message_ana )
      std::cout << "*******************************"
		<< " [ DSTRTAUNU ] "
		<< "expNo = "   << expNo
		<< ", runNo = " << runNo
		<< " (" << runNoforsigMC << ")"
		<< ", evtNo = " << evtNo
		<< ", mc = " << mc
		<< " ******************************"
		<< std::endl;

    if( fl_hadtag_dump ){
      FullRecon( fl_message_hadtag, fl_hadtag_dump );
      return 0;
    }
    
    if( flag_hadtag && FullReconN( B0_LUND, 0.005, 0)==0 ) return 0;

    Gen_info( child_id_map_D,   n_particle_set_D,
	      child_id_map_tau, n_particle_set_tau,
	      gen_b_decay_info,
	      gen_d_mode_info, gen_tau_mode_info,
	      fl_message_gen, fl_gen_dump );

    if( fl_message_gen ){
      std::cout << "*************** gen_b_decay_info[][] *****************" << std::endl;
      for(int i=0; i<40; i++ ){
	std::cout << std::setw(3) << std::right << i << " : "
		  << std::setw(7) << std::right << gen_b_decay_info[0][i] << ", "
		  << std::setw(7) << std::right << gen_b_decay_info[1][i] << std::endl;
      }
      std::cout << "******************************************************" << std::endl;
    }
    
    Tracks_cand( trk_list,                fl_message_funda );
    EE_cand    ( lep_list,    0.1, -10.0, fl_message_funda );
    Mu_cand    ( lep_list,    0.1, -10.0, fl_message_funda );
    Gamma_cand ( gamma_list,              fl_message_funda );
    if( !flag_single && !flag_hadtag && lep_list.size() < 2 ) return 0;

    Pion_cand ( pion_list,      fl_message_funda );
    Kaon_cand ( kaon_list, 0.1, fl_message_funda );
    Ks_cand   ( ks_list,     0, fl_message_funda );

    Pi0_cand  ( pi0_list,          masscut_pi0_L, masscut_pi0_H, 0.0, 0.050, 0.050, -10.0, fl_message_funda ); // for D  decay ( cos_cut, gam_e_cut_L, gam_e_cut_H, pi0_e_cut ), 50 MeV, cos>0.0
    Pi0_cand  ( pi0_for_dstr_list, masscut_pi0_L, masscut_pi0_H, -10, 0.020, 0.050, -10.0, fl_message_funda ); // for D* decay, 20 MeV and 50 MeV
    // pion study                                         <massL><massH> <cos> <EgamL><EgamH> <E_pi0>
    std::vector<Particle> pi0_test1; Pi0_cand( pi0_test1, -0.010, 0.010,   0.0, 0.050, 0.050, -10.0, fl_message_funda ); // default setting is added @ 20141127
    std::vector<Particle> pi0_test2; Pi0_cand( pi0_test2, -0.010, 0.010, -10.0, 0.050, 0.050, -10.0, fl_message_funda ); // cos theta
    std::vector<Particle> pi0_test3; Pi0_cand( pi0_test3, -0.015, 0.015,   0.0, 0.050, 0.050, -10.0, fl_message_funda ); // mass range
    std::vector<Particle> pi0_test4; Pi0_cand( pi0_test4, -0.030, 0.030,   0.0, 0.050, 0.050, -10.0, fl_message_funda ); // mass range
    std::vector<Particle> pi0_test5; Pi0_cand( pi0_test5, -0.010, 0.010,   0.0, 0.020, 0.050, -10.0, fl_message_funda ); // energy
    std::vector<Particle> pi0_test6; Pi0_cand( pi0_test6, -0.015, 0.015, -10.0, 0.050, 0.050, -10.0, fl_message_funda ); // cos theta and mass range
    std::vector<Particle> pi0_test7; Pi0_cand( pi0_test7, -0.025, 0.015, -10.0, 0.050, 0.050, -10.0, fl_message_funda ); // cos theta and mass range
    std::vector<Particle> pi0_test8; Pi0_cand( pi0_test8, -0.015, 0.015, -10.0, 0.020, 0.050, -10.0, fl_message_funda ); // cos theta and mass range and energy
    //

    // D0
    Rec_D_2body( recd0_list,    kaon_list, pion_list,  fl_message_rec ); // K+ pi-
    Rec_D_2body( recd0_list,    ks_list,   pi0_list,   fl_message_rec ); // Ks pi0
    Rec_D_2body( recd0_list,    kaon_list,  fl_message_rec ); // K+ K-
    Rec_D_2body( recd0_list,    pion_list,  fl_message_rec ); // pi+ pi-

    Rec_D_3body( recd0_list,   kaon_list, pion_list, pi0_list, fl_message_rec ); // K+ pi- pi0
    Rec_D_3body( recd0_list,     ks_list, pion_list,           fl_message_rec ); // Ks pi+ pi-
    Rec_D_3body( recd0_list,    pi0_list, pion_list,           fl_message_rec ); // pi+ pi- pi0
    Rec_D_3body( recd0_list,     ks_list, kaon_list,           fl_message_rec ); // Ks K+ K-

    Rec_D_4body( recd0_list,     ks_list, pion_list, pi0_list, fl_message_rec ); // Ks pi+ pi- pi0
    Rec_D_4body( recd0_list,   kaon_list, pion_list,           fl_message_rec ); // K+ pi- pi+ pi-
    // D+
    //Rec_D_2body( recdp_list,    kaon_list, pi0_list,   fl_message_rec ); // not allowed
    Rec_D_2body( recdp_list,     ks_list, pion_list,           fl_message_rec ); // Ks pi+
    Rec_D_3body( recdp_list,     ks_list, pion_list, pi0_list, fl_message_rec ); // Ks pi+ pi0
    Rec_D_3body( recdp_list,   kaon_list, pion_list,           fl_message_rec ); // K+ pi+ pi-
    Rec_D_3body( recdp_list,   pion_list, kaon_list,           fl_message_rec ); // K+ K- pi+
    //Rec_D_4body( recdp_list,   kaon_list, pion_list, pi0_list, fl_message_rec ); // K+ pi+ pi- pi0 // removed due to too bas S/N
    //Rec_D_4body( recdp_list,   pion_list, kaon_list, pi0_list, fl_message_rec ); // K+ K-  pi+ pi0 // removed due to too bas S/N
    Rec_D_4body( recdp_list,     ks_list, pion_list,           fl_message_rec ); // Ks pi+ pi- pi+

    // D*+
    //Rec_Dstr ( recdstrp_list, recd0_list, pion_list,         0.005,  fl_message_rec ); // D0 pi+
    //Rec_Dstr ( recdstrp_list, recdp_list, pi0_for_dstr_list, 0.005,  fl_message_rec ); // D+ pi0
    Rec_Dstr ( recdstrp_list, recd0_list, pion_list,         0.010,  fl_message_rec ); // D0 pi+ // expand to 10 MeV @20150130
    Rec_Dstr ( recdstrp_list, recdp_list, pi0_for_dstr_list, 0.010,  fl_message_rec ); // D+ pi0 // expand to 10 MeV @20150130


    // D*0
    //Rec_Dstr ( recdstr0_list, recd0_list, pi0_for_dstr_list, 0.005,  fl_message_rec ); // D0 pi0
    //Rec_Dstr ( recdstr0_list, recd0_list, gamma_list,        0.005,  fl_message_rec ); // D0 gamma

    //Rec_Dtau( recdtau_list, recd0_list,    lep_list, gamma_list, fl_message_rec ); // D0
    //Rec_Dtau( recdtau_list, recdstr0_list, lep_list, gamma_list, fl_message_rec ); // D*0
    //Rec_Dtau( recdtau_list, recdp_list,    lep_list, gamma_list, fl_message_rec ); // D+
    Rec_Dtau( recdtau_list, recdstrp_list, lep_list, gamma_list, fl_message_rec ); // D*+

    if( fl_message_ana )
      std::cout << " # of D0     candidates : " << recd0_list   .size()    << std::endl
		<< " # of D+     candidates : " << recdp_list   .size()    << std::endl
		<< " # of D*0    candidates : " << recdstr0_list.size()    << std::endl
		<< " # of D*+    candidates : " << recdstrp_list.size()    << std::endl
		<< " # of D(*)l  candidates : " <<  recdtau_list.size()    << std::endl;
    
    if( flag_single ){
      int cnt_cand = Ana_single( recdtau_list, gen_b_decay_info, gen_d_mode_info, gen_tau_mode_info, gamma_list, lep_list, Rec_single_dist, true, fl_message_ana );
      if( cnt_cand ) Skim_ana();
    }else if( flag_hadtag ){
      int cnt_cand = Ana_hadtag( recdtau_list, gen_b_decay_info, gen_d_mode_info, gen_tau_mode_info,
				 trk_list, gamma_list, lep_list, pi0_list, ks_list,
				 pi0_for_dstr_list, pi0_test1, pi0_test2, pi0_test3, pi0_test4, pi0_test5, pi0_test6, pi0_test7, pi0_test8,
				 Rec_hadtag_dist, fl_message_ana ); // which pi0 is better ? (pi0_list / pi0_dstr_list)
      if( cnt_cand ) Skim_ana();
    }else{
      int cnt_cand = Ana( recdtau_list, gen_b_decay_info, gen_d_mode_info, gen_tau_mode_info,
			  trk_list, gamma_list, pi0_list, ks_list,
			  pi0_for_dstr_list, pi0_test1, pi0_test2, pi0_test3, pi0_test4, pi0_test5, pi0_test6, pi0_test7, pi0_test8,
			  Rec_dist, fl_message_ana ); // which pi0 is better ? (pi0_list / pi0_dstr_list)
      if( cnt_cand ) Skim_ana();
    }
    return 0;
  }

  void DSTRTAUNU::Skim_ana( ){
    if( !flag_SkimFile ) return;
    Skim_dist->column("exp",   expNo);
    Skim_dist->column("run",   runNo);
    Skim_dist->column("evt",   evtNo);
    Skim_dist->column("skim",  numberOfSkim++);
    Skim_dist->column("event", numberOfEvent);
    Skim_dist->dumpData();
    SkimFile->write();
    
    return;
  }


  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
