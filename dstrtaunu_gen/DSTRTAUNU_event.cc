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
    
    std::vector<Particle> recd_list;
    std::vector<Particle> recdstr_list;
    std::vector<Particle> recdtau_list;
    
    std::map<int,int>  child_id_map_D    [2][2]; // < hep.get_ID(), hep.idhep() >
    std::multiset<int> n_particle_set_D  [2][2]; // < LUND >
    std::map<int,int>  child_id_map_tau  [2][2]; // < hep.get_ID(), hep.idhep() >
    std::multiset<int> n_particle_set_tau[2][2]; // < LUND >

    int gen_b_decay_info[2][30]={{0},{0}};
    // [0][] is signal-B, [1][] is normalization-B, if there is signal events.
    // [][ 0] B->Get_ID()
    // [][ 1] B semi-decay : 0(no lepton), 1(e-nu), 2(mu-nu), 3(tau-nu), -1(other two leptons), -2(others)
    // [][ 2] lepton->Get_ID(), lepton from semileptonic-B decay
    // [][ 3] # of B's (direct) children
    // [][ 4] # of gamma directly from B
    // [][ 5] # of rootD 'directly' from B
    // [][ 6] 1st-rootD->idhep()
    // [][ 7] 1st-rootD->get_ID()
    // [][ 8] 2nd-rootD->idhep()
    // [][ 9] 2nd-rootD->get_ID()
    // [][10] # of D in B decay, where D is [D0/D+/Ds]
    // [][11] 1st-D->Get_ID()
    // [][12] 2nd-D->Get_ID()
    // [][13] flag of D/D*/D**/DD
    // [][14] accompany particle with D* decay->idhep()
    // [][15] accompany particle with D* decay->get_ID()
    
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
    const int  fl_message_gen   = 0; // 0(no-message),1(gen),2(gen+table)
    const bool fl_message_funda = !true;
    const bool fl_message_rec   = !true;
    const bool fl_message_ana   = !true;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( fl_message_gen || fl_message_funda || fl_message_rec || fl_message_ana )
      std::cout << "*******************************"
		<< " [ DSTRTAUNU ] "
		<< "expNo = "   << expNo
		<< ", runNo = " << runNo
		<< ", evtNo = " << evtNo
		<< " ******************************"
		<< std::endl;

    Gen_info( child_id_map_D,   n_particle_set_D,
	      child_id_map_tau, n_particle_set_tau,
	      gen_b_decay_info,
	      gen_d_mode_info, gen_tau_mode_info,
	      fl_message_gen );
    /*
    if( fl_message_gen ){
      std::cout << "*************** gen_b_decay_info[][] *****************" << std::endl;
      for(int i=0; i<30; i++ ){
	std::cout << std::setw(3) << std::right << i << " : "
		  << std::setw(7) << std::right << gen_b_decay_info[0][i] << ", "
		  << std::setw(7) << std::right << gen_b_decay_info[1][i] << std::endl;
      }
      std::cout << "******************************************************" << std::endl;
    }

    Tracks_cand( trk_list,    fl_message_funda );
    EE_cand    ( lep_list,    fl_message_funda );
    Mu_cand    ( lep_list,    fl_message_funda );
    Gamma_cand ( gamma_list,  fl_message_funda );
    
    Pion_cand ( pion_list,    fl_message_funda );
    Kaon_cand ( kaon_list,    fl_message_funda );
    Pi0_cand  ( pi0_list,     fl_message_funda );
    Ks_cand   ( ks_list,   0, fl_message_funda );

    Rec_D_2body( recd_list,    kaon_list, pion_list,  fl_message_rec );
    //Rec_D_2body( recd_list,    kaon_list, pi0_list,   fl_message_rec ); // not allowed
    Rec_D_2body( recd_list,    ks_list,   pion_list,  fl_message_rec );
    //Rec_D_2body( recd_list,    ks_list,   pi0_list,   fl_message_rec );

    //Rec_D_2body( recd_list,    kaon_list,  fl_message_rec );
    //Rec_D_2body( recd_list,    pion_list,  fl_message_rec );

    //Rec_D_3body( recd_list,   kaon_list, pion_list, pi0_list, fl_message_rec );
    //Rec_D_3body( recd_list,     ks_list, pion_list, pi0_list, fl_message_rec );
    //Rec_D_3body( recd_list,     ks_list, pion_list,           fl_message_rec );
    //Rec_D_3body( recd_list,    pi0_list, pion_list,           fl_message_rec );
    //Rec_D_3body( recd_list,     ks_list, kaon_list,           fl_message_rec );
    //Rec_D_3body( recd_list,   kaon_list, pion_list,           fl_message_rec );
    //Rec_D_3body( recd_list,   pion_list, kaon_list,           fl_message_rec );

    //Rec_D_4body( recd_list,   kaon_list, pion_list, pi0_list, fl_message_rec );
    //Rec_D_4body( recd_list,     ks_list, pion_list, pi0_list, fl_message_rec );
    //Rec_D_4body( recd_list,   kaon_list, pion_list,           fl_message_rec );
    //Rec_D_4body( recd_list,     ks_list, pion_list,           fl_message_rec );

    Rec_Dstr ( recdstr_list, recd_list, pion_list,  fl_message_rec );
    //Rec_Dstr ( recdstr_list, recd_list, pi0_list,   fl_message_rec );
    //Rec_Dstr ( recdstr_list, recd_list, gamma_list, fl_message_rec );

    Rec_Dtau( recdtau_list, recd_list,    lep_list, gamma_list, fl_message_rec );
    Rec_Dtau( recdtau_list, recdstr_list, lep_list, gamma_list, fl_message_rec );

    if( fl_message_ana ) std::cout << " # of D     candidates : " << recd_list   .size()    << std::endl
				   << " # of D*    candidates : " << recdstr_list.size()    << std::endl
				   << " # of D(*)l candidates : " << recdtau_list.size()    << std::endl;
    
    Ana ( recdtau_list, gen_b_decay_info, gen_d_mode_info, gen_tau_mode_info,
	  trk_list, gamma_list, pi0_list, ks_list,
	  fl_message_ana );
    */
  }

  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
