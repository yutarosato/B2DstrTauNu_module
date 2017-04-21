#include "belle.h"
#include "DSTRTAUNU.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif


  void DSTRTAUNU::NP_Weight( std::map<int, int> child_id_map_d  [][2], std::multiset<int> n_particle_set_d  [][2],
			     std::map<int, int> child_id_map_tau[][2], std::multiset<int> n_particle_set_tau[][2],
			     int gen_b_decay_info[][40],
			     int gen_d_mode_info [][2][6], int gen_tau_mode_info[][6], 
			     const int fl_message, const bool fl_dump )
  {

    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();

    if( gen_b_decay_info[0][1]!=3 ) return; // not semitauonic decay

    double exprunNo = expNo*10000+runNo;
    if( runNo==0 ) exprunNo = expNo*10000+runNoforsigMC; // for signal MC
    
    NP_dist->column( "exprun",    exprunNo   );
    NP_dist->column( "evt_type",  event_type );

    // LUND
    NP_dist->column( "blund",   Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[0][0]) ) ).idhep() );
    NP_dist->column( "dlund",   gen_b_decay_info[0][6] );  // root-D 
    NP_dist->column( "acclund", gen_b_decay_info[0][14] ); // accompany particle in D* decay
    NP_dist->column( "d2lund",  gen_b_decay_info[0][27] ); // D from D* decay
    // Decay Mode
    int dec_d   = -10;
    if     ( abs(gen_b_decay_info[0][6])==421                                      ) dec_d = 0; // D0
    else if( abs(gen_b_decay_info[0][6])==411                                      ) dec_d = 1; // D+
    else if( abs(gen_b_decay_info[0][6])==423 &&     gen_b_decay_info[0][14] ==111 ) dec_d = 2; // D*0 -> D0 pi0
    else if( abs(gen_b_decay_info[0][6])==423 &&     gen_b_decay_info[0][14] == 22 ) dec_d = 3; // D*0 -> D0 gamma
    else if( abs(gen_b_decay_info[0][6])==413 && abs(gen_b_decay_info[0][14])==211 ) dec_d = 4; // D*+ -> D0 pi+
    else if( abs(gen_b_decay_info[0][6])==413 &&     gen_b_decay_info[0][14] ==111 ) dec_d = 5; // D*+ -> D+ pi0
    else if( abs(gen_b_decay_info[0][6])==413 &&     gen_b_decay_info[0][14] == 22 ) dec_d = 6; // D*+ -> D+ gamma
    NP_dist->column( "dec_d",   dec_d );
    
    int dec_tau = -10;

    if( gen_b_decay_info[0][1]==3 ){
      Gen_hepevt& Tau( genMgr( Panther_ID(gen_b_decay_info[0][2]) ) );
      int taunchild = Tau.daLast()==0 ? 0 : Tau.daLast() - Tau.daFirst() + 1;
      int cnt_dec_tau = 0;
      for( int i=0; i<taunchild; i++){
	if( genMgr( Panther_ID(Tau.daFirst()+i) ).idhep() == Gamma_LUND ) continue;
	cnt_dec_tau += abs( genMgr( Panther_ID(Tau.daFirst()+i) ).idhep() );
      }
      if     ( cnt_dec_tau == 11 + 12 + 16 ) dec_tau = 1; // e  nu nu
      else if( cnt_dec_tau == 13 + 14 + 16 ) dec_tau = 2; // mu nu nu
      else if( cnt_dec_tau ==     211 + 16 ) dec_tau = 3; // pi  nu
      else if( cnt_dec_tau ==     213 + 16 ) dec_tau = 4; // rho nu
      else if( cnt_dec_tau ==   20213 + 16 ) dec_tau = 5; // a_1 nu
      
      else if( cnt_dec_tau ==     321 + 16 ) dec_tau = 8; // K nu
      else if( cnt_dec_tau ==     323 + 16 ) dec_tau = 9; // K* nu

      else if( cnt_dec_tau == 16                         +   211 +   111 ) dec_tau = 11; //     pi  pi0
      else if( cnt_dec_tau == 16                         +   211 + 2*111 ) dec_tau = 12; //     pi 2pi0
      else if( cnt_dec_tau == 16                         +   211 + 3*111 ) dec_tau = 13; //     pi 3pi0
      else if( cnt_dec_tau == 16                         +   211 + 4*111 ) dec_tau = 14; //     pi 4pi0
      else if( cnt_dec_tau == 16                         + 3*211 +   111 ) dec_tau = 15; //    3pi  pi0
      else if( cnt_dec_tau == 16                         + 3*211 + 2*111 ) dec_tau = 16; //    3pi 2pi0
      else if( cnt_dec_tau == 16                         + 3*211 + 3*111 ) dec_tau = 17; //    3pi 3pi0
      else if( cnt_dec_tau == 16                         + 5*211         ) dec_tau = 18; //    5pi
      else if( cnt_dec_tau == 16                         + 5*211 +   111 ) dec_tau = 19; //    5pi  pi0
      else if( cnt_dec_tau == 16         +   321                 +   111 ) dec_tau = 21; //  K         pi0
      else if( cnt_dec_tau == 16         +   321                 + 2*111 ) dec_tau = 22; //  K        2pi0
      else if( cnt_dec_tau == 16         +   321                 + 3*111 ) dec_tau = 23; //  K        3pi0
      else if( cnt_dec_tau == 16         +   321         + 2*211         ) dec_tau = 24; //  K    2pi
      else if( cnt_dec_tau == 16         +   321         + 2*211 +   111 ) dec_tau = 25; //  K    2pi  pi0
      else if( cnt_dec_tau == 16         + 2*321         +   211         ) dec_tau = 26; // 2K     pi
      else if( cnt_dec_tau == 16         + 2*321         +   211 +   111 ) dec_tau = 27; // 2K     pi  pi0
      else if( cnt_dec_tau == 16         +   321 +   311                 ) dec_tau = 31; //  K K0 
      else if( cnt_dec_tau == 16                 +   311 +   211         ) dec_tau = 32; //    K0  pi
      else if( cnt_dec_tau == 16                 +   311 + 3*211         ) dec_tau = 33; //    K0 3pi
      else if( cnt_dec_tau == 16                 +   311 +   211 +   111 ) dec_tau = 34; //    K0  pi pi0
      else if( cnt_dec_tau == 16             + 310 + 130 +   211         ) dec_tau = 35; //    K0 K0 pi     (KS KL pi    )
      else if( cnt_dec_tau == 16             + 310 + 130 +   211 +   111 ) dec_tau = 36; //    K0 K0 pi pi0 (KS KL pi pi0)
      else if( cnt_dec_tau == 16             + 310 + 310 +   211         ) dec_tau = 37; //    K0 K0 pi     (KS KS pi    )
      else if( cnt_dec_tau == 16             + 310 + 310 +   211 +   111 ) dec_tau = 38; //    K0 K0 pi     (KS KS pi pi0)
      else if( cnt_dec_tau == 16 +   223                 +   211         ) dec_tau = 41; //  w      pi
      else if( cnt_dec_tau == 16 +   223                 +   211 +   111 ) dec_tau = 42; //  w      pi  pi0
      else if( cnt_dec_tau == 16 +   223                 +   211 + 2*111 ) dec_tau = 43; //  w      pi 2pi0
      else if( cnt_dec_tau == 16 +   223                 + 3*211         ) dec_tau = 44; //  w     3pi
      else if( cnt_dec_tau == 16 +   223                 + 3*211 +   111 ) dec_tau = 45; //  w     3pi  pi0
      else if( cnt_dec_tau == 16 +   223                 + 3*211 + 2*111 ) dec_tau = 46; //  w     3pi 2pi0
      else if( cnt_dec_tau == 16 +   221 +  321                          ) dec_tau = 51; //  K eta
      else                                                                 dec_tau =  0; // other

    }
    NP_dist->column( "dec_tau",   dec_tau );

    // Kinematics
    double q2        = -10;
    double cos_d     = -10;
    double cos_w     = -10;
    double cos_t     = -10;
    double cos_dau   = -10;
    double p_dau_tau = -10;
    double p_dau_w   = -10;
    double p_dau_b   = -10;
    double p_dau_lab = -10;
    double p_dau_cm  = -10; // added @ 20160221

    double theta_dau_tau = -10;
    double theta_dau_w   = -10;
    double theta_dau_b   = -10;
    double theta_dau_lab = -10;

    calKinematics_np( gen_b_decay_info[0][0], gen_b_decay_info[0][7], gen_b_decay_info[0][2],
		      dec_d, dec_tau,
		      q2, cos_d, cos_w, cos_t,
		      cos_dau,
		      p_dau_tau, p_dau_w, p_dau_b, p_dau_lab, p_dau_cm );

    NP_dist->column( "q2",        q2        );
    NP_dist->column( "cos_d",     cos_d     );
    NP_dist->column( "cos_w",     cos_w     );
    NP_dist->column( "cos_t",     cos_t     );
    NP_dist->column( "cos_dau",   cos_dau   );
    NP_dist->column( "p_dau_tau", p_dau_tau );
    NP_dist->column( "p_dau_w",   p_dau_w   );
    NP_dist->column( "p_dau_b",   p_dau_b   );
    NP_dist->column( "p_dau_lab", p_dau_lab );
    NP_dist->column( "p_dau_cm",  p_dau_cm  ); // added @20160221

    // momentum
    Gen_hepevt       Dstr1 = Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[0][7]) ) );
    HepLorentzVector Dstr1_4V  ( Dstr1.PX(), Dstr1.PY(), Dstr1.PZ(), Dstr1.E() );
    HepLorentzVector Dstr1_4Vcm( Dstr1.PX(), Dstr1.PY(), Dstr1.PZ(), Dstr1.E() );
    Dstr1_4Vcm.boost( cmboost );
    // dstp & dstpcm
    NP_dist->column( "dstp",   Dstr1_4V.vect().mag()   );
    NP_dist->column( "dstpcm", Dstr1_4Vcm.vect().mag() );

    // fitting variables
    if( gen_b_decay_info[0][ 1]!=3 || // not semitauonic decay
	gen_b_decay_info[0][22]==0 || // not leptonic tau decay
	!(gen_b_decay_info[1][1]==1 || gen_b_decay_info[1][1]==2) || // not semileptonic decay
	gen_b_decay_info[1][2]==0
	){
      NP_dist->column( "evis2",  -10 );
      NP_dist->column( "mmiss2", -10 );
      NP_dist->column( "cosdl1",  10 );
      NP_dist->column( "cosdl2",  10 );
      NP_dist->column( "cosdll",  10 );
      NP_dist->column( "cosdlh",  10 );
      if( fl_dump ) NP_dist->dumpData();
      return;
    }

    Gen_hepevt       Dstr2 = Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[1][7]) ) );
    HepLorentzVector Dstr2_4V  ( Dstr2.PX(), Dstr2.PY(), Dstr2.PZ(), Dstr2.E() );
    HepLorentzVector Dstr2_4Vcm( Dstr2.PX(), Dstr2.PY(), Dstr2.PZ(), Dstr2.E() );
    Dstr2_4Vcm.boost( cmboost );

    Gen_hepevt       lep1 = Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[0][22]) ) );
    HepLorentzVector lep1_4V  ( lep1.PX(), lep1.PY(), lep1.PZ(), lep1.E() );
    HepLorentzVector lep1_4Vcm( lep1.PX(), lep1.PY(), lep1.PZ(), lep1.E() );
    lep1_4Vcm.boost( cmboost );


    Gen_hepevt       lep2 = Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[1][ 2]) ) );
    HepLorentzVector lep2_4V  ( lep2.PX(), lep2.PY(), lep2.PZ(), lep2.E() );
    HepLorentzVector lep2_4Vcm( lep2.PX(), lep2.PY(), lep2.PZ(), lep2.E() );
    lep2_4Vcm.boost( cmboost );
    
    double     mmiss2 = 0;
    double     Evis   = 0;
    Hep3Vector Ptot(0);
    Evis += Dstr1_4Vcm.e();
    Evis += Dstr2_4Vcm.e();
    Evis +=  lep1_4Vcm.e();
    Evis +=  lep2_4Vcm.e();

    Ptot += Dstr1_4Vcm.vect();
    Ptot += Dstr2_4Vcm.vect();
    Ptot +=  lep1_4Vcm.vect();
    Ptot +=  lep2_4Vcm.vect();
    mmiss2 = (2*eb-Evis)*(2*eb-Evis) - Ptot.mag2();

    NP_dist->column( "evis2",  Evis                                       );
    NP_dist->column( "mmiss2", mmiss2 > 0 ? sqrt(mmiss2) : -sqrt(-mmiss2) );

    Gen_hepevt B1 = Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[0][0]) ) );
    Gen_hepevt B2 = Gen_hepevt( genMgr( Panther_ID(gen_b_decay_info[1][0]) ) );
    double Bmass1 = B1.M();
    double Bmass2 = B2.M();

    HepLorentzVector Dtau1_4V   =  Dstr1_4V   + lep1_4V;
    HepLorentzVector Dtau1_4Vcm =  Dstr1_4Vcm + lep1_4Vcm;
    HepLorentzVector Dtau2_4V   =  Dstr2_4V   + lep2_4V;
    HepLorentzVector Dtau2_4Vcm =  Dstr2_4Vcm + lep2_4Vcm;
    double cosBDl1 = ( 2.0*eb*Dtau1_4Vcm.e() - Bmass1*Bmass1-Dtau1_4V.m()*Dtau1_4V.m() ) / ( 2.0*sqrt(eb*eb-Bmass1*Bmass1)*Dtau1_4Vcm.vect().mag() );
    double cosBDl2 = ( 2.0*eb*Dtau2_4Vcm.e() - Bmass2*Bmass2-Dtau2_4V.m()*Dtau2_4V.m() ) / ( 2.0*sqrt(eb*eb-Bmass2*Bmass2)*Dtau2_4Vcm.vect().mag() );
    NP_dist->column( "cosdl1", cosBDl1 );
    NP_dist->column( "cosdl2", cosBDl2 );
    NP_dist->column( "cosdll", (cosBDl1 < cosBDl2 ? cosBDl1 : cosBDl2) );
    NP_dist->column( "cosdlh", (cosBDl1 < cosBDl2 ? cosBDl2 : cosBDl1) );

    if( fl_dump ) NP_dist->dumpData();
  }
  
#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
