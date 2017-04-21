#include "belle.h"
#include "DSTRTAUNU.h"
//#include "toolbox/FoxWolfr.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif
  
  int DSTRTAUNU::check_selfF(Particle& particle){ // final state particle(K,PI,E,MU,Proton)
    if( mc==0 ) return -1;
    int self=0;
    if( particle.mdstCharged() ){
      //const Gen_hepevt & hep( get_hepevt(particle.mdstCharged()) );
      const Gen_hepevt & hep( gen_level(get_hepevt(particle.mdstCharged())) ); // test @20150123
      particle.relation().genHepevt( hep );
      if ( hep && particle.pType().lund() == hep.idhep() ) self=1;
    }
    return self;
  }

  int DSTRTAUNU::check_selfF(Particle& particle, int self_LUND, int mother_LUND, int mothermother_LUND){ // particle originated from specific decay mode?
    if( mc==0 ) return -1;
    if( particle.mdstCharged() ){
      const Gen_hepevt & hep( get_hepevt(particle.mdstCharged()) );
      if ( !(hep && abs(particle.pType().lund())==self_LUND && particle.pType().lund()==hep.idhep()) ){ // self
	return 0;
      }else{
	if( mother_LUND==0 ){
	  return 1;
	}else{
	  if( !(hep.mother() && hep.mother().idhep()==mother_LUND) ){ // mother
	    return 0;
	  }else{
	    if( mothermother_LUND==0 ){
	      return 1;
	    }else{
	      if( !(hep.mother().mother() && hep.mother().mother().idhep() == mothermother_LUND) ){ // mother-mohter
		return 0;
	      }else{
		return 1;
	      }
	    }
	  }
	}
      }
    }
    return 0;
  }

  int DSTRTAUNU::check_idF(Particle& particle){ // final state particle(K,PI,E,MU,Proton)
    if( mc==0 ) return -1;
    int id=0;
    if( particle.mdstCharged() ){
      const Gen_hepevt & hep( get_hepevt(particle.mdstCharged()) );
      if( hep ) id = hep.idhep();
    }
    return id;
  }

  int DSTRTAUNU::check_motheridF(Particle& particle){ // final state particle(K,PI,E,MU,Proton)
    if( mc==0 ) return -1;
    int motherid=0;
    if( particle.mdstCharged() ){
      const Gen_hepevt & hep( get_hepevt(particle.mdstCharged()) );
      if( hep && hep.mother() ) motherid = hep.mother().idhep();
    }
    return motherid;
  }

  int DSTRTAUNU::check_selfG(Particle& particle){ // gamma
    if( mc==0 ) return -1;
    int self=0;
    if( particle.mdstGamma() ){
      const Gen_hepevt & hep( gen_level(get_hepevt(particle.mdstGamma())) );
      if( hep && particle.pType().lund() == hep.idhep() ){
	particle.relation().genHepevt( hep );
	self=1;
      }
    }
    return self;
  }
    
  int DSTRTAUNU::check_selfR( Particle& particle, int specified_nchild ){ // reconstructed particles
    if( mc==0 ) return -1;

    int nchild = particle.nChildren();
    if( specified_nchild ) nchild = specified_nchild;
    // check whether genHep references exist;
    for(int i=0; i<nchild; ++i){
      if( !particle.child(i).genHepevt() ) return 0;
      if( abs(particle.child(i).pType().lund()) != abs(particle.child(i).genHepevt().idhep()) ) return 0; // added 20131204
    }
    
    
    // check duplication and same mother(we do not distinguish B0-antiB0 and D0-antiD0)
    for(int i=0; i<nchild-1; ++i){
      for(int j=i+1; j<nchild; ++j){
	if( particle.child(i).genHepevt().get_ID() == particle.child(j).genHepevt().get_ID() ) return 0;
	if( particle.child(i).genHepevt().mother().get_ID() != particle.child(j).genHepevt().mother().get_ID()  ) return 0;
	if(  abs(particle.pType().lund()) == B0_LUND || abs(particle.pType().lund()) == D0_LUND ){
	  if( abs(particle.pType().lund()) != abs(particle.child(i).genHepevt().mother().idhep()) ) return 0;
	}else{
	  if( particle.pType().lund() != particle.child(i).genHepevt().mother().idhep() ) return 0;
	}
      }
    }

    // check the number of children
    int truenchild = particle.child(0).genHepevt().mother().daLast()==0 ? 0 : particle.child(0).genHepevt().mother().daLast() - particle.child(0).genHepevt().mother().daFirst() + 1;
    if( nchild != truenchild ){ // allow the non-reconstructed low-E gamma 
      Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
      Gen_hepevt& extra_child( genMgr(Panther_ID(particle.child(0).genHepevt().mother().daFirst() + nchild )) );
      if( !(extra_child.idhep() == 22 && extra_child.E() < 0.050) ) return 0;
    }
    
    particle.relation().genHepevt( particle.child(0).genHepevt().mother() );
    return 1;
  }

  int DSTRTAUNU::check_selfPi0( Particle& particle ){ // for pi0
    if( mc==0 ) return -1;

    Particle& gam0 = particle.child(0);
    Particle& gam1 = particle.child(1);

    if( gam0.genHepevt()          && gam1.genHepevt() &&
	gam0.genHepevt().mother() && gam1.genHepevt().mother() &&
	gam0.genHepevt().mother().idhep()  == PI0_LUND && 
	gam1.genHepevt().mother().idhep()  == PI0_LUND &&
	gam0.genHepevt().mother().get_ID() == gam1.genHepevt().mother().get_ID() 
	){ // perfect pi0
      particle.relation().genHepevt( particle.child(0).genHepevt().mother() );
      return 1;
    }else if( gam0.e() > gam1.e() &&
	      gam0.genHepevt() && gam0.genHepevt().mother() &
	      gam0.genHepevt().mother().idhep()==PI0_LUND
	      ){ // pseudo pi0
      particle.relation().genHepevt( gam0.genHepevt().mother() );      
      return 2;
    }else if( gam0.e() < gam1.e() &&
	      gam1.genHepevt() && gam1.genHepevt().mother() &
	      gam1.genHepevt().mother().idhep()==PI0_LUND
	      ){ // pseudo pi0
      particle.relation().genHepevt( gam1.genHepevt().mother() );      
      return 2;
    }
    return 0;
  }

  int DSTRTAUNU::check_selfR2( Particle& particle ){ // reconstructed particles, considering inter-mediate state
    if( mc==0 ) return -1;

    int nchild = particle.nChildren();

    // check that genHep references exist;
    for(int i=0; i<nchild; ++i){
      if( !particle.child(i).genHepevt() ) return 0;
      if( abs(particle.child(i).pType().lund()) != abs(particle.child(i).genHepevt().idhep()) ) return 0; // added 20131204
    }
    
    // check duplication
    for(int i=0; i<nchild-1; ++i){
      for(int j=i+1; j<nchild; ++j){
	if( particle.child(i).genHepevt().get_ID() == particle.child(j).genHepevt().get_ID() ) return 0;
      }
    }
    
    // check if daughter have mother with correct lund-code
    int id1 = 0;
    int id2 = 0;
    double energy = 0;
    for(int i=0; i<nchild; ++i){
      if( i==0 ){
	id1 = search_origin(particle.child(i).genHepevt(), particle.pType().lund());
	id2 = id1;
      }else id2 = search_origin(particle.child(i).genHepevt(), particle.pType().lund());
      energy += particle.child(i).genHepevt().E();
      if( !id1 || !id2 || id1!=id2 ) return 0;
    }

    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    Gen_hepevt& gen( genMgr(Panther_ID(id1)) );
    if( abs(gen.E()-energy)<0.05 ){
      particle.relation().genHepevt( gen );
      return 1;
    }else{
      return 0;
    }
  }

  void DSTRTAUNU::display_particle(Particle particle){
    std::cout << " ******* < display particle > *************************************** " << std::endl
	      << "   < position >   : " << particle.x()
	      << particle.momentum().dx()
	      << "   < momentum >   : " << particle.p()
	      << particle.momentum().dp()
	      << " < error matrix > "   << particle.momentum().dpx()
	      << " < prod. vertex > : " << particle.momentum().vertex()
	      << particle.momentum().dVertex()
	      << " < decay vertex > : " << particle.momentum().decayVertex()
	      << particle.momentum().dDecayVertex()
	      << "LUND = " << particle.pType().lund() << "(" << particle.pType().name() << ")" << std::endl
	      << "nChildren = " << particle.nChildren() << std::endl;
    for(int i=0; i<particle.nChildren(); i++)
      std::cout << "child " << i << " : " << particle.child(i).pType().lund() << "(" << particle.child(i).pType().name() << ")" << std::endl;
    //std::cout << "nFinalStateParticles = " << particle.relation().nFinalStateParticles() << std::endl;
    return;
  }
  
  void DSTRTAUNU::display_hepevt(Gen_hepevt gen){
    if( mc==0 ) return;
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast()==0 ? 0 : gen.daLast() - gen.daFirst() + 1;
    std::cout << " ******* < display hepevt > *************************************** " << std::endl
	      << " <    self ID    > : " << gen.idhep()
	      << " ( " << gen.get_ID() << " ) " << std::endl
	      << " < # of children > : " << nchild       << std::endl;
    for( int i=0; i<nchild; i++)
      std::cout << "    " << i << " : "
		<< std::setw(6) << genMgr( Panther_ID(gen.daFirst()+i) ).idhep()
		<< " ( "
		<< genMgr( Panther_ID(gen.daFirst()+i) ).get_ID()
		<< " )"
		<< std::endl;
    return;
  }

  void DSTRTAUNU::display_list( std::set<int> set, std::multiset<int> multiset ){
    std::cout << "=============================================================================" << std::endl
	      << "# of particles :  LUND code" << std::endl
	      << "-----------------------------------------------------------------------------" << std::endl;
    for (std::set<int>::iterator i = set.begin(); i != set.end(); i++){
      Ptype tmp( *i );
      std::cout << std::setw(14) << std::right << multiset.count( *i ) << " : "
		<< std::setw(10) << std::right << *i                   << " ( "
		<< std::setw(20) << std::left  << tmp.name()           << " ) "
		<< std::endl;
    }
    std::cout << "-----------------------------------------------------------------------------" << std::endl;
    return;
  }

  void DSTRTAUNU::rec_message( std::vector<Particle>::iterator p ){
    setUserInfo( *p );
    UserInfo& info_p = dynamic_cast<UserInfo&>( p->userInfo() );
    if( p->lund()==Ks_LUND ){ // Ks
      std::cout << " LUND = " << std::setw(4) << std::right << p->lund()       << " ("
		<< "cntID="   << std::setw(3) << std::right << info_p.cntid()  << ", "
		<< "self="    << std::setw(3) << std::right << info_p.self()   << ")";
    }else if( p->lund()==PI0_LUND ){ // pi0
      std::cout << " LUND = " << std::setw(4) << std::right << p->lund()      << " ("
		<< "cntID="   << std::setw(3) << std::right << info_p.cntid() << ", "
		<< "self="    << std::setw(3) << std::right << info_p.self()  << ")";
    }else if( p->lund()==Gamma_LUND ){ // gamma
      std::cout << " LUND = " << std::setw(4) << std::right << p->lund()             << " ("
		<< "gamID="   << std::setw(2) << std::right << p->mdstEcl().get_ID() << ", "
		<< "self="    << std::setw(2) << std::right << info_p.self()         << ")";
    }else if( abs(p->lund())==PIplus_LUND || abs(p->lund())==Kplus_LUND ){ // K+/pi+
      std::cout << " LUND = "      << std::setw(4) << std::right << p->lund()                 << " ("
		<< "genID="         << std::setw(3) << std::right << p->genHepevt().get_ID()   << ", "
		<< "mdstChargedID=" << std::setw(3) << std::right << p->mdstCharged().get_ID() << ", "
		<< "genLUND="       << std::setw(6) << std::right << p->genHepevt().idhep()    << ")";
    }else{ // not supported particles
      std::cerr << "[WARNING] Not supported particles for rec_message : " << p->lund() << std::endl;
    }

    return;
  }

  void DSTRTAUNU::display_rec_particle( Particle particle, const bool fl_message, int indent ){
    if( mc==0 ) return;
    if( fl_message && indent==1 ){
      std::cout << "[Reconstructed Particle]" << std::endl;
      std::cout << "==============================================================================" << std::endl
		<< " idhep() [ nchild ]    " << std::endl
		<< "------------------------------------------------------------------------------" << std::endl;
      std::cout << particle.lund()         << "("
		<< particle.pType().name() << ") ["
		<< particle.nChildren()    << "]"
		<< std::endl;
    }
    int nchild = particle.nChildren();
    for( int i=0; i<nchild; i++ ){
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << particle.child(i).lund()         << "("
		  << particle.child(i).pType().name() << ") ["
		  << particle.child(i).nChildren()    << "] ";


	if( particle.child(i).mdstCharged() ) std::cout << "(chgID)=("
							<< particle.child(i).mdstCharged().get_ID()
							<< "), ";
	if( particle.child(i).mdstEcl() ) std::cout << "(gamID)=("
							<< particle.child(i).mdstEcl().get_ID()
							<< "), ";
	if( particle.child(i).genHepevt() ){
	  std::cout << "(genHepevtID,genHepevtLUND)=("
		    << particle.child(i).genHepevt().get_ID() << ", "
		    << particle.child(i).genHepevt().idhep()  << "), ";
	}
	
	std::cout << std::endl;
      }
      if( particle.child(i).nChildren() ) display_rec_particle( particle.child(i), fl_message, indent+1 );
    } // end of child-loop
    return;
  }

  void DSTRTAUNU::find_fin_child( Gen_hepevt gen,
				  std::map<int, int>& child_id_map,
				  std::multiset<int>& n_particle_set,
				  bool fl_message, // if fl_message = true, messages are displayed.
				  int indent       // the number of indent. Users do not handle it.
				  ){
    if( mc==0 ) return;
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast()==0 ? 0 : gen.daLast() - gen.daFirst() + 1;
    if( fl_message && indent==1 ){
      std::cout << "==============================================================================" << std::endl
		<< " idhep() [ get_ID(), nchild ]    " << std::endl
		<< "------------------------------------------------------------------------------" << std::endl;
      std::cout << gen.idhep()  << "[ "
		<< gen.get_ID() << ", "
		<< nchild       << " ] "
		<< std::endl;
    }
    
    if( nchild==0 ){
      std::cout << " particles not having child" << std::endl;
      n_particle_set.insert( abs(gen.idhep()) );
    }

    for(int i=0; i<nchild; i++){
      Gen_hepevt& daughter( genMgr(Panther_ID(gen.daFirst() + i)) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
	  }
      
      if ( recon_set.count( daughter.idhep()) ){ // recon_particles
	if( fl_message ) std::cout << " --> recon_particles" << std::endl;
	if( child_id_map.count(gen.get_ID()) != 0) std::cout << "error ????? duplication" << std::endl;
	child_id_map.insert( std::pair<int, int> (daughter.get_ID(), daughter.idhep() ) );
	n_particle_set.insert( abs(daughter.idhep()) );
      }else if( nonrecon_set.count( daughter.idhep()) ){ // nonrecont_particles
	if( daughter.idhep() == Gamma_LUND && daughter.E() < 0.050 ){ // low energy gamma is neglected.
	  if( fl_message ) std::cout << "--> FSR gamma" << std::endl;
	  continue;
	}
	if( fl_message ) std::cout << "--> nonrecon_particles" << std::endl;
	n_particle_set.insert( abs(daughter.idhep()) );	
      }else{ // other particle --> continue iterations.
	if( fl_message ) std::cout << " --> other particles --> continue iterations" << std::endl;
	  find_fin_child( daughter, child_id_map, n_particle_set, fl_message, indent+1 );
	}
    }
    return;
  }

  int DSTRTAUNU::mapping_delete( Particle particle, std::map<int, int>& child_id_map,
			     int& n_recon, int gen_mode_bg,
			     char* name, const bool fl_message ){
    int fl=0;
    UserInfo& info_particle = dynamic_cast<UserInfo&>( particle.userInfo() );
    if( gen_mode_bg == 0 && particle.genHepevt() && child_id_map.count( particle.genHepevt().get_ID() ) ){
      if( child_id_map[ particle.genHepevt().get_ID() ] == particle.lund() ){
	n_recon -= 1;
	fl=1;
	if( fl_message ) std::cout << "[ ERASE : "
				   << std::setw(15) << std::right << name
				   << "] : "
				   << std::setw(10) << std::right << particle.lund() << "("
				   << std::setw(3)  << std::right << particle.genHepevt().get_ID() << ") : "
				   << n_recon     << std::endl;
      }
    }
    return fl;
  }

  int DSTRTAUNU::search_daughter_D( Gen_hepevt gen, int& id1, int& id2, bool fl_message, double threshold ){
    int indent=1;
    std::vector<int> d_id_list; // < gen(D).get_ID >
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;

    if( fl_message ){
      std::cout << "==============================================================================" << std::endl
		<< " idhep() [ get_ID(), nchild ]   --> < search_daughter_D > "                     << std::endl
		<< "------------------------------------------------------------------------------" << std::endl;
      std::cout << gen.idhep()  << "[ "
		<< gen.get_ID() << ", "
		<< nchild       << " ] "
		<< std::endl;
    }

    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
      }
      if( daughter.M() > threshold ){
	if( abs(daughter.idhep()) == D0_LUND || abs(daughter.idhep()) == Dplus_LUND || abs(daughter.idhep()) == Dsplus_LUND ){
	  d_id_list.push_back( daughter.get_ID() );
	  if( fl_message ) std::cout << " --> find D particle" << std::endl;
	}else{
	  if( fl_message ) std::cout << " --> not D particle --> continue iterations" << std::endl;
	  search_daughter_D_sub( daughter, d_id_list, fl_message, indent+1 );
	}
      }else{
	if( fl_message ) std::cout << " --> low mass particle" << std::endl;
      }
    }

    if( fl_message ) std::cout << " ===========> d_id_list.size() == " << d_id_list.size() << std::endl;
    if( d_id_list.size() > 0 ) id1 = d_id_list[0];
    if( d_id_list.size() > 1 ) id2 = d_id_list[1];
    return d_id_list.size();
  }

  int DSTRTAUNU::search_daughter_D_sub( Gen_hepevt gen, std::vector<int>& d_id_list,
					bool fl_message, int indent, double threshold ){
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;

    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
      }
      if( daughter.M() > threshold ){
	if( abs(daughter.idhep()) == D0_LUND || abs(daughter.idhep()) == Dplus_LUND || abs(daughter.idhep()) == Dsplus_LUND ){
	  d_id_list.push_back( daughter.get_ID() );
	  if( fl_message ) std::cout << " --> find D particle" << std::endl;
	}else{
	  if( fl_message ) std::cout << " --> not D particle --> continue iterations" << std::endl;
	  search_daughter_D_sub( daughter, d_id_list, fl_message, indent+1 );
	}
      }else{
	if( fl_message ) std::cout << " --> low mass particles" << std::endl;
      }
    }
  }

    int DSTRTAUNU::search_daughter_Dstr( Gen_hepevt gen, int& id1, int& id2, bool fl_message, double threshold ){
    int indent=1;
    std::vector<int> dstr_id_list; // < gen(D*).get_ID >
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;

    if( fl_message ){
      std::cout << "==============================================================================" << std::endl
		<< " idhep() [ get_ID(), nchild ]   --> < search_daughter_D* > "                    << std::endl
		<< "------------------------------------------------------------------------------" << std::endl;
      std::cout << gen.idhep()  << "[ "
		<< gen.get_ID() << ", "
		<< nchild       << " ] "
		<< std::endl;
    }

    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
      }
      if( daughter.M() > threshold ){
	if( abs(daughter.idhep()) == Dstr0_LUND || abs(daughter.idhep()) == Dstrp_LUND ){
	  dstr_id_list.push_back( daughter.get_ID() );
	  if( fl_message ) std::cout << " --> find D* particle" << std::endl;
	}else{
	  if( fl_message ) std::cout << " --> not D* particle --> continue iterations" << std::endl;
	  search_daughter_Dstr_sub( daughter, dstr_id_list, fl_message, indent+1 );
	}
      }else{
	if( fl_message ) std::cout << " --> low mass particle" << std::endl;
      }
    }
    
    if( fl_message ) std::cout << " ===========> dstr_id_list.size() == " << dstr_id_list.size() << std::endl;
    if( dstr_id_list.size() > 0 ) id1 = dstr_id_list[0];
    if( dstr_id_list.size() > 1 ) id2 = dstr_id_list[1];
    return dstr_id_list.size();
  }

  int DSTRTAUNU::search_daughter_Dstr_sub( Gen_hepevt gen, std::vector<int>& dstr_id_list,
					   bool fl_message, int indent, double threshold ){
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;

    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
      }
      if( daughter.M() > threshold ){
	if( abs(daughter.idhep()) == Dstr0_LUND || abs(daughter.idhep()) == Dstrp_LUND ){
	  dstr_id_list.push_back( daughter.get_ID() );
	  if( fl_message ) std::cout << " --> find D* particle" << std::endl;
	}else{
	  if( fl_message ) std::cout << " --> not D* particle --> continue iterations" << std::endl;
	  search_daughter_Dstr_sub( daughter, dstr_id_list, fl_message, indent+1 );
	}
      }else{
	if( fl_message ) std::cout << " --> low mass particles" << std::endl;
      }
    }
  }

    int DSTRTAUNU::search_daughter_charm( Gen_hepevt gen, int& id1, int& id2, bool fl_message, double threshold ){
    int indent=1;
    std::vector<int> d_id_list; // < gen(D).get_ID >
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;

    if( fl_message ){
      std::cout << "==============================================================================" << std::endl
		<< " idhep() [ get_ID(), nchild ]   --> < search_daughter_charm > "                 << std::endl
		<< "------------------------------------------------------------------------------" << std::endl;
      std::cout << gen.idhep()  << "[ "
		<< gen.get_ID() << ", "
		<< nchild       << " ] "
		<< std::endl;
    }

    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
      }
      if( daughter.M() > threshold ){
	//if( abs(daughter.idhep()) == D0_LUND || abs(daughter.idhep()) == Dplus_LUND || abs(daughter.idhep()) == Dsplus_LUND ){
	if( D_set.count(abs(daughter.idhep())) || ( abs(daughter.idhep())>4000 && abs(daughter.idhep())<4999 )  ){ // 4XXX is charmed baryons
	  d_id_list.push_back( daughter.get_ID() );
	  if( fl_message ) std::cout << " --> find D particle" << std::endl;
	}else{
	  if( fl_message ) std::cout << " --> not D particle --> continue iterations" << std::endl;
	  search_daughter_charm_sub( daughter, d_id_list, fl_message, indent+1 );
	}
      }else{
	if( fl_message ) std::cout << " --> low mass particle" << std::endl;
      }
    }

    if( fl_message ) std::cout << " ===========> d_id_list.size() == " << d_id_list.size() << std::endl;
    if( d_id_list.size() > 0 ) id1 = d_id_list[0];
    if( d_id_list.size() > 1 ) id2 = d_id_list[1];
    return d_id_list.size();
  }

  int DSTRTAUNU::search_daughter_charm_sub( Gen_hepevt gen, std::vector<int>& d_id_list,
					   bool fl_message, int indent, double threshold ){
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;

    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
      }
      if( daughter.M() > threshold ){
	//if( abs(daughter.idhep()) == D0_LUND || abs(daughter.idhep()) == Dplus_LUND || abs(daughter.idhep()) == Dsplus_LUND ){
	if( D_set.count(abs(daughter.idhep())) || ( abs(daughter.idhep())>4000 && abs(daughter.idhep())<4999 )  ){
	  d_id_list.push_back( daughter.get_ID() );
	  if( fl_message ) std::cout << " --> find D particle" << std::endl;
	}else{
	  if( fl_message ) std::cout << " --> not D particle --> continue iterations" << std::endl;
	  search_daughter_charm_sub( daughter, d_id_list, fl_message, indent+1 );
	}
      }else{
	if( fl_message ) std::cout << " --> low mass particles" << std::endl;
      }
    }
  }

  int DSTRTAUNU::search_origin( Gen_hepevt gen, int lund ){
    if( !gen ) return 0;
    //else if( gen.idhep()!= lund ){ // does not work for unknown D
    else if( abs(gen.idhep())!= abs(lund) ){
      int tmp = search_origin( gen.mother(), lund );
      if( tmp ) return tmp;
    }else{
      return gen.get_ID();
    }
    return 0;
  }

  int DSTRTAUNU::search_prompt_lepton( Gen_hepevt gen, int& lund, int& id, bool fl_message ){
    // # of prompt charged tracks
    // id is gen(lepton from semileptonic).get_ID()
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    std::multiset<int> dec_particle_set; // LUND
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;
    if( fl_message ) std::cout << "===========================================" << std::endl
			       << "LUND(" << gen.idhep() << ") has " << nchild  << " childeren" << std::endl
			       << "-------------------------------------------" << std::endl
			       << " idhep() [ get_ID()]  --> < prompt lepton >" << std::endl
			       << "-------------------------------------------" << std::endl;
    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      if( fl_message ) std::cout << i << ": "
				 << std::setw(7) << std::right << daughter.idhep()
				 << " [ "
				 << std::setw(3) << std::right << daughter.get_ID()
				 << " ] ";
      dec_particle_set.insert( daughter.idhep() );

      if( abs(daughter.idhep()) == Electron_LUND
	  || abs(daughter.idhep()) == MUminus_LUND
	  || abs(daughter.idhep()) == TAUminus_LUND ){
	lund = daughter.idhep();
	id   = daughter.get_ID();
	if( fl_message ) std::cout << " -> prompt lepton" << std::endl;
      }else{
	if( fl_message ) std::cout << std::endl;
      }
    }

    int n =0;
    n += dec_particle_set.count( Electron_LUND  );
    n += dec_particle_set.count( Positron_LUND  );
    n += dec_particle_set.count( MUplus_LUND    );
    n += dec_particle_set.count( MUminus_LUND   );
    return n;

  }

  int DSTRTAUNU::semilept_dec( Gen_hepevt gen, int& id, bool fl_message ){
    // 0 : no lepton, 1 : e, 2 : mu, 3 : tau, -1 : other two leptons, -2 : other
    // id is gen(lepton from semileptonic).get_ID()
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    std::multiset<int> dec_particle_set; // LUND
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;

    if( fl_message ) std::cout << "===========================================" << std::endl
			       << "LUND(" << gen.idhep() << ") has " << nchild  << " childeren" << std::endl
			       << "-------------------------------------------" << std::endl
			       << " idhep() [ get_ID()]  --> < semilept_dec > " << std::endl
			       << "-------------------------------------------" << std::endl;
    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      if( fl_message ) std::cout << i << ": "
				 << std::setw(7) << std::right << daughter.idhep()
				 << " [ "
				 << std::setw(3) << std::right << daughter.get_ID()
				 << " ] ";
      dec_particle_set.insert( daughter.idhep() );

      if( abs(daughter.idhep()) == Electron_LUND
	  || abs(daughter.idhep()) == MUminus_LUND
	  || abs(daughter.idhep()) == TAUminus_LUND ){
	id = daughter.get_ID();
	if( fl_message ) std::cout << " -> charged lepton" << std::endl;
      }else{
	if( fl_message ) std::cout << std::endl;
      }
    }

    int n =0;
    n += dec_particle_set.count( Electron_LUND  );
    n += dec_particle_set.count( antiNu_E_LUND  );
    n += dec_particle_set.count( Positron_LUND  );
    n += dec_particle_set.count( Nu_E_LUND      );
    n += dec_particle_set.count( MUplus_LUND    );
    n += dec_particle_set.count( Nu_MU_LUND     );
    n += dec_particle_set.count( MUminus_LUND   );
    n += dec_particle_set.count( antiNu_MU_LUND );
    n += dec_particle_set.count( TAUplus_LUND    );
    n += dec_particle_set.count( Nu_TAU_LUND     );
    n += dec_particle_set.count( TAUminus_LUND   );
    n += dec_particle_set.count( antiNu_TAU_LUND );
    int flag=0;
    if( n==0 ) flag=0;
    else if( n==2 ){
      if(      (dec_particle_set.count(   Electron_LUND ) == 1 || dec_particle_set.count( Positron_LUND ) == 1)
	       && (dec_particle_set.count(   antiNu_E_LUND ) == 1 || dec_particle_set.count(     Nu_E_LUND ) == 1) ) flag = 1; // e
      else if( (dec_particle_set.count(    MUminus_LUND ) == 1 || dec_particle_set.count(   MUplus_LUND ) == 1)
	       && (dec_particle_set.count(  antiNu_MU_LUND ) == 1 || dec_particle_set.count(    Nu_MU_LUND ) == 1) ) flag = 2; // mu
      else if( (dec_particle_set.count(   TAUminus_LUND ) == 1 || dec_particle_set.count(  TAUplus_LUND ) == 1)
	       && (dec_particle_set.count( antiNu_TAU_LUND ) == 1 || dec_particle_set.count(   Nu_TAU_LUND ) == 1) ) flag = 3; // tau
      else{
	id   = 0;
	flag = -1; // other two leptons
      }
    }else{
      id   = 0;
      flag = -2; // other
    }
    if( fl_message ) std::cout << " ===========> semileptonic flag = " << flag << "( " << id << " )" << std::endl;


    return flag;
  }

  int DSTRTAUNU::search_daughter_CC( Gen_hepevt gen, bool fl_message, int indent, double threshold ){
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int nchild = gen.daLast() ==0 ? 0 : gen.daLast() - gen.daFirst() + 1;
    
    if( fl_message && indent==1 ){
      std::cout << "==============================================================================" << std::endl
		<< " idhep() [ get_ID(), nchild ]    --> search_daughter_CC " << std::endl
		<< "------------------------------------------------------------------------------" << std::endl;
      std::cout << gen.idhep()  << "[ "
		<< gen.get_ID() << ", "
		<< nchild       << " ] "
		<< std::endl;
    }
    
    for( int i=0; i<nchild; i++ ){
      Gen_hepevt& daughter( genMgr( Panther_ID(gen.daFirst()+i ) ) );
      int nchildchild = daughter.daLast()==0 ? 0 : daughter.daLast() - daughter.daFirst() + 1;
      if( fl_message ){
	std::cout << std::setw(indent*6) << std::right << " +-> "
		  << daughter.idhep()  << " [ "
		  << daughter.get_ID() << " , "
		  << nchildchild       << " ] ";
      }
      if( daughter.M() > threshold ){
	if( daughter.idhep() == JPsi_LUND       || daughter.idhep() == Psi2S_LUND   || daughter.idhep() == Psi3770_LUND
	    || daughter.idhep() == Psi4040_LUND || daughter.idhep() == Psi4160_LUND || daughter.idhep() == Psi4415_LUND
	    || daughter.idhep() == -7776773     || daughter.idhep() == 7766773 // [bug LUND code] psi(4040)=-7776773, psi(4160)=-7766773
	    ){
	  if( fl_message ) std::cout << " --> find CC(->ll) particle" << std::endl;
	  if     ( daughter.idhep() == -7776773 ) return 9000443; // psi(4040)
	  else if( daughter.idhep() == -7766773 ) return 9010443; // psi(4160)
	  else                                    return daughter.idhep();
	}else{
	  if( fl_message ) std::cout << " --> not CC particle --> continue iterations" << std::endl;
	  int tmp_code = search_daughter_CC( daughter, fl_message, indent+1 );
	  if( tmp_code ) return tmp_code;
	}
      }else{
	if( fl_message ) std::cout << " --> low mass particle" << std::endl;
      }
    }
    return 0;
  }

  int DSTRTAUNU::make_mode_digit( std::multiset<int>& n_particle_set, int digit[] ){
    
    // Judgement of gen_mode(B1)
    int n_total      = n_particle_set.size();
    int n_charged_pi = n_particle_set.count( PIplus_LUND   );
    int n_neutral_pi = n_particle_set.count( PI0_LUND      );
    int n_charged_k  = n_particle_set.count( Kplus_LUND    );
    int n_ks         = n_particle_set.count( Ks_LUND       );   
    int n_electron   = n_particle_set.count( Electron_LUND );
    int n_muon       = n_particle_set.count( MUminus_LUND  );
    
    int n_kl         = n_particle_set.count( Kl_LUND       );
    int n_gamma      = n_particle_set.count( Gamma_LUND    );
    int n_proton     = n_particle_set.count( Proton_LUND   );
    int n_neutron    = n_particle_set.count( Neutron_LUND  );
    int n_lambda     = n_particle_set.count( Lambda_LUND  );
    int n_nu_e       = n_particle_set.count( Nu_E_LUND     );
    int n_nu_mu      = n_particle_set.count( Nu_MU_LUND    );
    int n_nu_tau     = n_particle_set.count( Nu_TAU_LUND   );
    
    int n_rest       = n_total - (n_charged_pi + n_neutral_pi ) - ( n_charged_k + n_ks ) - ( n_electron + n_muon )
      - ( n_kl + n_gamma + n_proton + n_neutron + n_lambda ) - ( n_nu_e + n_nu_mu + n_nu_tau );
    int n_prompt     = n_charged_pi + n_neutral_pi + n_charged_k + n_ks + n_kl + n_proton + n_neutron + n_lambda + n_electron + n_muon;
    
    if( n_charged_pi > 9 ) n_charged_pi = 9;
    if( n_neutral_pi > 9 ) n_neutral_pi = 9;
    if( n_charged_k  > 9 ) n_charged_k  = 9;
    if( n_ks         > 9 ) n_ks         = 9;
    if( n_kl         > 9 ) n_kl         = 9;
    if( n_gamma      > 9 ) n_gamma      = 9;
    if( n_proton     > 9 ) n_proton     = 9;
    if( n_neutron    > 9 ) n_neutron    = 9;
    if( n_lambda     > 9 ) n_lambda     = 9;
    if( n_electron   > 9 ) n_electron   = 9;
    if( n_muon       > 9 ) n_muon       = 9;
    if( n_nu_e       > 9 ) n_nu_e       = 9;
    if( n_nu_mu      > 9 ) n_nu_mu      = 9;
    if( n_nu_tau     > 9 ) n_nu_tau     = 9;
    if( n_rest       > 9 ) n_rest       = 9;
    
    
    int gen_mode_had_r  =                                      1000*n_neutral_pi + 100*n_charged_pi + 10*n_ks       + n_charged_k;
    int gen_mode_had_nr = 100000*n_rest + 10000*n_kl         + 1000*n_lambda + 100*n_proton         + 10*n_neutron  + n_gamma;
    int gen_mode_lepton =                                                                             10*n_electron + n_muon;
    int gen_mode_nu     =                                                      100*n_nu_e           + 10*n_nu_mu    + n_nu_tau;

    
    digit[0] = gen_mode_had_r;
    digit[1] = gen_mode_had_nr;
    digit[2] = gen_mode_lepton;
    digit[3] = gen_mode_nu;
    digit[4] = n_gamma;
    digit[5] = n_prompt;
    return n_rest;
  }

  int DSTRTAUNU::make_flag_semi( int fl1, int fl2 ){
    int flag = 0; // other
    if     (  fl1==1 && fl2==1 )  flag =  -1; //   e-nu :   e-nu
    else if(  fl1==2 && fl2==2 )  flag =  -2; //  mu-nu :  mu-nu
    else if( (fl1==1 && fl2==2) ||
	     (fl1==2 && fl2==1) ) flag =  -3; //   e-nu :  mu-nu
    else if(  fl1==3 && fl2==3 )  flag =   1; // tau-nu : tau-nu
    else if( (fl1==1 && fl2==3)
	     ||
	     (fl1==3 && fl2==1) ) flag =   2; //   e-nu : tau-nu
    else if( (fl1==2 && fl2==3) ||
	     (fl1==3 && fl2==2) ) flag =   3; //  mu-nu : tau-nu

    return flag;
  }

  int DSTRTAUNU::make_flag_dd( int fl1, int fl2 ){
    int flag = 0; // other
    if     (  fl1==1 && fl2==1  ) flag = 1; // DD
    else if( (fl1==1 && fl2==2) ||
	     (fl1==2 && fl2==1) ) flag = 2; // D*D
    else if(  fl1==2 && fl2==2  ) flag = 3; // D*D*
    
    return flag;
  }

  int DSTRTAUNU::which_B( Gen_hepevt gen, int id1, int id2 ){ // 1(B1), 2(B2), -1(other)
    if( mc!=1 ) return -1; // for real-data
    if( !gen ) return -1;
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    if( abs(gen.idhep()) == B0_LUND || abs(gen.idhep()) == Bplus_LUND ){
      if( gen.get_ID() == id1 ){
	return 1;
      }else if( gen.get_ID() == id2 ){
	return 2;
      }else{
	return -1;
      }
    }else{
      int flag = 0;
      if( gen.moFirst() ) flag = which_B( genMgr(Panther_ID(gen.moFirst())), id1, id2 );
      else flag = -1;
      if( flag ) return flag;
    }
 
  }

bool DSTRTAUNU::myCheckSame(const Particle& p , const Mdst_charged& chg)
{
  if(p.nChildren()){
    for(int i=0;i<p.nChildren();++i){
      if(myCheckSame(p.child(i),chg))return true;
    }
  }
  else {
    if(p.mdstCharged()&&p.mdstCharged().get_ID()==chg.get_ID())return true;
  }
  return false;
}

bool DSTRTAUNU::myCheckSame(const Particle& p, const Mdst_gamma& gam)
{
  if(p.nChildren()){
    for(int i=0;i<p.nChildren();++i){
      if(myCheckSame(p.child(i),gam))return true;
    }
  }
  else {
    if(p.mdstGamma()&&p.mdstGamma().get_ID()==gam.get_ID())return true;
  }
  return false;
}

bool DSTRTAUNU::myCheckSame(const Particle& p, const Mdst_pi0& pi0)
{
  if(p.nChildren()){
    for(int i=0;i<p.nChildren();++i){
      if(myCheckSame(p.child(i),pi0.gamma(0))
	 ||myCheckSame(p.child(i),pi0.gamma(1)))return true;
    }
  }
  else {
    if(myCheckSame(p,pi0.gamma(0))
       ||myCheckSame(p,pi0.gamma(1)))return true;
  }
  return false;
}

bool DSTRTAUNU::myCheckSame(const std::vector<Particle>& plist,const Mdst_charged& p){
  for(std::vector<Particle>::const_iterator it=plist.begin();
      it!=plist.end();++it){
    if(myCheckSame(*it,p))return true;
  }
  return false;
}
  
bool DSTRTAUNU::myCheckSame(const Particle& p, const Mdst_ecl& ecl)
{
  if(p.nChildren()){
    for(int i=0;i<p.nChildren();++i){
      if(myCheckSame(p.child(i),ecl))return true;
    }
  }
  else {
    if(p.mdstGamma()&&p.mdstGamma().ecl().get_ID()==ecl.get_ID())return true;
  }
  return false;
}

  double DSTRTAUNU::radgam_energy( Gen_hepevt gen ){
    if( mc!=1 ) return 0; // for real-data
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    double energy = 0;
    if( !gen          ) return 0;
    if( !gen.mother() ) return 0;
    int nchild = gen.mother().daLast()==0 ? 0 : gen.mother().daLast() - gen.mother().daFirst() + 1;
    for( int i=0; i<nchild; i++ ){
      if( genMgr( Panther_ID( gen.mother().daFirst() + i ) ).idhep()==Gamma_LUND ) energy += genMgr( Panther_ID(gen.mother().daFirst() + i) ).E();
    }
    return energy;
  }
  
  int DSTRTAUNU::radgam_cnt( Gen_hepevt gen ){
    if( mc!=1 ) return 0; // for real-data
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    int cnt = 0;
    if( !gen          ) return 0;
    if( !gen.mother() ) return 0;
    int nchild = gen.mother().daLast()==0 ? 0 : gen.mother().daLast() - gen.mother().daFirst() + 1;
    for( int i=0; i<nchild; i++ ){
      if( genMgr( Panther_ID( gen.mother().daFirst() + i ) ).idhep()==Gamma_LUND ) cnt++;
    }
    return cnt;
  }


  void DSTRTAUNU::calKinematics( int id_b, int id_d, int id_lep,
				 double& q2, double& plep, double& w, double& coslep ){
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    Gen_hepevt& B  ( genMgr(Panther_ID(id_b  ) )); HepLorentzVector   B_4V(   B.PX(),   B.PY(),   B.PZ(),   B.E() );
    Gen_hepevt& D  ( genMgr(Panther_ID(id_d  ) )); HepLorentzVector   D_4V(   D.PX(),   D.PY(),   D.PZ(),   D.E() );
    Gen_hepevt& lep( genMgr(Panther_ID(id_lep) )); HepLorentzVector lep_4V( lep.PX(), lep.PY(), lep.PZ(), lep.E() );
    int nchild = B.daLast()==0 ? 0 : B.daLast() - B.daFirst() + 1;
    for( int ichild=0; ichild<nchild; ichild++ ){
      Gen_hepevt& ch( genMgr( Panther_ID(B.daFirst()+ichild) ) );
      if( ch.idhep() == Gamma_LUND ){ // add bremsstrahlug photon to lepton
	HepLorentzVector G_4V( ch.PX(), ch.PY(), ch.PZ(), ch.E() );
	lep_4V += G_4V;
      }
    }
    HepLorentzVector W_4V = B_4V - D_4V;

    //std::cout << "[B] " <<   B_4V << std::endl;
    //std::cout << "[D] " <<   D_4V << std::endl;
    //std::cout << "[L] " << lep_4V << std::endl;
    //std::cout << "[W] " <<   W_4V << std::endl;
    //std::cout << "nchild = " << nchild << std::endl;
    //std::cout << "leplund = " << lep.idhep() << std::endl;
    
    // calculation
    // [q2-plep*]
    q2 = W_4V.m2(); // [q2]
    lep_4V.boost( -(B_4V.boostVector()) ); // boost to B-rest frame
    plep = lep_4V.vect().mag(); // [p_lep*]
    lep_4V.boost(  (B_4V.boostVector()) ); // boost back to lab. frame for w&cos calculation

    // [w-cos]
    w = (B_4V*D_4V)/(B_4V.m()*D_4V.m()); // [w]
    lep_4V.boost( -(W_4V.boostVector()) ); // boost to W-rest frame
    D_4V.boost  ( -(W_4V.boostVector()) ); // boost to W-rest frame
    coslep =  cos( lep_4V.vect().angle(D_4V.vect()) ); // [cos-theta_lep]

    //std::cout << "q2(w) ~ " << B_4V.m2() + D_4V.m2() - 2*B_4V.m()*D_4V.m()*w  << " : " << q2   << std::endl;
    //std::cout << "plep* ~ " << 0.5*(B_4V.m()-D_4V.m()*(w+sqrt(w*w-1)*coslep)) << " : " << plep << std::endl;
    //std::cout << "plep* ~ " << sqrt( pow(0.5*(B_4V.m()-D_4V.m()*(w+sqrt(w*w-1)*coslep)),2) - lep_4V.mag2()) << " : " << plep << std::endl;
    
    return;
  }
#if defined(BELLE_NAMESPACE)
      } //namespace Belle
#endif
