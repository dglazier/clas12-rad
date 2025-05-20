#pragma once

//!  Derived class to configure CLAS12 root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for CLAS12 files with fixed particle order
*/
#include "ElectroIonReaction.h"
#include "CLAS12Utilities.h"
#include "ReactionUtilities.h"
#include "hipo4/RHipoDS.hxx"

namespace rad{
  namespace clas12 {
    using rad::names::data_type::Rec;
    using rad::names::data_type::Truth;

    //! Class definition

    class CLAS12Reaction : public rad::config::ElectroIonReaction {

     

    public:

      CLAS12Reaction(const std::string_view fileNameGlob ) :
	rad::config::ElectroIonReaction{ ROOT::RDataFrame{std::move(std::make_unique<RHipoDS>(fileNameGlob))} } {

      }

      CLAS12Reaction(const std::vector<std::string> &filenames ) : rad::config::ElectroIonReaction{ROOT::RDataFrame{std::move(std::make_unique<RHipoDS>(filenames))}} {
 
      }

      void AliasColumns(Bool_t IsEnd=kTRUE);
      void AliasColumnsFTB(Bool_t IsEnd=kTRUE);
      void AliasColumnsMC(Bool_t IsEnd=kTRUE);
      void AliasColumnsAndMC(Bool_t IsEnd=kTRUE);
      void AliasColumnsAndMatchWithMC(Bool_t IsEnd=kTRUE);
      void PostParticles() override;
      void AddAdditionalComponents();
      template<typename T> 
      void RedefineFundamental( const string& name );

 

      void UseFTB(){_isFTBased=true;}
      
      bool IsTruthMatched()const {return _truthMatched;}

    private:

      bool _isFTBased=false;     
      bool _truthMatched =false;
     
    };//CLAS12Reaction


      /////////Class method implementations below
      /**
       * Only alias REC::Particles columns
       */ 
    void CLAS12Reaction::AliasColumns(Bool_t IsEnd){
      if(_isFTBased) return AliasColumnsFTB();
	 
      AddType(Rec());

      setBranchAlias("REC_Particle_px",Rec()+"px");
      setBranchAlias("REC_Particle_py",Rec()+"py");
      setBranchAlias("REC_Particle_pz",Rec()+"pz");
      setBranchAlias("REC_Particle_pid",Rec()+"pid");
      //create a column for particle masses
      Define("REC_Particle_m",rad::clas12::AssignMasses,{Rec()+"pid"});
      //need to alias it for redefines
      setBranchAlias("REC_Particle_m",Rec()+"m");

      setBranchAlias("REC_Particle_status",Rec()+"status");
      setBranchAlias("REC_Particle_vt",Rec()+"vt");
      setBranchAlias("REC_Particle_vx",Rec()+"vx");
      setBranchAlias("REC_Particle_vy",Rec()+"vy");
      setBranchAlias("REC_Particle_vz",Rec()+"vz");
      setBranchAlias("REC_Particle_beta",Rec()+"beta");
      setBranchAlias("REC_Particle_chi2pid",Rec()+"chi2pid");

      //Make a list of good particles
      //INPROCESS
      // Define("REC_Particle_good"[](){});
	 
      reaction::util::CountParticles(this,Rec());
    }
    /**
     * Only alias RECFT::Particles columns
     */ 
    void CLAS12Reaction::AliasColumnsFTB(Bool_t IsEnd){
      //Actually would like to switch to
      //REC::Particle if no RECFT entry, i.e. no FT electron
      //Could make intermediate Defines instead,
      //returning FTB by default
      //and normal if no FT electron
      //then alias the defines to rec_px etc
	 
      AddType(Rec());
      
      Define(Rec()+"n",[](const ROOT::RVecD& px){return px.size();},{"RECFT_Particle_px"});
   
      setBranchAlias("RECFT_Particle_px",Rec()+"px");
      setBranchAlias("RECFT_Particle_py",Rec()+"py");
      setBranchAlias("RECFT_Particle_pz",Rec()+"pz");
      setBranchAlias("RECFT_Particle_pid",Rec()+"pid");
      //create a column for particle masses
      Define("RECFT_Particle_m",rad::clas12::AssignMasses,{Rec()+"pid"});
      //need to alias it for redefines
      setBranchAlias("RECFT_Particle_m",Rec()+"m");

      setBranchAlias("RECFT_Particle_status",Rec()+"status");
      setBranchAlias("RECFT_Particle_vt",Rec()+"vt");
      setBranchAlias("REC_Particle_vx",Rec()+"vx");
      setBranchAlias("REC_Particle_vy",Rec()+"vy");
      setBranchAlias("REC_Particle_vz",Rec()+"vz");
      setBranchAlias("RECFT_Particle_beta",Rec()+"beta");
      setBranchAlias("RECFT_Particle_chi2pid",Rec()+"chi2pid");

      //Make a list of good particles
      //INPROCESS
      // Define("REC_Particle_good"[](){});
    
      reaction::util::CountParticles(this,Rec());
    }

    /**
     * Only alias MC::Lund columns
     */ 
    void CLAS12Reaction::AliasColumnsMC(Bool_t IsEnd){
      AddType(Truth());
	
      setBranchAlias("MC_Lund_px",Truth()+"px");
      setBranchAlias("MC_Lund_py",Truth()+"py");
      setBranchAlias("MC_Lund_pz",Truth()+"pz");
      setBranchAlias("MC_Lund_pid",Truth()+"pid");
      setBranchAlias("MC_Lund_mass",Truth()+"m");
	
      Define(Truth()+"n",Form("rad::helpers::Count(MC_Lund_type,static_cast<short>(1))") );
    }
    /**
     * Alias ReconstructedParticles and MCParticle columns
     */ 
    void CLAS12Reaction::AliasColumnsAndMC(Bool_t IsEnd){
      AliasColumns(kFALSE);
      AliasColumnsMC(kFALSE);

      //columns for mc matching
      setBranchAlias("MC_GenMatch_pindex",Rec()+"match_id");
      setBranchAlias("MC_GenMatch_mcindex",Truth()+"match_id");
      setBranchAlias("MC_GenMatch_quality",Truth()+"qual");
    }
    /**
     * Alias the columns and rearrange entries 
     * according to MC::Match
     * this reorders reconstructed to match mc
     */
    void CLAS12Reaction::AliasColumnsAndMatchWithMC(Bool_t IsEnd){

      AliasColumnsAndMC(kFALSE);
      _truthMatched = true;
	
      if(IsEnd){
	reaction::util::RedefineFundamentalAliases(this);

      }
    }
    /**
     * Particles have been reordered etc, ready to do some calculations
     * here we add some useful momentum components and resolutions
     */
    void CLAS12Reaction::PostParticles() {
      //once particle are added to vectors
      //we can calculate additional components
      AddAdditionalComponents();

      if(IsTruthMatched()){
	reaction::util::ResolutionFraction(this,"pmag");
	reaction::util::Resolution(this,"theta");
	reaction::util::Resolution(this,"phi");
      }
    }
    /**
     * Calculate spherical momentum components
     */
    void CLAS12Reaction::AddAdditionalComponents(){
      //and add some additional columns
      DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
      DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
      DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
	 
	 
    }

    /**
     * Reorder REC::Particles to match MC::Lund
     */
    template<typename T> 
    void CLAS12Reaction::RedefineFundamental( const string& name ){
	
      auto contains = [](const std::string&  s1,const std::string& s2){
	return (s1.find(s2) != std::string::npos);
      };
      if(contains(name,Rec()+"match_id") ){return;}//don't reorder our order!
      if(contains(name,Truth()+"match_id") ){return;}//don't reorder our order!
	  
      if(contains(name,Rec()) ){
	RedefineViaAlias(name,helpers::Reorder<T,short,short>,{name.data(),Rec()+"match_id",Truth()+"match_id",Truth()+"n"});
      }

    }
      
  }//config
}//rad
