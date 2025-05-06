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
  namespace config{
    using rad::names::data_type::Rec;
    using rad::names::data_type::Truth;

    //! Class definition

    class CLAS12Reaction : public ElectroIonReaction {

 

    public:

      CLAS12Reaction(const std::string_view fileNameGlob ) :
      ElectroIonReaction{ ROOT::RDataFrame{std::move(std::make_unique<RHipoDS>(fileNameGlob))} } {
	/* auto ds  = std::make_unique<RHipoDS>(fileNameGlob); */
        /* auto total_events = ds->GetEntries(); */
 
	/* auto rdf = ROOT::RDataFrame{std::move(ds)}; */
	/* setBaseFrame(rdf); */
	/* setCurrFrame(rdf); */
	/*   std::cout<<"CLAS12Reaction "<<CurrFrame().GetNFiles()<<" "<<rdf.GetNFiles()<<" "<<total_events<<std::endl; */
      }
    CLAS12Reaction(const std::vector<std::string> &filenames ) : ElectroIonReaction{ROOT::RDataFrame{std::move(std::make_unique<RHipoDS>(filenames))}} {
	std::cout<<"CLAS12Reaction "<<CurrFrame().GetNFiles()<<std::endl;
 
      }

      /**
       * Only alias REC::Particles columns
       */ 
       void AliasColumns(Bool_t IsEnd=kTRUE){
	AddType(Rec());

	 setBranchAlias("REC_Particle_px",Rec()+"px");
	 setBranchAlias("REC_Particle_py",Rec()+"py");
	 setBranchAlias("REC_Particle_pz",Rec()+"pz");
	 setBranchAlias("REC_Particle_pid",Rec()+"pid");
	 //create a column for particle masses
	 Define(Rec()+"m",rad::clas12::AssignMasses,{Rec()+"pid"});

	 //reaction::util::CountParticles(this,Rec());
       }


      void PostParticles() override{
	//once particle are added to vectors
	//we can calculate additional components
	AddAdditionalComponents();

	/* if(IsTruthMatched()){ */
	/*   //add resolution functions */
	/*   ResolutionFraction<float>("pmag"); */
	/*   Resolution("theta"); */
	/*   Resolution("phi"); */
	/*   Resolution("eta"); */

	/* } */
      }
      void AddAdditionalComponents(){
	 //and add some additional columns
	 DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
	 DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
	 DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
	 
	 Filter([](const ROOT::RVecF& rth,const ROOT::RVecF& rpmag,const ROOT::RVecF& rph){return true;},{Rec()+"pmag",Rec()+"theta",Rec()+"phi"});
	 
       }
       
    };//CLAS12Reaction

  }
}
