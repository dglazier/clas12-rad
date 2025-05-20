#pragma once

//!  Derived class to configure CLAS12 root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for CLAS12 files with fixed particle order
*/
#include "CLAS12Reaction.h"
#include "CLAS12Utilities.h"
#include "clas12defs.h"
#include "ReactionUtilities.h"
#include "hipo4/RHipoDS.hxx"

namespace rad{
  namespace clas12 {
    
    using rad::names::data_type::Rec;
    //using rad::names::data_type::Truth;

    /**
     *  Function to align detector element with particle entry = index
     *  needs modification for case multiple entries per particle
     */
    template<typename T>
      T ParticleDetInfo(const int index, const ROOT::RVec<ROOT::RVec<Short_t>>& indices, const ROOT::RVec<T>& vals){
      // std::cout<<" ParticleDetInfo "<<index<<" "<<indices<<" "<<vals<<std::endl;
      if(index>indices.size()) return 0;
      if(indices[index].empty()==true) return 0;
      return vals[ indices[index].front() ];
    }
    
    /**
     *  Function to align detector element with particle entry = index
     *  needs modification for case multiple entries per particle
     */
    template<typename T>
      T ParticleSubDetInfo(const int index, const ROOT::RVec<ROOT::RVec<Short_t>>& indices, const ROOT::RVec<T>& vals, const ROOT::RVec<Int_t>& detector, Int_t detid){
      //std::cout<<" ParticleSubDetInfo "<<index<<" "<<indices<<" "<<vals<<detector<<detid<<std::endl;
      if(index>indices.size()) return 0;
      if(indices[index].empty()==true) return 0;
      //only return value if for requested sub detector
      for(auto pentry: indices[index]){
	if(detector[  pentry ] != detid ) continue;
	return vals[ pentry ];
      }
      return 0;
    }
    
    //! Class definition

    class CLAS12DetectorReaction : public CLAS12Reaction {

 

    public:

      CLAS12DetectorReaction(const std::string_view fileNameGlob ) :
	CLAS12Reaction{ fileNameGlob } {

	}
      CLAS12DetectorReaction(const std::vector<std::string> &filenames ) : CLAS12Reaction{ filenames } {
 
      }

	void AssociateDetector(const string& det, const std::vector<int>& subdets, const std::vector<string>& particles, const std::vector<string>& info){

	  std::string det_col{"REC_"};
	  det_col+=det+"_";

	  // Define mapping from detector index to REC::Particle index
	  // many detector hits may go to a single particle
	  auto det_to_rec = det+"_to_rec" + DoNotWriteTag();
  	  Define(det_to_rec ,rad::clas12::ReverseIndexN<short,unsigned long>,{det_col+"pindex",Rec()+"n"});

	  /* //Now map sub-detector label to REC::Particle using det_to_rec */
	  /* [](int subdet, const ROOT::RVec<Short_t>& indices, const ROOT::RVec<T>& vals){ */
	  /*   ROOT::VecOps::Take(vals,indices); */
	  /* } */
	  /* Define(subdet+"_index", */
	  /* 	 Form("rad::clas12::ParticleDetInfo(%s,%s,%s)",particle.data(),(det+"_to_rec").data(), (det_col+"detector").data()) ); */
	  
	  if(IsTruthMatched()){
	    Redefine(det_to_rec,helpers::Rearrange<ROOT::RVec<short>,short>,{det_to_rec,Rec()+"match_id"});
	  }

	  for(const auto& particle:particles){
	    for(const auto& item:info){
	      // Use JIT compilation, can't use templates as need to know at run time
	      // could get the column type at runtime, but this is easier
	      // and PArticleDetInfo is simple so I do not expect slow down
	      // just some JIT overhead
	      
	      /* Define(particle+"_"+det+"_"+item, */
	      /* 	     Form("rad::clas12::ParticleDetInfo(%s,%s,%s)",particle.data(),(det+"_to_rec").data(), (det_col+item).data()) ); */

	      for(const auto& subdet:subdets){
		Define(particle+"_"+ _detectors.DetName(subdet)+"_"+item,
		       Form("rad::clas12::ParticleSubDetInfo(%s,%s,%s,%s,%d)",particle.data(),(det_to_rec).data(), (det_col+item).data(),(det_col+"detector").data(), subdet) );
		std::cout<<"Define particle/detector column : "<<particle+"_"+ _detectors.DetName(subdet)+"_"+item<<std::endl;
	      }
	      
	    }
	  }
	  
	  
	}

    private :
	clas12::DetId2Name _detectors;
	
    };//class def



  }
}
