#pragma once
#include <ROOT/RVec.hxx>


namespace rad{
  namespace clas12 {
    /**
     * Rearrange vec in order of its own elements, with additional size for missing elements
     * missing elements will have index -1 and users will have to deal with this 
     * in subsequent functions
     * e.g. [1,3,2,0,0] , 6 -> [[3,4],[0],[2],[1],[],[]]
     */
    //detector matrix [det][subdet][layer]
    using detector_matrix_t = ROOT::VecOps::RVec<ROOT::VecOps::RVec<short>>;
    template<typename T,typename Tn>
     detector_matrix_t ReverseIndexN(const ROOT::VecOps::RVec<short>& vec,Tn nentries){
      //std::cout<<"ReverseIndexN "<< vec <<" "<<vec.size()<<" "<<nentries<<std::endl;
      if(nentries==0)return detector_matrix_t();
      if( nentries<vec.size() ) nentries = vec.size() ;//no procedure for choosing what to remove
      detector_matrix_t  result(nentries);//unfilled elements will be = -1
      //std::cout<<"done init"<<std::endl;
      T entry = 0;
      for(auto idx:vec){
	//	result[idx] = entry;
	result[idx].push_back(entry);
	++entry;
      }
      // std::cout<<result<<std::endl;
      return result;
    }
    
    ///////////////////////////////////////////////////////
    constexpr double PdgToMass(int pdg){

      switch ( pdg ) {
      case 11 :
	return 0.00051099900;
      case -11 :
	return 0.00051099900;
      case 211 :
	return 0.13957040;
      case -211 :
	return 0.13957040;
      case 321 :
	return 0.49367700;
      case -321 :
	return 0.49367700;
      case 2212 :
	return 0.93827210;
      case -2212 :
	return 0.93827210;
      case 2112 :
	return 0.93956540;
      case 22 :
	return 0.;
      case 45: //CLAS12 deuteron
	return 1.875612;

	
      default :
	return 0.;
      }

      
    }//PdgToMass

    ///////////////////////////////////////////////////////
    ROOT::RVecD AssignMasses( const ROOT::RVecI &pid){
      auto n = pid.size();
      ROOT::RVecD masses(n);
      for(size_t i=0;i<n;++i){
	masses[i] = PdgToMass(pid[i]);
      }
      return masses;
    }

    
  }//clas12
}//rad
