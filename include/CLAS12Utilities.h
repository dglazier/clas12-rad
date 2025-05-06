#pragma once
#include <ROOT/RVec.hxx>


namespace rad{
  namespace clas12 {

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
