#include "CLAS12Reaction.h"
#include "ParticleCreator.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RLogger.hxx>
#include <chrono>

inline constexpr std::array<double,4>  rad::beams::InitBotComponents() {return {0,0,0,0.938272};}
inline constexpr std::array<double,4>  rad::beams::InitTopComponents() {return {0,0,10.4,0.000510999};}

void Process_eppippim(){
  
  using namespace rad::names::data_type; //for Rec(), Truth()
  auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
  
  ROOT::EnableImplicitMT(4);
  auto filename = "~/Jlab/clas12/data/hipo/DVPipPimP_006733.hipo";
  std::vector<std::string> files = {filename};
  
  rad::config::CLAS12Reaction rf{files};
  
  rf.AliasColumns();
  std::cout<<" pid type "<< rf.CurrFrame().GetColumnType("rec_pid")<<std::endl;

  // rf.setBeamIonIndex(rad::beams::InitBotFix());
  // rf.setBeamElectronIndex(rad::beams::InitTopFix());
  rf.FixBeamElectronMomentum(0,0,10.4);
  rf.FixBeamIonMomentum(0,0,0);
  
  rf.setScatElectronIndex(rad::indice::useNthOccurance(1,11),{"rec_pid"});
  //rf.DefineVirtualPhoton();
  
  rf.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{"rec_pid"});
  rf.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{"rec_pid"});
  rf.setParticleIndex("proton",rad::indice::useNthOccurance(1,2212),{"rec_pid"});

  rf.Particles().Sum("rho",{"pip","pim"});
  
  rf.setBaryonParticles({"proton"});
  rf.setMesonParticles({"pip","pim"});
  
  //must call this after all particles are configured
  rf.makeParticleMap();
 

//Minimum momentum cut on reconstructed particles
  rf.Filter("(rec_pmag[scat_ele]>0.1)*(rec_pmag[pip]>0.1)*(rec_pmag[pim]>0.1)*(rec_pmag[proton]>0.1)","rec_cut");

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(rf,"W","{scat_ele}");
  rad::rdf::Mass(rf,"Whad","{rho,proton}");
  rad::rdf::Mass(rf,"RhoMass","{rho}");

  //t distribution, column name
  rad::rdf::TBot(rf,"tb");
  rad::rdf::TPrimeBot(rf,"tbp");
  rad::rdf::TTop(rf,"tt");
  rad::rdf::TPrimeTop(rf,"ttp");

  //CM production angles
  rad::rdf::CMAngles(rf,"CM");
  rad::rdf::Q2(rf,"Q2");

  //decay angles
  rad::rdf::gn2s0s0s12::HelicityAngles(rf,"Heli");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",rf};
  histo.Init({Rec()});//will create histograms for mc

  histo.Create<TH1D,double>({"Q2","Q2",500,0,5},{"Q2"});
  histo.Create<TH1D,double>({"W","W",100,0,20.},{"W"});
  histo.Create<TH1D,double>({"RhoMass","M(2#pi) [GeV]",100,0,3},{"RhoMass"});
  histo.Create<TH1D,double>({"tb","t(p,p') [GeV^{2}]",100,-2,5},{"tb"});
  histo.Create<TH1D,double>({"tt","t(g,X) [GeV^{2}]",100,-2,5},{"tt"});
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,float>({"cthHeli","cos(#theta_{hel})",100,-1,1},{"Heli_CosTheta"});
  histo.Create<TH1D,float>({"phCHeli","#phi_{hel})",100,-TMath::Pi(),TMath::Pi()},{"Heli_Phi"});
  histo.Create<TH1D,float>({"EleP","p_{e}",100,0,20},{"pmag[scat_ele]"});
 
  gBenchmark->Start("processing");
  //save all histograms to file
  histo.File("histos/eppippim_histos.root");
 
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
  
}
