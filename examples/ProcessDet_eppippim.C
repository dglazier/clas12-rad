#include "CLAS12DetectorReaction.h"
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

void ProcessDet_eppippim(){
  ///////////////////////////////////////////////////////////
  // Some Preliminaries
  ///////////////////////////////////////////////////////////
 
  using namespace rad::names::data_type; //for Rec(), Truth()
  // auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);  // Log timing etc
  // ROOT::EnableImplicitMT(4); // run multi-core, needs most recent hipo master


  ///////////////////////////////////////////////////////////
  // Setup files to process
  ///////////////////////////////////////////////////////////
  //  auto filename = "~/Jlab/clas12/data/hipo/DVPipPimP_006733.hipo"; //my real data file
  auto filename = "~/Jlab/clas12/data/simulation/RhoFeb24/rho-7221-9*.hipo"; //my simulated file
  std::vector<std::string> files = {filename}; //can add as many files as you wish
  
  ///////////////////////////////////////////////////////////
  // Setup RAD dataframe object. Initialise with files
  ///////////////////////////////////////////////////////////
  rad::clas12::CLAS12DetectorReaction rf{files};
  rf.UseFTB(); //Use ForwardTagger based REC::Particle
  rf.AliasColumnsAndMatchWithMC(); //when using simulated data, mc-match
  //auto pidtype = "tru_pid";

  //rf.AliasColumns(); //when using real data just use REC::Particles
  auto pidtype = "rec_pid";
 
  //Set beam energy. Will eventually remove this whn get rcdb interface
  rf.FixBeamElectronMomentum(0,0,10.4); //default e- mass
  rf.FixBeamIonMomentum(0,0,0); //default p mass

  
  ///////////////////////////////////////////////////////////////
  // Setup final state particles and where to get them
  // useNthOccurance selects the nth particle of a particular pid
  // note the final argument is a PID which creates a particle_OK flag
  // this ==1 if REC::Particle has correct PID ==0 if not
  ///////////////////////////////////////////////////////////////
  rf.setScatElectronIndex(rad::indice::useNthOccurance(1,11),{pidtype});//n=1, pid == 11
  rf.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{pidtype},211);
  rf.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{pidtype},-211);
  rf.setParticleIndex("proton",rad::indice::useNthOccurance(1,2212),{pidtype},2212);


  ///////////////////////////////////////////////////////////
  // Create intermediate particles
  ///////////////////////////////////////////////////////////
  rf.Particles().Sum("rho",{"pip","pim"}); // rho -> pi+ + pi- 

  ///////////////////////////////////////////////////////////
  // Set the particles associated with
  // top (meson) and bottom (baryon)vertices
  ///////////////////////////////////////////////////////////
  rf.setBaryonParticles({"proton"}); //recoil proton
  rf.setMesonParticles({"pip","pim"}); //intermediate meson

  ////  note I could analyse pi- + Delta++ instead
  // rf.setBaryonParticles({"proton","pip"}); //recoil proton
  // rf.setMesonParticles({"pim"}); //intermediate meson

  //must call this after all particles are configured
  rf.makeParticleMap();

  ///////////////////////////////////////////////////////////
  // For debugging I can output particle info to terminal
  ///////////////////////////////////////////////////////////
  //Print mc particles for each event
  //rad::rdf::PrintParticles(rf);
  //Print rec particle for each event
  //rad::rdf::PrintParticles(rf,Rec());
  
  ///////////////////////////////////////////////////////////
  // Create some filter to act on this final state
  // I have columns for indicing the particles in the record
  // these are given by the particle names set above
  // e.g. pip, pim ... also note the scattered electron is scat_ele
  ///////////////////////////////////////////////////////////
  //Minimum momentum cut on reconstructed particles
  //note columns rec_pmag, rec_theta, rec_phi have been created by default
  //I use the particle index to access its value
  rf.Filter("(rec_pmag[scat_ele]>0.1)*(rec_pmag[pip]>0.5)*(rec_pmag[pim]>0.5)*(rec_pmag[proton]>0.1)","rec_cut");
  //EB PID cut on all particles
  rf.Filter("(scat_ele_OK==1) *(pip_OK==1) * (pim_OK==1) * (proton_OK==1)","pid_cut");

  ///////////////////////////////////////////////////////////
  // Specific to CLAS12DetectorReaction
  // I can associate the particles detector information
  // AssociateDetector("HipoBankName",{specific detector IDs},{particle names},{item in bank to associate})
  ///////////////////////////////////////////////////////////
  rf.AssociateDetector("Scintillator",{rad::clas12::FTOF,rad::clas12::CTOF},{"pip"},{"energy","time"});
  rf.AssociateDetector("ForwardTagger",{rad::clas12::FTCAL},{"scat_ele"},{"energy"});
  
  ///////////////////////////////////////////////////////////
  // PErform some kinematic calculations.
  // These come from predefined functions in RAD header files
  // I can also add my own header and functions via a #include
  // at the top of this macro
  ///////////////////////////////////////////////////////////
  //masses column name, {+ve particles}, {-ve particles}
  // #include "BasicKinematicsRDF.h"
  rad::rdf::MissMass(rf,"W","{scat_ele}");
  rad::rdf::Mass(rf,"Whad","{rho,proton}");
  rad::rdf::Mass(rf,"RhoMass","{rho}");

  //t distribution, column name
  //#include "ReactionKinematicsRDF.h"
  rad::rdf::TBot(rf,"tb");
  rad::rdf::TPrimeBot(rf,"tbp");
  rad::rdf::TTop(rf,"tt");
  rad::rdf::TPrimeTop(rf,"ttp");

  //CM production angles
  rad::rdf::CMAngles(rf,"CM");
  rad::rdf::Q2(rf,"Q2");

  //decay angles
  // #include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
  rad::rdf::gn2s0s0s12::HelicityAngles(rf,"Heli");

  ///////////////////////////////////////////////////////////
  // Define histograms using rad::histo::Histogrammer
  // Will create histograms for Rec and Truth variables
  // can split in many kinemtic bins and have histo for each bin
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",rf};
  // we can create many histograms by splitting events into
  // bins, where the same histogram is produced for the given bin
  // e.g. create 10 bins in tru_W between 4 and 54 GeV 
  // histo.Splitter().AddRegularDimension(Truth()+"W", rad::histo::RegularSplits(10,4,14) );
  
  histo.Init({Rec(),Truth()});//will create same histograms for rec and truth variables

  //just got to watch the type given for each variable
  //if you get it wrong it will complain at run time
  //and tell you which type you should have used 
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

  //I have associated detectors so I can plot that information too
  histo.Create<TH2D,float,float>({"ElePvECal","p_{e} versus calorimeter E",100,0,20,100,0,10},{"pmag[scat_ele]","scat_ele_FTCAL_energy"});

  //create another histogrammer for resolutions
  rad::histo::Histogrammer histo_res{"res",rf};
  histo_res.Init();// just create untyped histograms (i.e. will not append Rec or Truth, you must do this)
  histo_res.Create<TH1D,float>({"resEleP","#Delta P_{e'}",100,-1,1},{"res_pmag[scat_ele]"});
  histo_res.Create<TH2D,float,float>({"PVresEleP","P_{e'} v #Delta P_{e'}",100,-1,1,100,0,20},{"res_pmag[scat_ele]","rec_pmag[scat_ele]"});

  ///////////////////////////////////////////////////////////
  // Process by saving all histograms to file
  ///////////////////////////////////////////////////////////
  //save all histograms to file
  histo.File("histos/eppippim_histos.root");
  histo_res.File("histos/eppippim_res_histos.root");
 

  ///////////////////////////////////////////////////////////
  // Process by saving all columns to a tree
  ///////////////////////////////////////////////////////////
  //save tree with all defined branches
  rf.Snapshot("trees/det_eppippim_trees.root");

  
}
