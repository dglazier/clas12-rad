#include "AnalysisManager.h"
#include "CLAS12Reaction.h"
#include "CLAS12DetectorBuilder.h" // New inclusion
#include "clas12defs.h"
#include "KinematicsProcElectro.h"
#include "ElectronScatterKinematics.h"
#include "BasicKinematicsRDF.h"
#include <TBenchmark.h>

/**
 * @brief Analysis Example: CLAS12 e p -> e' p' pi+ pi- (rho)
 * Updated to use the modern CLAS12DetectorBuilder architecture.
 */
void ProcessDet_eppippim() {
  // ROOT::EnableImplicitMT();
  
  using namespace rad;
  using namespace rad::consts::data_type; 
  using Reaction = rad::clas12::CLAS12Reaction;
  using Processor = KinematicsProcElectro;

  gBenchmark->Start("df");

  // =================================================================================
  // 1. SETUP & MATCHING
  // =================================================================================
  std::string filename = "~/Jlab/clas12/data/simulation/RhoFeb24/rho-7221-9*.hipo";
  
  AnalysisManager<Reaction,Processor> mgr{"Rho", "events", filename};
  mgr.SetOutputDir("histos");
  
  auto& clas12_df = mgr.Reaction();

  // clas12_df.InspectBanks({"MC_Lund", "REC_Particle"}, 3);

  // Define Beam Kinematics
  clas12_df.SetBeamEnergy(10.4);
  

  // Turn on FTB Ambiguity Resolution (Prefers RECFT over REC)
  clas12_df.UseFTB(); 
  
  // Setup MC Matching (maps MC_GenMatch into dense SoA arrays)
  clas12_df.SetupMatching();

   // --- Particle Candidates +2 for beams---
  const int Role_ScatEle = 3 + 2; 
  const int Role_Proton  = 2 + 2; 
  const int Role_PiP     = 0 + 2; 
  const int Role_PiM     = 1 + 2; 

  clas12_df.SetParticleCandidates(consts::ScatEle(), Role_ScatEle, rad::index::FilterIndices(11), {"rec_pid"});
  clas12_df.SetParticleCandidates("pip", Role_PiP, rad::index::FilterIndices(211), {"rec_pid"}); 
  clas12_df.SetParticleCandidates("pim", Role_PiM, rad::index::FilterIndices(-211), {"rec_pid"}); 
  clas12_df.SetParticleCandidates("proton", Role_Proton, rad::index::FilterIndices(2212), {"rec_pid"}); 
  // clas12_df.SetParticleCandidates(consts::ScatEle(), Role_ScatEle, rad::index::FilterIndices(11), {"rec_true_pid"});
  // clas12_df.SetParticleCandidates("pip", Role_PiP, rad::index::FilterIndices(211), {"rec_true_pid"}); 
  // clas12_df.SetParticleCandidates("pim", Role_PiM, rad::index::FilterIndices(-211), {"rec_true_pid"}); 
  // clas12_df.SetParticleCandidates("proton", Role_Proton, rad::index::FilterIndices(2212), {"rec_true_pid"}); 

  // --- Generate Combinations ---
  // Must happen before detector building for AutoMap to work
  clas12_df.MakeCombinations();

  // --- Detector Building ---
  // This automatically Extracts banks, synthesizes Regions (FD/CD/FT), 
  // and projects them onto the candidates created above.
  rad::clas12::CLAS12DetectorBuilder det_builder(clas12_df);
  det_builder.BuildAll(); 

  // Add processing streams
  mgr.AddStream(Rec(), "base");
  mgr.AddStream(Truth(), "base");

  // =================================================================================
  // 2. ANALYSIS CONFIGURATION 
  // =================================================================================
   
  auto topology_recipe = [](Processor& p) {
    using namespace consts;
    
    p.Creator().Sum("rho", {{"pip", "pim"}});       
    p.Creator().Sum("Whad_sys", {{"rho", "proton"}});
    p.Creator().Diff("Miss", {{BeamEle(), BeamIon()}, {"rho", "proton", ScatEle()}});
    p.Creator().Diff("W_sys", {{BeamEle(), BeamIon()}, {ScatEle()}});
  
    p.SetMesonParticles({"pip", "pim"}); 
    p.SetBaryonParticles({"proton"});
    
    p.Mass("RhoMass", {"rho"});             
    p.Mass("Whad",    {"Whad_sys"});             
    p.Mass("W",       {"W_sys"});             
    p.Mass2("MissMass2", {"Miss"});             
    p.Q2();
    p.CosThetaCM(); 
    p.PhiCM();       
  
    p.RegisterCalc("tb", rad::physics::TBot);
     
    p.ParticleP({ScatEle(), "pip", "pim", "proton"});
    p.ParticleTheta({ScatEle(), "pip", "pim", "proton"});
    p.ParticlePhi({ScatEle(), "pip", "pim", "proton"});

    // // >>> PASS THROUGH DETECTOR ARRAYS TO THE COMBINATORIAL STREAM <<<
    // // Syntax: PassThrough(ParticleName, RawDetectorArray, OutputSuffix)
    p.PassThrough(ScatEle(), "rec_FT_DetEnergy", "_FT_DetEnergy");
    p.PassThrough("proton", "rec_CD_Time", "_CD_Time");
  };

  mgr.ConfigureKinematics(topology_recipe);

  auto selection_recipe = [](PhysicsSelection& s) {
    using namespace consts;
    s.AddCutMin("ele_p_cut", ScatEle() + "_pmag", 0.1); 
    s.AddCutMin("pip_p_cut", "pip_pmag", 0.3); 
    s.AddCutMin("pim_p_cut", "pim_pmag", 0.3); 
    s.AddCutMin("pro_p_cut", "proton_pmag", 0.1); 
  };

  //mgr.ConfigureSelection(Rec(), selection_recipe);

  auto histogram_recipe = [](histo::Histogrammer& h) {
    using namespace consts;
    h.Create("hQ2",        "Q^{2}; [GeV^2]", 500, 0, 5, "Q2");
    h.Create("hRhoMass",   "M(2#pi) [GeV]", 100, 0, 3, "RhoMass");
    
    // Using the automatically mapped high-level Region variables
    // No brackets, no manual mapping required!
    h.Create2D("hElePvFTCal", "P_{e} vs FT Calorimeter E", 
               100, 0, 10, 100, 0, 10, 
               ScatEle() + "_pmag", ScatEle() + "_FT_DetEnergy");

    h.Create2D("hProtonPvFDTime", "P_{p} vs FD Best Time", 
               100, 0, 10, 100, 0, 50, 
               "proton_pmag", "proton_FD_Time");
  };
 
  // rad::rdf::PrintParticles(clas12_df, Rec());
  // rad::rdf::PrintParticles(clas12_df, Truth());
  //mgr.ConfigureHistograms(histogram_recipe);
  mgr.Snapshot({consts::TruthMatchedCombi()});

  // =================================================================================
  // 3. RUN EVENT LOOP
  // =================================================================================
  // Tell RDataFrame to keep track of the count (Lazy Action)
  // GetBaseFrame() accesses the raw un-filtered HIPO tree
  auto total_events = clas12_df.GetBaseFrame().Count();

  std::cout << "\n=== STARTING CLAS12 ANALYSIS ===\n" << std::endl;
  gBenchmark->Start("analysis");
  mgr.Run();
  gBenchmark->Stop("analysis");
  gBenchmark->Print("analysis");

  std::cout << "\n>>> Total Events Processed: " << *total_events << " <<<\n" << std::endl;

}
