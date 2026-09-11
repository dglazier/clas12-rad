// CLAS12-Specific Headers
#include "CLAS12Reaction.h"
#include "clas12defs.h"

// Core RAD Headers
#include "AnalysisManager.h"
#include "KinematicsProcElectro.h"

void HistogramAnalysis(const std::string& input_file = "my_clas12_data.hipo") {
    
    // Enable multithreading for fast processing
    ROOT::EnableImplicitMT(4);

    using namespace rad;
    using namespace rad::consts;
    using namespace rad::consts::data_type; 
    using Reaction = clas12::CLAS12Reaction<RHipoDS>;
    using Processor = KinematicsProcElectro;

    // =================================================================================
    // 1. SETUP & INITIALIZATION
    // =================================================================================
    AnalysisManager<Reaction, Processor> mgr{"HistoSplit", "events", input_file};
    mgr.SetOutputDir("out_histograms");
    
    auto& clas12_df = mgr.Reaction();
    clas12_df.SetBeamEnergy(10.6);

    // Initialize ONLY the reconstructed stream (No MC Truth needed here)
    clas12_df.SetupReconstructed();

    // =================================================================================
    // 2. PARTICLE CANDIDATE DEFINITIONS
    // =================================================================================
    clas12_df.SetParticlePID(ScatEle(), Rec(), 11);
    clas12_df.SetParticlePID("proton",  Rec(), 2212);
    clas12_df.SetParticlePID("pip",     Rec(), 211);
    clas12_df.SetParticlePID("pim",     Rec(), -211);

    clas12_df.MakeCombinations();
    mgr.AddStream(Rec());

    // =================================================================================
    // 3. TOPOLOGY & KINEMATICS RECIPE
    // =================================================================================
    auto topology_recipe = [](Processor& p) {
        
        // A. Combinatorial Sums and Differences
        p.Creator().Sum("rho", {{"pip", "pim"}});       
        p.Creator().Sum("Whad_sys", {{"rho", "proton"}});
        p.Creator().Diff("Miss", {{BeamEle(), BeamIon()}, {"rho", "proton", ScatEle()}});
        p.Creator().Diff("W_sys", {{BeamEle(), BeamIon()}, {ScatEle()}});
        
        // B. Physics Groupings
        p.SetMesonParticles({"pip", "pim"}); 
        p.SetBaryonParticles({"proton"});
        
        // C. Calculate Variables (These become our base column names)
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
    };

   
    // =================================================================================
    // 5. HISTOGRAM RECIPE (With Splitting!)
    // =================================================================================
    auto histogram_recipe = [](rad::histo::Histogrammer& h) {
        
        // A. DEFINE THE SPLIT AXIS
        // Syntax: AddSplit(SplitName, TargetCol, Bins, Min, Max)
        // This sets up 4 discrete Q2 bins between 1.0 and 5.0 GeV^2.
        // NOTE: We drop the "rec_" prefix for "Q2" here!
      h.AddSplit("Q2_Bin", "Q2", 4, 0, 1.0);

        // B. BOOK THE HISTOGRAMS
        // Every histogram below will automatically be generated 4 times (once per Q2 bin).
        // Because of the selection recipe, these will ONLY fill for events passing our cuts.
        
        // Exclusivity and Resonance
        h.Create("hMissMass2", "Missing Mass^{2}; MM^{2} [GeV^{2}];", 100, -0.1, 0.1, "MissMass2");
        h.Create("hRhoMass", "Rho Mass; M_{#pi#pi} [GeV];", 100, 0.2, 1.5, "RhoMass");
        
        // High-Level Kinematics
        h.Create("hW", "W; W [GeV];", 100, 1.5, 4.0, "W");
        h.Create("htb", "Momentum Transfer; -t [GeV^{2}];", 100, 0.0, 3.0, "tb");
        
        // Track Kinematics
        h.Create("hPipP", "#pi^{+} Momentum; P [GeV/c];", 100, 0, 10, "pip_pmag");
        
        // 2D Profiles
        h.Create2D("hElePvTheta", "Scattered Electron P vs #theta; #theta [rad]; P [GeV/c]", 
                   100, 0, 1.0, 100, 0, 10.6, "scat_ele_theta", "scat_ele_pmag");
    };

    // =================================================================================
    // 6. CONFIGURATION & EXECUTION
    // =================================================================================
    mgr.ConfigureKinematics(topology_recipe);
    mgr.ConfigureHistograms(histogram_recipe);

    // Optional: Print the routing diagnostics to verify the setup
    // mgr.PrintDiagnostics(); 

    std::cout << "\n=== RUNNING HISTOGRAM ENGINE ===\n";
    
    // We strictly call Run() without Snapshot()! 
    // This executes the graph solely to populate and save the split histograms.
    mgr.Run();
}
