// CLAS12-Specific Headers
#include "CLAS12Reaction.h"
#include "clas12defs.h"

// Core RAD Headers
#include "AnalysisManager.h"
#include "KinematicsProcElectro.h"

void FilterAnalysis(const std::string& input_file = "my_clas12_data.hipo") {
    
    // Enable ROOT's implicit multithreading
    ROOT::EnableImplicitMT(4);

    using namespace rad;
    using namespace rad::consts;
    using namespace rad::consts::data_type; 
    using Reaction = clas12::CLAS12Reaction<RHipoDS>;
    using Processor = KinematicsProcElectro;

    // =================================================================================
    // 1. SETUP & INITIALIZATION
    // =================================================================================
    AnalysisManager<Reaction, Processor> mgr{"Filter", "events", input_file};
    mgr.SetOutputDir("out_filter");
    
    auto& clas12_df = mgr.Reaction();
    clas12_df.SetBeamEnergy(10.6);

    // Initialize ONLY the reconstructed stream (No MC Truth Matching)
    clas12_df.SetupReconstructed();

    // =================================================================================
    // 2. PARTICLE CANDIDATE DEFINITIONS
    // =================================================================================
    // Using the rationalized API: SetParticlePID(Name, Stream, PDG_Code)
    clas12_df.SetParticlePID(ScatEle(), Rec(), 11);
    clas12_df.SetParticlePID("proton",  Rec(), 2212);
    clas12_df.SetParticlePID("pip",     Rec(), 211);
    clas12_df.SetParticlePID("pim",     Rec(), -211);

    // Build the Combinatorial arrays
    clas12_df.MakeCombinations();

    // Attach the processing stream to the Manager
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
        
        // C. Invariant Masses & Kinematics
        p.Mass("RhoMass", {"rho"});             
        p.Mass("Whad",    {"Whad_sys"});             
        p.Mass("W",       {"W_sys"});             
        p.Mass2("MissMass2", {"Miss"}); // Creates the variable for the exclusivity cut         
        p.Q2();
        p.CosThetaCM(); 
        p.PhiCM();       
        
        p.RegisterCalc("tb", rad::physics::TBot);
         
        p.ParticleP({ScatEle(), "pip", "pim", "proton"});
        p.ParticleTheta({ScatEle(), "pip", "pim", "proton"});
        p.ParticlePhi({ScatEle(), "pip", "pim", "proton"});
    };

    // =================================================================================
    // 4. SELECTION RECIPE (The Filters)
    // =================================================================================
    auto selection_recipe = [](PhysicsSelection& s) {
        // Minimum momentum cuts to remove low-energy background
        s.AddCutMin("ele_p_cut", ScatEle() + "_pmag", 0.1); 
        s.AddCutMin("pip_p_cut", "pip_pmag", 0.3); 
        s.AddCutMin("pim_p_cut", "pim_pmag", 0.3); 
        s.AddCutMin("pro_p_cut", "proton_pmag", 0.1); 
        
        // Exclusivity Cut: Missing Mass Squared strictly around 0
        s.AddCutRange("mm2_cut", "MissMass2", -0.05, 0.05);
    };

    // =================================================================================
    // 5. HISTOGRAM RECIPE
    // =================================================================================
    auto histogram_recipe = [](rad::histo::Histogrammer& h) {
        // Because of the selection mask, these histograms will ONLY fill for events 
        // that passed the momentum and missing mass squared cuts!
        h.Create("hMissMass2", "Missing Mass^{2}; MM^{2} [GeV^{2}];", 200, -0.2, 0.2, "MissMass2");
        h.Create("hRhoMass", "Rho Mass; M_{#pi#pi} [GeV];", 200, 0.2, 1.5, "RhoMass");
        h.Create("hQ2", "Q^{2}; Q^{2} [GeV^{2}];", 200, 0, 5, "Q2");
        
        h.Create2D("hQ2_vs_W", "Q^{2} vs W; W [GeV]; Q^{2} [GeV^{2}]", 
                   100, 0, 5, 100, 0, 5, "W", "Q2");
    };

    // =================================================================================
    // 6. CONFIGURATION & EXECUTION
    // =================================================================================
    mgr.ConfigureKinematics(topology_recipe);
    mgr.ConfigureSelection(selection_recipe);
    mgr.ConfigureHistograms(histogram_recipe);

    // Optional: Print diagnostics before running
    mgr.PrintDiagnostics(2); 

    std::cout << "\n=== RUNNING FILTER ANALYSIS ===\n";
    
    // Snapshot the output. Since we have no truth matching, we leave the args empty!
    mgr.Snapshot(); 
    mgr.Run();
}
