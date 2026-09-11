// CLAS12-Specific Headers
#include "CLAS12Reaction.h"
#include "clas12defs.h"

// Core RAD Headers
#include "AnalysisManager.h"
#include "KinematicsProcElectro.h"
#include "ElectronScatterKinematics.h"


void MinimalExample(const std::string& input_file = "my_clas12_data.hipo") {
    
    // Enable ROOT's implicit multithreading for fast, parallel execution
    ROOT::EnableImplicitMT(4);

    using namespace rad;
    using namespace rad::consts;
    using namespace rad::consts::data_type; 
 
    // Define the Reaction and Processor types for a CLAS12 Electroproduction analysis
    using Reaction = clas12::CLAS12Reaction<RHipoDS>;
    using Processor = KinematicsProcElectro;

    // =================================================================================
    // 1. SETUP & INITIALIZATION
    // =================================================================================
    AnalysisManager<Reaction, Processor> mgr{"eppippim", "events", input_file};
    mgr.SetOutputDir("out_minimal");
    
    auto& clas12_df = mgr.Reaction();
    clas12_df.SetBeamEnergy(10.6);

    // 1. Explicitly register the "rec_" data type to the Reaction object
    clas12_df.SetupReconstructed();


    // =================================================================================
    // 2. PARTICLE CANDIDATE DEFINITIONS
    // =================================================================================
    clas12_df.SetParticleRecPID(ScatEle(), 11);
    clas12_df.SetParticleRecPID("pip", 211 );
    clas12_df.SetParticleRecPID("pim",-211);
    clas12_df.SetParticleRecPID("proton",2212 );
    
    // 2. Build the Combinatorial arrays
    clas12_df.MakeCombinations();

    // 3. Attach the processing stream to the Manager.
    // This instantiates the KinematicsProcElectro, which requires the types to be 
    // registered and the combinations to already be made!
    mgr.AddStream(Rec());
   
  // =================================================================================
    // 3. TOPOLOGY & KINEMATICS RECIPE
    // =================================================================================
    // Define the physics topology inside a lambda. The framework injects the KinematicsProcessor.
    auto topology_recipe = [](Processor& p) {
        
        // A. Define "Missing" 4-Vector
        // Equation: P_miss = (P_beam_ele + P_target_ion) - (P_scat_ele + P_proton + P_pip + P_pim)
        p.Creator().Diff("miss", {
            {BeamEle(), BeamIon()},                      // Initial State
            {ScatEle(), "proton", "pip", "pim"}          // Final State
        });
        
        // B. Calculate Variables
        // Calculate the invariant mass of the newly created "miss" 4-vector.
        // This registers a new column named "MissMass" in the flat output tree.
        p.Mass("MissMass", {"miss"});
        
        // (Optional) Calculate standard electroproduction variables (Q2, W, xbj)
        p.Q2();
        
        // (Optional) Calculate the invariant mass of the rho(770) meson candidate
        p.Mass("RhoMass", {"pip", "pim"});
    };

    // =================================================================================
    // 4. HISTOGRAM RECIPE
    // =================================================================================
    // Define histograms using the registered variable names from the Topology Recipe.
    auto histogram_recipe = [](rad::histo::Histogrammer& h) {
        
        // 1D Histogram: Missing Mass
        // Signature: Name, Title; X-axis; Y-axis, Bins, Min, Max, Target Column
        h.Create("hMissMass", "Missing Mass; MM_{e'p'#pi^{+}#pi^{-}} [GeV/c^{2}]", 
                 200, -0.5, 1.5, "MissMass");
                 
        // 1D Histogram: Rho Mass
        h.Create("hRhoMass", "Invariant Mass(#pi^{+}#pi^{-}); M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", 
                 200, 0.2, 1.5, "RhoMass");
                 
        // 2D Histogram: Q2 vs Missing Mass
        h.Create2D("hQ2_vs_MM", "Q^{2} vs Missing Mass; MM [GeV/c^{2}]; Q^{2} [GeV^{2}]", 
                   200, -0.5, 1.5, 200, 0, 10, "MissMass", "Q2");
    };

    // =================================================================================
    // 5. CONFIGURATION & EXECUTION
    // =================================================================================
    // Inject the recipes into the Analysis Manager
    mgr.ConfigureKinematics(topology_recipe);
    mgr.ConfigureHistograms(histogram_recipe);

    // (Optional) Tell the framework to save the calculated variables to a flat TTree
    mgr.Snapshot(); 

    //(Optional) Print some configuration diagnostics
    mgr.PrintDiagnostics(0); //0=silent, 1 = Verifies data streams and input types, 2 = Dumps all 66+ variables, alias maps, and calculations

    // Run the RDataFrame event loop (Lazy Evaluation triggered here)
    std::cout << "\n[Minimal_eppippim] Commencing Event Loop..." << std::endl;
    mgr.Run();
    std::cout << "[Minimal_eppippim] Analysis Complete!" << std::endl;
}
