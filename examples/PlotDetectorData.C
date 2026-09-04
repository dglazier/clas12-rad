// CLAS12-Specific Headers
#include "CLAS12Reaction.h"
#include "CLAS12DetectorBuilder.h"
#include "clas12defs.h"

// Core RAD Headers
#include "AnalysisManager.h"
#include "KinematicsProcElectro.h"
#include <TBenchmark.h>

void PlotDetectorData(const std::string& input_file = "my_clas12_data.hipo") {
    
    ROOT::EnableImplicitMT(4 );

    using namespace rad;
    using namespace rad::consts;
    using namespace rad::consts::data_type; 
    
    using Reaction = clas12::CLAS12Reaction<RHipoDS>;
    using Processor = KinematicsProcElectro;

    // =================================================================================
    // 1. SETUP & INITIALIZATION
    // =================================================================================
    AnalysisManager<Reaction, Processor> mgr{"DetectorPlots", "events", input_file};
    mgr.SetOutputDir("out_detectors");
    
    auto& clas12_df = mgr.Reaction();
    clas12_df.SetBeamEnergy(10.6);
    clas12_df.SetupReconstructed(); 

    // =================================================================================
    // 2. PARTICLE CANDIDATES & COMBINATORICS
    // =================================================================================
    clas12_df.SetParticleRecPID(ScatEle(),11);
    clas12_df.SetParticleRecPID("proton", 2212); 
    clas12_df.SetParticleRecPID("pip",    211);  
    clas12_df.SetParticleRecPID("pim",   -211); 

    clas12_df.MakeCombinations();

    // =================================================================================
    // 3. DETECTOR BUILDING 
    // =================================================================================
    // Automatically synthesize all detector mappings and regions natively
    rad::clas12::CLAS12DetectorBuilder det_builder(clas12_df);
    det_builder.BuildAll(); 

    mgr.AddStream(Rec());
    // =================================================================================
    // 4. TOPOLOGY RECIPE 
    // =================================================================================
    auto topology_recipe = [](Processor& p) {
        // Protect particles from memory pruning
        p.SetBaryonParticles({"proton"});
        p.SetMesonParticles({"pip", "pim"});
        
        // Calculate the base kinematics for our 2D plots
        p.ParticleP({ScatEle(), "proton", "pip", "pim"});
        p.ParticleTheta({ScatEle(), "proton", "pip", "pim"});
    };

    // =================================================================================
    // 5. HISTOGRAM RECIPE (DETECTOR EXAMPLES)
    // =================================================================================
    auto histogram_recipe = [](rad::histo::Histogrammer& h) {
        
        // -----------------------------------------------------------------------------
        // CATEGORY 1: Synthesized Regions (FD, CD, FT)
        // -----------------------------------------------------------------------------
        // 1D: Forward Tagger Energy
        h.Create("hEleFTEnergy", "Electron FT Energy; E_{FT} [GeV]", 
                 100, 0, 10, ScatEle() + "_FT_DetEnergy");
        
        // 2D: Momentum vs Central Detector (CTOF/CND) Time
        h.Create2D("hProtonPvCDTime", "Proton: P vs CD Time; P [GeV]; Time_{CD} [ns]", 
                   100, 0, 5, 100, 0, 50, 
                   "proton_pmag", "proton_CD_Time");

        // -----------------------------------------------------------------------------
        // CATEGORY 2: Granular Calorimetry (ECAL Layers)
        // -----------------------------------------------------------------------------
        // 2D: Momentum vs PCAL Energy (Layer 1)
        h.Create2D("hPipPvPCAL", "#pi^{+}: P vs PCAL Energy; P [GeV]; E_{PCAL} [GeV]", 
                   100, 0, 6, 100, 0, 2, 
                   "pip_pmag", "pip_ECAL_energy_L1");

        // -----------------------------------------------------------------------------
        // CATEGORY 3: Granular Scintillators (FTOF Layers)
        // -----------------------------------------------------------------------------
        // 2D: Momentum vs FTOF1B Energy Loss (Layer 2)
        h.Create2D("hPimPvFTOF", "#pi^{-}: P vs FTOF1B dE/dx; P [GeV]; E_{FTOF1B} [MeV]", 
                   100, 0, 6, 100, 0, 25, 
                   "pim_pmag", "pim_FTOF_energy_L2");

        // -----------------------------------------------------------------------------
        // CATEGORY 4: Tracking & Cherenkov (DC, HTCC)
        // -----------------------------------------------------------------------------
        // 1D: Drift Chamber Fit Quality (Chi2)
        h.Create("hProtonDCChi2", "Proton DC Track Fit #chi^{2}; #chi^{2}", 
                 100, 0, 15, "proton_DC_chi2");

        // 2D: Theta Angle vs Number of Photoelectrons (HTCC)
        h.Create2D("hEleThetaVHTCC", "Electron: #theta vs HTCC N_{phe}; #theta [rad]; N_{phe}", 
                   100, 0, 0.6, 100, 0, 50, 
                   ScatEle() + "_theta", ScatEle() + "_HTCC_nphe");

        // -----------------------------------------------------------------------------
        // CATEGORY 5: Covariance Matrix (CovMat)
        // -----------------------------------------------------------------------------
        // 1D: C11 (Variance in tracking parameter 1)
        h.Create("hPipC11", "#pi^{+} CovMat C11; C11 Variance", 
                 100, 0, 0.05, "pip_CovMat_C11");
    };

    // =================================================================================
    // 6. EXECUTION
    // =================================================================================
    mgr.ConfigureKinematics(topology_recipe);
    mgr.ConfigureHistograms(histogram_recipe);

    std::cout << "\n=== STARTING DETECTOR PLOT ANALYSIS ===\n" << std::endl;
    gBenchmark->Start("analysis");
    mgr.Run();
    gBenchmark->Stop("analysis");
    gBenchmark->Print("analysis");
}
