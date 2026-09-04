// ResolutionAnalysis.C

// CLAS12-Specific Headers
#include "CLAS12Reaction.h"
#include "CLAS12DetectorBuilder.h"
#include "clas12defs.h"

// Core RAD Headers
#include "AnalysisManager.h"
#include "KinematicsProcElectro.h"
#include <TBenchmark.h>

void ResolutionAnalysis(const std::string& input_file = "my_clas12_data.hipo") {
    
    ROOT::EnableImplicitMT(4);

    using namespace rad;
    using namespace rad::consts;
    using namespace rad::consts::data_type; 
    using Reaction = clas12::CLAS12Reaction<RHipoDS>;
    using Processor = KinematicsProcElectro;

    // 1. SETUP
    AnalysisManager<Reaction, Processor> mgr{"Res", "events", input_file};
    mgr.SetOutputDir("out_resolutions");
    
    auto& clas12_df = mgr.Reaction();
    clas12_df.SetBeamEnergy(10.6);

    // Initialize both REC and TRUTH streams
    clas12_df.SetupMatching();      

    // 2. PARTICLE CANDIDATES (WITH MC ROLES)
    // The framework needs these integers to map the MC_Lund truth banks to your names!
    const int Role_ScatEle = 3 + 2; 
    const int Role_Proton  = 2 + 2; 
    const int Role_PiP     = 0 + 2; 
    const int Role_PiM     = 1 + 2; 

    clas12_df.SetParticleCandidates(ScatEle(), Role_ScatEle, rad::index::FilterIndices(11), {"rec_pid"});
    clas12_df.SetParticleCandidates("proton",  Role_Proton,  rad::index::FilterIndices(2212), {"rec_pid"});
    clas12_df.SetParticleCandidates("pip",     Role_PiP,     rad::index::FilterIndices(211), {"rec_pid"});
    clas12_df.SetParticleCandidates("pim",     Role_PiM,     rad::index::FilterIndices(-211), {"rec_pid"});  

    clas12_df.MakeCombinations();

    mgr.AddStream(Rec());
    mgr.AddStream(Truth());

    // 2. KINEMATICS RECIPE (Executed on both streams)
    auto topology_recipe = [](Processor& p) {
      //must add particles to reaction map
      p.SetMesonParticles({"pip", "pim"}); 
      p.SetBaryonParticles({"proton"});
   
        std::vector<std::string> tracks = {ScatEle(), "proton", "pip", "pim"};
        p.ParticleP(tracks);
        p.ParticleTheta(tracks);
        p.ParticlePhi(tracks);
    };
    mgr.ConfigureKinematics(topology_recipe);

    // 3. DEFINE CROSS-STREAM RESOLUTIONS (Rec - Truth)
    mgr.CrossStreamDifferences(Rec(), Truth(), {"scat_ele", "proton", "pip", "pim"}, {"pmag", "theta", "phi"});
    
   
    // 5. HISTOGRAM RECIPE (Plotting the Difference Branches)
    auto histogram_recipe = [](rad::histo::Histogrammer& h) {
        // 1D: Raw Resolutions
        h.Create("hResProtonP", "Proton #DeltaP; P_{Rec} - P_{Tru} [GeV]", 100, -0.5, 0.5, "res_proton_pmag");
        h.Create("hResPipTheta", "#pi^{+} #Delta#theta; #theta_{Rec} - #theta_{Tru} [rad]", 100, -0.05, 0.05, "res_pip_theta");
        
        // 2D: Profile against truth momenta
        h.Create2D("hResPvTruP", "Proton Resolution vs True P; P_{Tru} [GeV]; #DeltaP [GeV]", 
                   100, 0, 10, 100, -0.5, 0.5, 
                   "proton_pmag", "res_proton_pmag");
    };
    
    // Plot to the Rec stream to utilize the TrueMatch selection mask
    mgr.ConfigureHistograms(Rec(), histogram_recipe);

    //Print some configuration diagnostics
    mgr.PrintDiagnostics(2); //0=silent, 1 = Verifies data streams and input types, 2 = Dumps all 66+ variables, alias maps, and calculations

    std::cout << "\n=== RUNNING RESOLUTION ANALYSIS ===\n";
    mgr.Snapshot({TruthMatchedCombi()}); 
    mgr.Run();
}
