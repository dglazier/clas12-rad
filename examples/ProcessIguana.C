R__LOAD_LIBRARY(libhipo4)
R__LOAD_LIBRARY(libHipoDataFrame)
R__LOAD_LIBRARY(libIguanaServices)
R__LOAD_LIBRARY(libIguanaAlgorithms)

#include "AnalysisManager.h"
#include "CLAS12Reaction.h"
#include "KinematicsProcElectro.h"
#include "hipo4/ThreadedAlgo.hxx"

#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TBenchmark.h>

// IGUANA Algorithms (Same 5 algorithms used in Ex11_Iguana.C)
#include <iguana/algorithms/clas12/ZVertexFilter/Algorithm.h>
#include <iguana/algorithms/clas12/SectorFinder/Algorithm.h>
#include <iguana/algorithms/clas12/rga/MomentumCorrection/Algorithm.h>
#include <iguana/algorithms/clas12/rga/FiducialFilterPass2/Algorithm.h>
#include <iguana/algorithms/physics/InclusiveKinematics/Algorithm.h>

void ProcessIguana() {
  gSystem->Setenv("RCDB_CONNECTION", "sqlite:////home/dglazier/Dropbox/clas12/databases/rcdb_latest.sqlite");
  ROOT::EnableImplicitMT(1);

    using namespace rad;
    using namespace rad::consts;
    using namespace rad::consts::data_type; 
    
    // =================================================================================
    // 1. INITIALIZATION
    // =================================================================================
    // Inject RIguanaDS into the RAD framework
    using Reaction = rad::clas12::CLAS12Reaction<RIguanaDS>;
    using Processor = KinematicsProcElectro;

    std::string filename = "~/Jlab/clas12/data/hipo/DVPipPimP_006733.hipo";
    std::vector<std::string> files={filename,filename,filename,filename};
    AnalysisManager<Reaction, Processor> mgr{"IguanaAnalysis", "events", files};
    
    auto& reaction = mgr.Reaction();
    auto* ds = reaction.GetDataSource();

    // Snapshot the original kinematics BEFORE IGUANA mutates them
    ds->CloneBankTo("REC::Particle", "REC::OriginalParticle");

    // =================================================================================
    // 2. IGUANA CONFIGURATION
    // =================================================================================
    iguana::ThreadedAlgo<iguana::clas12::ZVertexFilter>            algo_vz;
    iguana::ThreadedCreatorAlgo<iguana::clas12::SectorFinder>      algo_sec(ds);
    iguana::ThreadedAlgo<iguana::clas12::rga::FiducialFilterPass2> algo_fidu;
    iguana::ThreadedAlgo<iguana::clas12::rga::MomentumCorrection>  algo_mom;
    iguana::ThreadedCreatorAlgo<iguana::physics::InclusiveKinematics> algo_inc(ds);

    algo_vz.Start();
    algo_sec.Start();
    algo_fidu.Start();
    algo_mom.Start();
    algo_inc.Start();

 
    // Cache indices for extreme performance in the hot loop
    int i_part = ds->GetBankIndex("REC::Particle");
    int i_conf = ds->GetBankIndex("RUN::config");
    int i_cal  = ds->GetBankIndex("REC::Calorimeter");
    int i_traj = ds->GetBankIndex("REC::Traj");
    int i_trk  = ds->GetBankIndex("REC::Track");
    int i_scin = ds->GetBankIndex("REC::Scintillator");
    int i_sec  = ds->GetBankIndex("REC::Particle::Sector");
    int i_inc  = ds->GetBankIndex("physics::InclusiveKinematics");

    ds->SetEventCallback([&](unsigned int slot, hipo::banklist& banks) {
        auto& b_part = banks[i_part];
        auto& b_conf = banks[i_conf];
        auto& b_cal  = banks[i_cal];
        auto& b_traj = banks[i_traj];
        auto& b_trk  = banks[i_trk];
        auto& b_scin = banks[i_scin];
        auto& b_sec  = banks[i_sec];
        auto& b_inc  = banks[i_inc];

        if (!algo_vz[slot].Run(b_part, b_conf)) return false;
        if (!algo_fidu[slot].Run(b_part, b_conf, b_cal, b_traj)) return false;
        if (!algo_sec[slot].Run(b_part, b_trk, b_cal, b_scin, b_sec)) return false;
        if (!algo_mom[slot].Run(b_part, b_sec, b_conf)) return false;
        if (!algo_inc[slot].Run(b_part, b_conf, b_inc)) return false;

        return true;
    });

    // Track which electrons survived the Fiducial cuts
    ds->AddIguanaMask("REC::Particle", "IGUANA_Particle_mask");

    // =================================================================================
    // 3. RAD SETUP (REAL DATA MODE)
    // =================================================================================
    reaction.SetBeamEnergy(10.6);
    reaction.SetupReconstructed();

    // Define the scattered electron candidate without a role ID
    reaction.SetParticleCandidates(ScatEle(), rad::index::FilterIndices(11), {"rec_pid"});
    reaction.MakeCombinations();
    mgr.AddStream(Rec(), "");

    // =================================================================================
    // 4. VARIABLE DEFINITIONS
    // =================================================================================
    // Replicate clas12root inclusive mapping
    reaction.Define("ele_pindex", "physics_InclusiveKinematics_pindex[0]");
    
    reaction.Define("Q2_iguana", "physics_InclusiveKinematics_Q2[0]");
    reaction.Define("W_iguana",  "physics_InclusiveKinematics_W[0]");
    reaction.Define("x_iguana",  "physics_InclusiveKinematics_x[0]");
    reaction.Define("y_iguana",  "physics_InclusiveKinematics_y[0]");
    reaction.Define("vz_iguana", "REC_Particle_vz[ele_pindex]");
  
    // Calculate original momentum and the momentum correction (delta p)
    reaction.Define("ele_p_orig", "sqrt(REC_OriginalParticle_px[ele_pindex]*REC_OriginalParticle_px[ele_pindex] + REC_OriginalParticle_py[ele_pindex]*REC_OriginalParticle_py[ele_pindex] + REC_OriginalParticle_pz[ele_pindex]*REC_OriginalParticle_pz[ele_pindex])");
    reaction.Define("ele_p_mod", "sqrt(rec_px[rec_scat_ele[0]]*rec_px[rec_scat_ele[0]] + rec_py[rec_scat_ele[0]]*rec_py[rec_scat_ele[0]] + rec_pz[rec_scat_ele[0]]*rec_pz[rec_scat_ele[0]])");
    reaction.Define("delta_p", "ele_p_mod - ele_p_orig");

    // =================================================================================
    // 5. STANDARD RDATAFRAME HISTOGRAMMING
    // =================================================================================
    
    // Get the base RDataFrame node and apply filters so arrays don't go out-of-bounds
    auto df = reaction.CurrFrame()
      .Filter("physics_InclusiveKinematics_pindex.size() > 0", "Valid Inclusive Kinematics")
                    .Filter("IGUANA_Particle_mask[ele_pindex] == 1", "Electron Passed Fiducial Cuts");

    // Book the standard ROOT histograms lazily
    auto hQ2_x   = df.Histo2D({"hQ2_x", "Q^{2} vs. x;x;Q^{2} [GeV^{2}]", 100, 0, 1, 100, 0, 12}, "x_iguana", "Q2_iguana");
    auto hQ2_W   = df.Histo2D({"hQ2_W", "Q^{2} vs. W;W [GeV];Q^{2} [GeV^{2}]", 100, 0, 5, 100, 0, 12}, "W_iguana", "Q2_iguana");
    auto hy      = df.Histo1D({"hy", "y distribution;y", 100, 0, 1}, "y_iguana");
    auto hvz     = df.Histo1D({"hvz", "electron v_{z};v_{z} [cm]", 100, -30, 30}, "vz_iguana");
    auto hdeltaP = df.Histo2D({"hdeltaP", "electron momentum correction;p_{meas} [GeV];p_{corr}-p_{meas} [GeV]", 100, 0, 12, 100, -0.2, 0.2}, "ele_p_orig", "delta_p");

    // =================================================================================
    // 6. EXECUTION & DRAWING
    // =================================================================================
    std::cout << "\n=== STARTING IGUANA-RAD ANALYSIS ===\n" << std::endl;
  
    gBenchmark->Start("clas12rad");
    // This will trigger the event loop
    mgr.Run(); 
    // Once the event loop finishes, we can draw the histograms
    gStyle->SetPalette(kRainBow);

    auto* c1 = new TCanvas("c1", "IGUANA Kinematics", 1200, 800);
    c1->Divide(3, 2);

    c1->cd(1); hQ2_x->DrawCopy("COLZ");
    c1->cd(2); hQ2_W->DrawCopy("COLZ");
    c1->cd(3); hy->DrawCopy();
    c1->cd(4); hvz->DrawCopy();
    c1->cd(5); hdeltaP->DrawCopy("COLZ");
   gBenchmark->Stop("clas12rad");
    gBenchmark->Print("clas12rad");

 
    c1->Update();
    c1->SaveAs("IguanaHistograms.pdf");
}
