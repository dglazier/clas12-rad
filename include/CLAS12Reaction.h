/**
 * @file CLAS12Reaction.h
 * @brief Unified Reaction class for CLAS12 analysis using ParticleInjector synchronization.
 * @details 
 * Adheres to the RAD Style Guide: 
 * - Indices/PIDs use 'Indices_t' (ROOT::RVecI).
 * - Kinematics/Detector values use 'ResultType_t' (Double_t).
 * - Explicitly prepends beam particles to indices 0 and 1.
 */

#pragma once

#include "ElectroIonReaction.h"
#include "CLAS12Utilities.h"
#include "CLAS12Names.h"
//#include "clas12defs.h"
#include "ReactionUtilities.h"
#include "Constants.h"
#include "ParticleInjector.h" 
#include "hipo4/RIguanaDS.hxx"

#include <memory>
#include <string>
#include <stdexcept>
#include <iostream>

namespace rad {
namespace clas12 {
    
    using rad::consts::data_type::Rec;
    using rad::consts::data_type::Truth;
    using rad::Indices_t;
    using rad::ResultType_t;

    // =================================================================================
    // CLASS DEFINITION
    // =================================================================================

    template<typename DS_t = RHipoDS>
    class CLAS12Reaction : public rad::ElectroIonReaction {

      private:
        DS_t* _customDS = nullptr; // Our hijacked backdoor pointer

        // Private delegating constructor to intercept the pointer safely
        CLAS12Reaction(DS_t* ds_ptr) 
            : ElectroIonReaction{ ROOT::RDataFrame{std::unique_ptr<DS_t>(ds_ptr)} }, 
              _customDS{ds_ptr} {}
      
    public:
        CLAS12Reaction(const std::string_view treeName, const std::string_view fileNameGlob);
        CLAS12Reaction(const std::string_view treeName, const std::vector<std::string>& filenames);
        CLAS12Reaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames); 
        CLAS12Reaction(ROOT::RDataFrame rdf);
        CLAS12Reaction(ROOT::RDF::RNode rdf);

       // --- Setup Methods ---
        
        /** * @brief Enables Forward Tagger Based (FTB) PID and kinematics resolution. */
        void UseFTB(bool use = true);

        /** * @brief Sets electron beam energy and initializes fixed target kinematics. */
        void SetBeamEnergy(double val);
      
        /** * @brief Sets up unified Reconstructed vectors using ParticleInjector. */
        void SetupReconstructed(Bool_t isEnd = kTRUE);
        
        /** * @brief Sets up unified Truth vectors using ParticleInjector. */
        void SetupTruth(Bool_t isEnd = kTRUE);
        
        /** * @brief Sets up matching metadata between Rec and Truth with +2 offset. */
        void SetupMatching(Bool_t isEnd = kTRUE);

        // --- Detector Association API ---
        
        void DefineDetectorAssociation(const std::string& det, const std::string& item, int subdet, int layer = -1, const std::string& val_bank = "");
        void DefineSimpleAssociation(const std::string& det, const std::string& item);

        bool IsTruthMatched() const { return _truthMatched; }

      /** * @brief Dynamically discovers and prints all columns for the requested HIPO banks. */
      void InspectBanks(const std::vector<std::string>& bankPrefixes, int nEvents = 5);
      /**
       * @brief Retrieves the underlying data source pointer.
       * @details Uses our stored poiner
       */
      DS_t* GetDataSource() {
            if (!_customDS) {
                throw std::runtime_error("CLAS12Reaction::GetDataSource - Data source not initialized via hijacking!");
            }
            return _customDS;
        }
      
    private:
      Int_t _idxBeamEle = 0;       
      Int_t _idxBeamIon = 1;
      bool _truthMatched = false;
      bool _useFTB = false;
      clas12::DetId2Name _detectors; 
    };


    // =================================================================================
    // IMPLEMENTATION
    // =================================================================================

 // =================================================================================
    // IMPLEMENTATION
    // =================================================================================

    template<typename DS_t>
    inline CLAS12Reaction<DS_t>::CLAS12Reaction(const std::string_view treeName, const std::string_view fileNameGlob)
        : CLAS12Reaction(new DS_t(fileNameGlob)) {} // Delegate to the private hijacker

    template<typename DS_t>
    inline CLAS12Reaction<DS_t>::CLAS12Reaction(const std::string_view treeName, const std::vector<std::string>& filenames)
        : CLAS12Reaction(new DS_t(filenames)) {} // Delegate to the private hijacker
  
  template<typename DS_t>
  inline CLAS12Reaction<DS_t>::CLAS12Reaction(const std::string_view treeName, const ROOT::RVec<std::string>& filenames)
    : CLAS12Reaction(new DS_t(std::vector<std::string>(filenames.begin(), filenames.end()))) {}

  template<typename DS_t>
  inline CLAS12Reaction<DS_t>::CLAS12Reaction(ROOT::RDataFrame rdf) 
        : ElectroIonReaction{rdf}, _customDS{nullptr} {}
    
    template<typename DS_t>
    inline CLAS12Reaction<DS_t>::CLAS12Reaction(ROOT::RDF::RNode rdf) 
        : ElectroIonReaction{rdf}, _customDS{nullptr} {}

  
    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::UseFTB(bool use) { _useFTB = use; }

    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::SetBeamEnergy(double val) {
        // Bypass the base class helper methods to avoid redundant index registration.
        // Directly set the protected kinematic members (just like ePICReaction does).
        _p4el_beam = PxPyPzMVector(0.0, 0.0, val, consts::M_ele());
        _p4ion_beam = PxPyPzMVector(0.0, 0.0, 0.0, consts::M_pro());
	_useBeamsFromMC = false;
    }

    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::SetupReconstructed(Bool_t isEnd) {
        AddType(Rec());
        DefineBeamComponents(Rec()); 
        SetBeamElectronIndex(_idxBeamEle, Rec());
        SetBeamIonIndex(_idxBeamIon, Rec());

        rad::ParticleInjector injector(this);
        
        ROOT::RVec<std::string> suffixes = {
            "double px", "double py", "double pz", "double m", 
            "int pid", "int status", "double beta", "double chi2pid", "double vt"
        };
        if (_truthMatched) { suffixes.push_back("int match_id"); }
        injector.DefineParticleInfo(suffixes);

        std::string bEle = Rec() + consts::BeamEle() + "_src_";
        std::string bIon = Rec() + consts::BeamIon() + "_src_";

        ROOT::RVec<std::string> ele_src = {
            bEle+"px", bEle+"py", bEle+"pz", bEle+"m", bEle+"pid", 
            "rad::Indices_t{0}", "rad::RVecResultType{1.0}", "rad::RVecResultType{0.0}", "rad::RVecResultType{0.0}"
        };
        ROOT::RVec<std::string> ion_src = {
            bIon+"px", bIon+"py", bIon+"pz", bIon+"m", bIon+"pid", 
            "rad::Indices_t{0}", "rad::RVecResultType{1.0}", "rad::RVecResultType{0.0}", "rad::RVecResultType{0.0}"
        };
        if(_truthMatched) { ele_src.push_back("rad::Indices_t{0}"); ion_src.push_back("rad::Indices_t{1}"); }
        
        injector.AddSource(Rec(), ele_src);
        injector.AddSource(Rec(), ion_src);

        // --- FIX: RESOLVE FINAL PID FIRST ---
        std::string pid_col = "REC_Particle_pid";
        std::string beta_col = "REC_Particle_beta";

        if (_useFTB && ColumnExists("RECFT_Particle_pid")) {
            Define("REC_Particle_pid_ftb", "rad::clas12::util::MergeFTB(REC_Particle_pid, RECFT_Particle_pid)");
            Define("REC_Particle_beta_ftb", "rad::clas12::util::MergeFTB(REC_Particle_beta, RECFT_Particle_beta)");
            pid_col = "REC_Particle_pid_ftb";
            beta_col = "REC_Particle_beta_ftb";
        }

        // --- FIX: ASSIGN MASSES BASED ON FINAL RESOLVED PID ---
        std::string m_col = Rec() + "m_pdg" + DoNotWriteTag();
        Define(m_col, "rad::util::AssignMasses(" + pid_col + ")");

        ROOT::RVec<std::string> track_src = {
            "REC_Particle_px", 
            "REC_Particle_py", 
            "REC_Particle_pz", 
            m_col,  // <--- Uses the perfectly synchronized PDG mass!
            pid_col, 
            "REC_Particle_status", 
            beta_col, 
            "REC_Particle_chi2pid", 
            "REC_Particle_vt"
        };
        
        if (_truthMatched) { track_src.push_back(Rec() + "match_id_raw" + DoNotWriteTag()); }

        injector.AddSource(Rec(), track_src);
        injector.CreateUnifiedVectors();

        Define(Rec() + "n", Rec() + "px.size()");
        rad::util::CountParticles(this, Rec());
    } 

    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::SetupTruth(Bool_t isEnd) {
        AddType(Truth());
        DefineBeamComponents(Truth());
	SetBeamElectronIndex(_idxBeamEle, Truth());
        SetBeamIonIndex(_idxBeamIon, Truth());

        rad::ParticleInjector injector(this);
        injector.DefineParticleInfo({"double px", "double py", "double pz", "double m", "int pid"});

        std::string bEle = Truth() + consts::BeamEle() + "_src_";
        std::string bIon = Truth() + consts::BeamIon() + "_src_";

        injector.AddSource(Truth(), {bEle+"px", bEle+"py", bEle+"pz", bEle+"m", bEle+"pid"});
        injector.AddSource(Truth(), {bIon+"px", bIon+"py", bIon+"pz", bIon+"m", bIon+"pid"});
        injector.AddSource(Truth(), {"MC_Lund_px", "MC_Lund_py", "MC_Lund_pz", "MC_Lund_mass", "MC_Lund_pid"});
        
        injector.CreateUnifiedVectors();

	Define(Truth() + "n", Truth() + "px.size()");

        rad::util::CountParticles(this, Truth());
    }

    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::SetupMatching(Bool_t isEnd) {
        _truthMatched = true;

        // Map raw HIPO match (Rec track -> MC track) and apply +2 offset for unified arrays
        std::string rawMatch = Rec() + "match_id_raw" + DoNotWriteTag();
        Define(rawMatch,
            [](const ROOT::RVec<short>& pindex, const ROOT::RVec<short>& mcindex, const rad::RVecResultType& raw_px) {
                rad::Indices_t match_id(raw_px.size(), rad::consts::InvalidIndex());
                for (size_t i = 0; i < pindex.size(); ++i) {
                    if (pindex[i] >= 0 && (size_t)pindex[i] < match_id.size()) {
                        // mcindex + 2 maps raw MC track to injected Truth array position
                        match_id[pindex[i]] = mcindex[i] >= 0 ? (int)mcindex[i] + 2 : -1; 
                    }
                }
                return match_id;
            },
            {"MC_GenMatch_pindex", "MC_GenMatch_mcindex", "REC_Particle_px"}
        );

        SetupReconstructed(kFALSE);
        SetupTruth(kFALSE);
        DefineTruePID(Rec());
    }

    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::DefineDetectorAssociation(const std::string& det, const std::string& item, int subdet, int layer, const std::string& val_bank) {
        std::string det_col = "REC_" + det + "_";
        std::string val_col = "REC_" + (val_bank.empty() ? det : val_bank) + "_"; 
        
        std::string outName = Rec() + _detectors.DetName(subdet) + "_" + item; 
        if (layer >= 0) outName += "_L" + std::to_string(layer); 

        // Write the physical column names explicitly into the JIT string expression
        std::string pindex = det_col + "pindex";
        std::string vals   = val_col + item;
        std::string dets   = det_col + "detector";
        std::string n_col  = Rec() + "n";

        if (ColumnExists(det_col + "layer")) {
            std::string layers = det_col + "layer";
            std::string func_call = Form("rad::clas12::util::FilterSubDetectorInfo(%s, %s, %s, %s, %d, %d, %s)", 
                                         pindex.c_str(), vals.c_str(), dets.c_str(), layers.c_str(), subdet, layer, n_col.c_str());
            Define(outName, func_call);
        } else {
            std::string func_call = Form("rad::clas12::util::FilterDetectorInfo(%s, %s, %s, %d, %s)", 
                                         pindex.c_str(), vals.c_str(), dets.c_str(), subdet, n_col.c_str());
            Define(outName, func_call);
        }
    }

    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::DefineSimpleAssociation(const std::string& det, const std::string& item) {
        std::string det_col = "REC_" + det + "_";
        std::string outName = Rec() + det + "_" + item; 

        std::string pindex = det_col + "pindex";
        std::string vals   = det_col + item;
        std::string n_col  = Rec() + "n";

        std::string func_call = Form("rad::clas12::util::FilterSimpleInfo(%s, %s, %s)", 
                                     pindex.c_str(), vals.c_str(), n_col.c_str());
        Define(outName, func_call);
    }
  
    template<typename DS_t>
    inline void CLAS12Reaction<DS_t>::InspectBanks(const std::vector<std::string>& bankPrefixes, int nEvents) {
        std::cout << "\n==================================================================\n";
        std::cout << "=== Inspecting Banks (First " << nEvents << " Events) ===\n";
        
        auto all_cols = GetBaseFrame().GetColumnNames();

        for (const auto& prefix : bankPrefixes) {
            std::vector<std::string> bank_cols;
            for (const auto& col : all_cols) {
                std::string col_str(col); 
                // If the column starts with the requested bank name, grab it!
                if (col_str.find(prefix) == 0) { 
                    bank_cols.push_back(col_str);
                }
            }

            if (bank_cols.empty()) {
                std::cout << "\n[CLAS12Reaction] WARNING: No columns found for bank '" << prefix << "'.\n";
                continue;
            }

            std::cout << "\n>>> DUMPING BANK: " << prefix << " <<<\n";

            // Chunk the columns into groups of 5 to completely bypass ROOT's display truncation
            const size_t chunkSize = 5;
            for (size_t i = 0; i < bank_cols.size(); i += chunkSize) {
                std::vector<std::string> chunk;
                for (size_t j = i; j < i + chunkSize && j < bank_cols.size(); ++j) {
                    chunk.push_back(bank_cols[j]);
                }
                
                auto display = GetBaseFrame().Range(nEvents).Display(chunk, nEvents);
                display->Print();
            }
        }
        std::cout << "==================================================================\n\n";
    }
  
} // namespace clas12
} // namespace rad
