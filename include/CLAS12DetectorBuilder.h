/**
 * @file CLAS12DetectorBuilder.h
 * @brief Automates the extraction and combinatorial synchronization of all CLAS12 detector banks.
 * @details 
 * Implements the Fluent / Lazy Evaluation design pattern. It provides modular access to every 
 * leaf of the CLAS12 tracking, PID, and calorimetry systems. 
 * * Features:
 * 1. **Dynamic Geometry:** Uses the info::DetectorLayers() registry to automatically iterate over 
 * the correct layout of the spectrometer without hardcoded layer strings.
 * 2. **Zero-Cost Abstraction:** Because RDataFrame is lazy, defining all possible detector leaves 
 * incurs zero CPU cost unless the analyst explicitly uses them in a plot or cut.
 * 3. **Combinatorial Syncing:** Automatically projects fundamental detector arrays (N_Tracks) 
 * onto combinatorial candidates (N_Combinations) using `ROOT::VecOps::Take`.
 * 4. **Region Abstractions:** Synthesizes `region_particle`-style variables (e.g. `FD_Time`, `CD_DetEnergy`) 
 * by replicating the object-oriented fallback logic dynamically using variadic SIMD vectorization.
 */

#pragma once

#include "CLAS12Reaction.h"
#include "CLAS12Utilities.h"
#include "clas12defs.h"
#include "CLAS12Names.h"
#include "Constants.h"

#include <vector>
#include <string>
#include <iostream>

namespace rad {
namespace clas12 {

    // =================================================================================
    // CLASS DECLARATION
    // =================================================================================

    template<typename DS_t = RHipoDS>
    class CLAS12DetectorBuilder {
    public:
        /**
         * @brief Constructor.
         * @param rxn Reference to the main CLAS12Reaction object.
         */
        explicit CLAS12DetectorBuilder(CLAS12Reaction<DS_t>& rxn);

        // --- High-Level Builders ---

        /**
         * @brief The master method: Extracts ALL detector banks and builds ALL high-level regions.
         * @details MUST be called AFTER `MakeCombinations()`.
         */
        void BuildAll();

        /**
         * @brief Synthesizes `region_particle`-style arrays (FD, CD, FT, BAND).
         * @details Applies physical fallback logic (e.g. CTOF > CND1 > CND2) dynamically using 
         * string builders to generate compound variables like `CD_Time` and `FD_DetEnergy`.
         */
        void BuildRegions();

        /**
         * @brief Aliases global scalars (Event info, Run config, FTB start time).
         */
        void BuildEventBanks();

        // --- Low-Level Granular Builders ---

        void BuildCalorimeter(int subdet, int layer = -1);
        void BuildScintillator(int subdet, int layer = -1);
        void BuildTracker(int subdet);
        void BuildTrajectories(int subdet, int layer = -1);
        void BuildCherenkov(int subdet);
        void BuildForwardTagger(int subdet, int layer = -1);
        void BuildCovMatrix();

    private:
        CLAS12Reaction<DS_t>& _rxn;
        clas12::DetId2Name _detectorNames;

        /** @brief Safely registers a column association without throwing on re-definitions. */
        void SafeDefineAssoc(const std::string& det, const std::string& item, int subdet, int layer = -1, const std::string& val_bank = "");
        
        /** @brief Safely registers simple index-only banks (e.g. CovMatrix). */
        void SafeDefineSimple(const std::string& det, const std::string& item);

        /** @brief Synchronizes a fundamental subdetector array to combinatorial candidates. */
        void AutoMap(int subdet, const std::string& item, int layer);

        /** @brief Synchronizes a generic auxiliary array to combinatorial candidates. */
        void AutoMapSimple(const std::string& prefix, const std::string& item);
    };

} // namespace clas12
} // namespace rad

// =================================================================================
// IMPLEMENTATION
// =================================================================================

namespace rad {
namespace clas12 {

    template<typename DS_t>
    inline CLAS12DetectorBuilder<DS_t>::CLAS12DetectorBuilder(CLAS12Reaction<DS_t>& rxn) : _rxn(rxn) {}

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildAll() {
        BuildEventBanks();
        BuildCovMatrix();
        
        // 1. Dynamic extraction driven entirely by the geometry registry!
        for (int lay : info::DetectorLayers().at(clas12::ECAL)) BuildCalorimeter(clas12::ECAL, lay);
        for (int lay : info::DetectorLayers().at(clas12::FTOF)) BuildScintillator(clas12::FTOF, lay);
        for (int lay : info::DetectorLayers().at(clas12::CND))  BuildScintillator(clas12::CND, lay);
        for (int lay : info::DetectorLayers().at(clas12::DC))   BuildTrajectories(clas12::DC, lay);
        for (int lay : info::DetectorLayers().at(clas12::CVT))  BuildTrajectories(clas12::CVT, lay);
        
        // Initialize BAND if it exists in the registry
        if (info::DetectorLayers().count(clas12::BAND)) {
            for (int lay : info::DetectorLayers().at(clas12::BAND)) BuildScintillator(clas12::BAND, lay);
        }

        BuildTracker(clas12::DC);
        BuildTracker(clas12::CVT);
        
        // 2. Layerless / Non-segmented detectors
        BuildScintillator(clas12::CTOF);
        BuildCherenkov(clas12::HTCC);
        BuildCherenkov(clas12::LTCC);
        BuildCherenkov(clas12::RICH);
        
        BuildForwardTagger(clas12::FTCAL);
        BuildForwardTagger(clas12::FTHODO);
        
        // 3. Build compound physics abstractions
        BuildRegions(); 
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildRegions() {
        std::string rec = rad::consts::data_type::Rec();
        std::string dnw = rad::DoNotWriteTag();

        // Dynamically resolve names
        std::string ecal = _detectorNames.DetName(clas12::ECAL);
        std::string ftof = _detectorNames.DetName(clas12::FTOF);
        std::string cnd  = _detectorNames.DetName(clas12::CND);
        std::string ctof = _detectorNames.DetName(clas12::CTOF);
        std::string ftcal = _detectorNames.DetName(clas12::FTCAL);
        std::string fthodo = _detectorNames.DetName(clas12::FTHODO);
        std::string trkName = _detectorNames.DetName(clas12::DC);

        // --- FORWARD DETECTOR (FD) ---
        if (!_rxn.ColumnExists(rec + "FD_DetEnergy")) {
            
            // Energy Sums
            _rxn.Define(rec + "FD_DetEnergy", 
                util::BuildLayerFunctionString("rad::clas12::util::SumValid", rec, ecal, "energy", 
                    info::DetectorLayers().at(clas12::ECAL))); 

            _rxn.Define(rec + "FD_DeltaEnergy", 
                util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ftof, "energy", 
                    {clas12::FTOF1B, clas12::FTOF1A, clas12::FTOF2})); 
            
            // Timing (Fallback arrays)
            std::string tFtof = rec + "FTOF_time_best" + dnw;
            std::string tEcal = rec + "ECAL_time_best" + dnw;

            _rxn.Define(tFtof, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ftof, "time", 
                               {clas12::FTOF1B, clas12::FTOF1A, clas12::FTOF2}));
            _rxn.Define(tEcal, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ecal, "time", 
                               {clas12::PCAL, clas12::ECOUT, clas12::ECIN}));

            // Explicit logic to enforce charge!=0 FTOF physics (with InvalidEntry handling)
            _rxn.Define(rec + "FD_Time", 
                [](const rad::RVecResultType& tof, const rad::RVecResultType& cal, const rad::Indices_t& pids) {
                    rad::RVecResultType res(tof.size(), rad::consts::InvalidEntry<rad::ResultType_t>());
                    for (size_t i = 0; i < tof.size(); ++i) {
                        bool isNeutral = (pids[i] == 22 || pids[i] == 2112 || pids[i] == 111 || pids[i] == 0);
                        if (!isNeutral && !rad::consts::IsInvalidEntry(tof[i])) res[i] = tof[i];
                        else if (!rad::consts::IsInvalidEntry(cal[i]))          res[i] = cal[i];
                    }
                    return res;
                }, {tFtof, tEcal, rec + "pid"}
            );

            // Path 
            std::string pFtof = rec + "FTOF_path_best" + dnw;
            std::string pEcal = rec + "ECAL_path_best" + dnw;

            _rxn.Define(pFtof, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ftof, "path", 
                               {clas12::FTOF1B, clas12::FTOF1A, clas12::FTOF2}));
            _rxn.Define(pEcal, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ecal, "path", 
                               {clas12::PCAL, clas12::ECOUT, clas12::ECIN}));

            _rxn.Define(rec + "FD_Path", 
                [](const rad::RVecResultType& pTof, const rad::RVecResultType& pCal, const rad::Indices_t& pids) {
                    rad::RVecResultType res(pTof.size(), rad::consts::InvalidEntry<rad::ResultType_t>());
                    for (size_t i = 0; i < pTof.size(); ++i) {
                        bool isNeutral = (pids[i] == 22 || pids[i] == 2112 || pids[i] == 111 || pids[i] == 0);
                        if (!isNeutral && !rad::consts::IsInvalidEntry(pTof[i])) res[i] = pTof[i];
                        else if (!rad::consts::IsInvalidEntry(pCal[i]))          res[i] = pCal[i];
                    }
                    return res;
                }, {pFtof, pEcal, rec + "pid"}
            );

            // Sector
            _rxn.Define(rec + "FD_Sector", util::BuildFunctionString("rad::clas12::util::Fallback", 
                {rec + trkName + "_sector", util::ColName(rec, ftof, "sector", clas12::FTOF1B), util::ColName(rec, ecal, "sector", clas12::PCAL)}));
        }

        // --- CENTRAL DETECTOR (CD) ---
        if (!_rxn.ColumnExists(rec + "CD_DetEnergy")) {
            
            _rxn.Define(rec + "CD_DetEnergy", 
                util::BuildLayerFunctionString("rad::clas12::util::SumValid", rec, cnd, "energy", 
                    info::DetectorLayers().at(clas12::CND)));

            _rxn.Define(rec + "CD_DeltaEnergy", rec + ctof + "_energy");
            
            std::string tCnd = rec + "CND_time_best" + dnw;
            _rxn.Define(tCnd, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, cnd, "time", 
                              info::DetectorLayers().at(clas12::CND)));
            
            // Explicit, Regex-proof C++ Callable for CD_Time
            _rxn.Define(rec + "CD_Time", 
                [](const rad::RVecResultType& tCtof, const rad::RVecResultType& tCnd) {
                    rad::RVecResultType res(tCtof.size(), rad::consts::InvalidEntry<rad::ResultType_t>());
                    for (size_t i = 0; i < tCtof.size(); ++i) {
                        if (!rad::consts::IsInvalidEntry(tCtof[i]))      res[i] = tCtof[i];
                        else if (!rad::consts::IsInvalidEntry(tCnd[i]))  res[i] = tCnd[i];
                    }
                    return res;
                }, {rec + ctof + "_time", tCnd}
            );

            std::string pCnd = rec + "CND_path_best" + dnw;
            _rxn.Define(pCnd, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, cnd, "path", 
                              info::DetectorLayers().at(clas12::CND)));
            
            // Explicit, Regex-proof C++ Callable for CD_Path
            _rxn.Define(rec + "CD_Path", 
                [](const rad::RVecResultType& pCtof, const rad::RVecResultType& pCnd) {
                    rad::RVecResultType res(pCtof.size(), rad::consts::InvalidEntry<rad::ResultType_t>());
                    for (size_t i = 0; i < pCtof.size(); ++i) {
                        if (!rad::consts::IsInvalidEntry(pCtof[i]))      res[i] = pCtof[i];
                        else if (!rad::consts::IsInvalidEntry(pCnd[i]))  res[i] = pCnd[i];
                    }
                    return res;
                }, {rec + ctof + "_path", pCnd}
            );
        }

        // --- FORWARD TAGGER (FT) ---
        if (!_rxn.ColumnExists(rec + "FT_DetEnergy")) {
            _rxn.Define(rec + "FT_DetEnergy", rec + ftcal + "_energy");
            _rxn.Define(rec + "FT_DeltaEnergy", rec + fthodo + "_energy");
            _rxn.Define(rec + "FT_Time", rec + ftcal + "_time");
            
            _rxn.Define(rec + "FT_Path", "sqrt(" + rec + ftcal + "_x*" + rec + ftcal + "_x + " + 
                                                   rec + ftcal + "_y*" + rec + ftcal + "_y + " + 
                                                   rec + ftcal + "_z*" + rec + ftcal + "_z)");
        }

        // Map regions directly to the combinatorial candidates
        const std::vector<std::string> regions = {"FD", "CD", "FT"};
        const std::vector<std::string> vars = {"DetEnergy", "DeltaEnergy", "Time", "Path", "Sector"};
        for (const auto& reg : regions) {
            for (const auto& v : vars) AutoMapSimple(reg, v); 
        }
    }
/* template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildRegions() { */
/* std::string rec = rad::consts::data_type::Rec(); */
/* std::string dnw = rad::DoNotWriteTag(); */

/* // Dynamically resolve names */
/* std::string ecal = _detectorNames.DetName(clas12::ECAL); */
/* std::string ftof = _detectorNames.DetName(clas12::FTOF); */
/* std::string cnd  = _detectorNames.DetName(clas12::CND); */
/* std::string ctof = _detectorNames.DetName(clas12::CTOF); */
/* std::string ftcal = _detectorNames.DetName(clas12::FTCAL); */
/* std::string fthodo = _detectorNames.DetName(clas12::FTHODO); */
/* std::string trkName = _detectorNames.DetName(clas12::DC); */
/* std::string band = _detectorNames.DetName(clas12::BAND); */

/* // --- FORWARD DETECTOR (FD) --- */
/* if (!_rxn.ColumnExists(rec + "FD_DetEnergy")) { */
            
/* // Energy Sums */
/* _rxn.Define(rec + "FD_DetEnergy",  */
/* util::BuildLayerFunctionString("rad::clas12::util::SumValid", rec, ecal, "energy",  */
/* info::DetectorLayers().at(clas12::ECAL))); // Automatically uses PCAL, ECIN, ECOUT */

/* _rxn.Define(rec + "FD_DeltaEnergy",  */
/* util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ftof, "energy",  */
/* {clas12::FTOF1B, clas12::FTOF1A, clas12::FTOF2})); // Custom order */
            
/* // Timing (Fallback arrays) */
/* std::string tFtof = rec + "FTOF_time_best" + dnw; */
/* std::string tEcal = rec + "ECAL_time_best" + dnw; */

/* _rxn.Define(tFtof, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ftof, "time",  */
/* {clas12::FTOF1B, clas12::FTOF1A, clas12::FTOF2})); */
/* _rxn.Define(tEcal, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ecal, "time",  */
/* {clas12::PCAL, clas12::ECOUT, clas12::ECIN})); */

/* // FIX: Mirror clas12root exactly - Charge != 0 uses FTOF, else ECAL */
/* _rxn.Define(rec + "FD_Time",  */
/* [](const rad::RVecResultType& tof, const rad::RVecResultType& cal, const rad::Indices_t& pids) { */
/* rad::RVecResultType res(tof.size(), 0.0); */
/* for (size_t i = 0; i < tof.size(); ++i) { */
/* // Check if particle is neutral (Gamma, Neutron, Pi0, Unknown) */
/* bool isNeutral = (pids[i] == 22 || pids[i] == 2112 || pids[i] == 111 || pids[i] == 0); */
/* if (!isNeutral && tof[i] > 0.0) res[i] = tof[i]; */
/* else                            res[i] = cal[i]; */
/* } */
/* return res; */
/* }, {tFtof, tEcal, rec + "pid"} */
/* ); */

/* // Path  */
/* std::string pFtof = rec + "FTOF_path_best" + dnw; */
/* std::string pEcal = rec + "ECAL_path_best" + dnw; */

/* _rxn.Define(pFtof, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ftof, "path",  */
/* {clas12::FTOF1B, clas12::FTOF1A, clas12::FTOF2})); */
/* _rxn.Define(pEcal, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, ecal, "path",  */
/* {clas12::PCAL, clas12::ECOUT, clas12::ECIN})); */

/* // FIX: Mirror clas12root exactly - Charge != 0 uses FTOF, else ECAL */
/* _rxn.Define(rec + "FD_Path",  */
/* [](const rad::RVecResultType& pTof, const rad::RVecResultType& pCal, const rad::Indices_t& pids) { */
/* rad::RVecResultType res(pTof.size(), 0.0); */
/* for (size_t i = 0; i < pTof.size(); ++i) { */
/* bool isNeutral = (pids[i] == 22 || pids[i] == 2112 || pids[i] == 111 || pids[i] == 0); */
/* if (!isNeutral && pTof[i] > 0.0) res[i] = pTof[i]; */
/* else                             res[i] = pCal[i]; */
/* } */
/* return res; */
/* }, {pFtof, pEcal, rec + "pid"} */
/* ); */

/* // Sector */
/* _rxn.Define(rec + "FD_Sector", util::BuildFunctionString("rad::clas12::util::Fallback",  */
/* {rec + trkName + "_sector", util::ColName(rec, ftof, "sector", clas12::FTOF1B), util::ColName(rec, ecal, "sector", clas12::PCAL)})); */
/* } */

/* // --- CENTRAL DETECTOR (CD) --- */
/* if (!_rxn.ColumnExists(rec + "CD_DetEnergy")) { */
            
/* _rxn.Define(rec + "CD_DetEnergy",  */
/* util::BuildLayerFunctionString("rad::clas12::util::SumValid", rec, cnd, "energy",  */
/* info::DetectorLayers().at(clas12::CND))); */

/* _rxn.Define(rec + "CD_DeltaEnergy", rec + ctof + "_energy"); */
            
/* // Timing */
/* std::string tCnd = rec + "CND_time_best" + dnw; */
/* _rxn.Define(tCnd, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, cnd, "time",  */
/* info::DetectorLayers().at(clas12::CND))); */
            
/* /\* _rxn.Define(rec + "CD_Time", util::BuildFunctionString("rad::clas12::util::Fallback",  *\/ */
/* /\* {rec + ctof + "_time", tCnd})); *\/ */
/* // FIX: Debug-enabled CD_Time evaluator */
/* _rxn.Define(rec + "CD_Time",  */
/* [](const rad::RVecResultType& tCtof, const rad::RVecResultType& tCnd, const rad::Indices_t& pids) { */
/* rad::RVecResultType res(tCtof.size(), 0.0); */
                    
/* // --- SAFE DEBUG OUTPUT --- */
/* static std::atomic<int> debug_calls{0}; */
/* bool do_debug = false; */
/* if (debug_calls < 20) {  */
/* do_debug = true; */
/* debug_calls++; */
/* std::cout << "\n[DEBUG] CD_Time Evaluator | N_Tracks: " << tCtof.size() << "\n"; */
/* } */

/* for (size_t i = 0; i < tCtof.size(); ++i) { */
/* // CD Logic: Prefer CTOF, fallback to CND */
/* if (tCtof[i] != 0.0) res[i] = tCtof[i]; */
/* else                 res[i] = tCnd[i]; */
                        
/* // Print the exact state of the arrays for this track */
/* if (do_debug) { */
/* std::cout << "  -> Track " << i << ": PID=" << pids[i]  */
/* << " | tCtof=" << tCtof[i]  */
/* << " | tCnd=" << tCnd[i]  */
/* << " | RESULT=" << res[i] << "\n"; */
/* } */
/* } */
/* return res; */
/* }, {rec + ctof + "_time", tCnd, rec + "pid"} */
/* ); */
/* // Path */
/* std::string pCnd = rec + "CND_path_best" + dnw; */
/* _rxn.Define(pCnd, util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, cnd, "path",  */
/* info::DetectorLayers().at(clas12::CND))); */
            
/* _rxn.Define(rec + "CD_Path", util::BuildFunctionString("rad::clas12::util::Fallback",  */
/* {rec + ctof + "_path", pCnd})); */
/* } */

/* // --- FORWARD TAGGER (FT) --- */
/* if (!_rxn.ColumnExists(rec + "FT_DetEnergy")) { */
/* _rxn.Define(rec + "FT_DetEnergy", rec + ftcal + "_energy"); */
/* _rxn.Define(rec + "FT_DeltaEnergy", rec + fthodo + "_energy"); */
/* _rxn.Define(rec + "FT_Time", rec + ftcal + "_time"); */
            
/* _rxn.Define(rec + "FT_Path", "sqrt(" + rec + ftcal + "_x*" + rec + ftcal + "_x + " +  */
/* rec + ftcal + "_y*" + rec + ftcal + "_y + " +  */
/* rec + ftcal + "_z*" + rec + ftcal + "_z)"); */
/* } */

/* // --- BACKWARD ANGLE NEUTRON DETECTOR (BAND) --- */
/* if (info::DetectorLayers().count(clas12::BAND) && !_rxn.ColumnExists(rec + "BAND_DetEnergy")) { */
            
/* // Assuming layers 1-5 are TOF, and layer 6 is VETO (Standard clas12 BAND configuration) */
/* std::vector<int> tofLayers = {1, 2, 3, 4, 5}; */
            
/* _rxn.Define(rec + "BAND_DetEnergy",  */
/* util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, band, "energy", tofLayers)); */
                
/* _rxn.Define(rec + "BAND_DeltaEnergy", util::ColName(rec, band, "energy", 6));  */

/* _rxn.Define(rec + "BAND_Time",  */
/* util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, band, "time", tofLayers)); */
                
/* _rxn.Define(rec + "BAND_Path",  */
/* util::BuildLayerFunctionString("rad::clas12::util::Fallback", rec, band, "path", tofLayers)); */
/* } */

/* // Map regions directly to the combinatorial candidates */
/* const std::vector<std::string> regions = {"FD", "CD", "FT", "BAND"}; */
/* const std::vector<std::string> vars = {"DetEnergy", "DeltaEnergy", "Time", "Path", "Sector"}; */
/* for (const auto& reg : regions) { */
/* for (const auto& v : vars) AutoMapSimple(reg, v);  */
/* } */
/* } */

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildEventBanks() {
        if(!_rxn.ColumnExists("rec_event_category")) {
            _rxn.SetBranchAlias("REC_Event_category", "rec_event_category");
            _rxn.SetBranchAlias("REC_Event_topology", "rec_event_topology");
            _rxn.SetBranchAlias("REC_Event_beamCharge", "rec_event_beamCharge");
            _rxn.SetBranchAlias("REC_Event_helicity", "rec_event_helicity");
            _rxn.SetBranchAlias("REC_Event_helicityRaw", "rec_event_helicityRaw");
            _rxn.SetBranchAlias("REC_Event_startTime", "rec_event_startTime");
            _rxn.SetBranchAlias("REC_Event_RFTime", "rec_event_RFTime");
            _rxn.SetBranchAlias("REC_Event_procTime", "rec_event_procTime");
            _rxn.SetBranchAlias("REC_Event_liveTime", "rec_event_liveTime");

            if (_rxn.ColumnExists("RECFT_Event_startTime")) {
                _rxn.SetBranchAlias("RECFT_Event_startTime", "rec_ftbevent_startTime");
                _rxn.Define("rec_event_bestStartTime", "rec_ftbevent_startTime != 0 ? rec_ftbevent_startTime : rec_event_startTime");
            } else {
                _rxn.SetBranchAlias("REC_Event_startTime", "rec_event_bestStartTime");
            }

            _rxn.SetBranchAlias("RUN_config_run", "rec_run_num");
            _rxn.SetBranchAlias("RUN_config_event", "rec_event_num");
            _rxn.SetBranchAlias("RUN_config_unixtime", "rec_run_unixtime");
            _rxn.SetBranchAlias("RUN_config_trigger", "rec_run_trigger");
            _rxn.SetBranchAlias("RUN_config_timestamp", "rec_run_timestamp");
            _rxn.SetBranchAlias("RUN_config_type", "rec_run_type");
            _rxn.SetBranchAlias("RUN_config_mode", "rec_run_mode");
            _rxn.SetBranchAlias("RUN_config_torus", "rec_run_torus");
            _rxn.SetBranchAlias("RUN_config_solenoid", "rec_run_solenoid");
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildCalorimeter(int subdet, int layer) {
        const std::vector<std::string> cal_vars = {
            "time", "energy", "path", "chi2", "x", "y", "z", "hx", "hy", "hz",
            "lu", "lv", "lw", "du", "dv", "dw", "m2u", "m2v", "m2w", "m3u", "m3v", "m3w", "status", "sector"
        };
        const std::vector<std::string> cal_extras = {
            "dbstU", "dbstV", "dbstW", "rawEU", "rawEV", "rawEW", "recEU", "recEV", "recEW",
            "recDTU", "recDTV", "recDTW", "recFTU", "recFTV", "recFTW"
        };

        for (const auto& var : cal_vars) {
            SafeDefineAssoc(bank::Calorimeter(), var, subdet, layer);
            AutoMap(subdet, var, layer);
        }
        for (const auto& var : cal_extras) {
            SafeDefineAssoc(bank::Calorimeter(), var, subdet, layer, bank::CalExtras());
            AutoMap(subdet, var, layer);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildScintillator(int subdet, int layer) {
        const std::vector<std::string> sci_vars = {
            "time", "energy", "path", "chi2", "x", "y", "z", "hx", "hy", "hz", "sector", "status", "component"
        };
        const std::vector<std::string> sci_extras = {"dedx", "size", "layermulti"};

        for (const auto& var : sci_vars) {
            SafeDefineAssoc(bank::Scintillator(), var, subdet, layer);
            AutoMap(subdet, var, layer);
        }
        for (const auto& var : sci_extras) {
            SafeDefineAssoc(bank::Scintillator(), var, subdet, layer, bank::ScintExtras());
            AutoMap(subdet, var, layer);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildTracker(int subdet) {
        const std::vector<std::string> trk_vars = {"NDF", "sector", "status", "q", "chi2"};
        for (const auto& var : trk_vars) {
            SafeDefineAssoc(bank::Track(), var, subdet);
            AutoMap(subdet, var, -1);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildTrajectories(int subdet, int layer) {
        const std::vector<std::string> traj_vars = {"cx", "cy", "cz", "x", "y", "z", "edge", "path"};
        for (const auto& var : traj_vars) {
            SafeDefineAssoc(bank::Traj(), var, subdet, layer);
            AutoMap(subdet, var, layer);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildCherenkov(int subdet) {
        const std::vector<std::string> cher_vars = {"nphe", "time", "path", "sector", "chi2", "x", "y", "z", "dtheta", "dphi", "status"};
        for (const auto& var : cher_vars) {
            SafeDefineAssoc(bank::Cherenkov(), var, subdet);
            AutoMap(subdet, var, -1);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildForwardTagger(int subdet, int layer) {
        const std::vector<std::string> ft_vars = {"time", "energy", "path", "status", "x", "y", "z", "dx", "dy", "radius", "size", "chi2"};
        for (const auto& var : ft_vars) {
            SafeDefineAssoc(bank::ForwardTagger(), var, subdet, layer);
            AutoMap(subdet, var, layer);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::BuildCovMatrix() {
        const std::vector<std::string> cov_vars = {
            "C11", "C12", "C13", "C14", "C15", "C22", "C23", "C24", "C25",
            "C33", "C34", "C35", "C44", "C45", "C55"
        };
        for (const auto& var : cov_vars) {
            SafeDefineSimple(bank::CovMat(), var);
            AutoMapSimple(bank::CovMat(), var);
        }
    }

    // --- Private Synchronization Helpers ---

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::SafeDefineAssoc(const std::string& det, const std::string& item, int subdet, int layer, const std::string& val_bank) {
        std::string rawCol = "REC_" + (val_bank.empty() ? det : val_bank) + "_" + item;
        if (!_rxn.ColumnExists(rawCol)) return;
        
        std::string auxBaseName = _detectorNames.DetName(subdet) + "_" + item;
        if (layer >= 0) auxBaseName += "_L" + std::to_string(layer);
        
        if (!_rxn.ColumnExists(rad::consts::data_type::Rec() + auxBaseName)) {
            if(val_bank.empty()) _rxn.DefineDetectorAssociation(det, item, subdet, layer);
            else                 _rxn.DefineDetectorAssociation(det, item, subdet, layer, val_bank);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::SafeDefineSimple(const std::string& det, const std::string& item) {
        std::string rawCol = "REC_" + det + "_" + item;
        if (!_rxn.ColumnExists(rawCol)) return;

        if (!_rxn.ColumnExists(rad::consts::data_type::Rec() + det + "_" + item)) {
            _rxn.DefineSimpleAssociation(det, item);
        }
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::AutoMap(int subdet, const std::string& item, int layer) {
        std::string auxBaseName = _detectorNames.DetName(subdet) + "_" + item;
        if (layer >= 0) auxBaseName += "_L" + std::to_string(layer);
        AutoMapSimple(auxBaseName, "");
    }

    template<typename DS_t>
    inline void CLAS12DetectorBuilder<DS_t>::AutoMapSimple(const std::string& prefix, const std::string& item) {
        std::string auxBaseName = item.empty() ? prefix : prefix + "_" + item;
        std::string recType = rad::consts::data_type::Rec();
        std::string baseCol = recType + auxBaseName; 

        if (!_rxn.ColumnExists(baseCol)) return;

        for (const auto& pName : _rxn.ParticleNames()) {
            std::string candCol = recType + pName; 
            std::string outCol  = recType + pName + "_" + auxBaseName; 

            if (_rxn.ColumnExists(candCol) && !_rxn.ColumnExists(outCol)) {
                _rxn.Define(outCol, Form("ROOT::VecOps::Take(%s, %s)", baseCol.c_str(), candCol.c_str()));
            }
        }
    }

} // namespace clas12
} // namespace rad
