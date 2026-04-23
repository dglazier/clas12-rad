/**
 * @file CLAS12Utilities.h
 * @brief High-performance vector math and utility functions for CLAS12.
 * @details 
 * Provides SIMD-compatible RDataFrame vector manipulations, string builders 
 * for JIT compilation, and detector mapping logic. Strictly header-only and 
 * interface/implementation separated.
 */

#pragma once

#include "Constants.h"

#include <ROOT/RVec.hxx>
#include <string>
#include <vector>

namespace rad {
namespace clas12 {
namespace util {

    // =================================================================================
    // Type Definitions
    // =================================================================================
    
    /** @brief Matrix storing associated detector hit indices: [track_idx][list_of_hits] */
    using detector_matrix_t = ROOT::VecOps::RVec<ROOT::VecOps::RVec<short>>;

    // =================================================================================
    // INTERFACE DECLARATIONS
    // =================================================================================

    // --- Indexing & Matrix Flattening ---

    /**
     * @brief Creates a reverse lookup matrix for detector hits.
     * @details Maps HIPO detector banks (where index = hit, value = track) 
     * to a track-aligned matrix (where index = track, value = list of hits).
     */
    template<typename T, typename Tn>
    detector_matrix_t ReverseIndexN(const ROOT::VecOps::RVec<short>& vec, Tn nentries);

    /**
     * @brief Projects simple detector indices into a flat array aligned with rec_px.
     * @details Used for banks with no subdetector ID (e.g., CovMatrix).
     */
    template<typename T>
    ROOT::RVec<T> FlattenDetectorInfo(const detector_matrix_t& indices, const ROOT::RVec<T>& vals);

    /**
     * @brief Projects nested Sub-Detector indices into a flat array (NO LAYER CHECK).
     * @details Used for banks like Cherenkov that do not have a layer column.
     */
    template<typename T, typename TDet>
    ROOT::RVec<T> FlattenSubDetectorInfo(const detector_matrix_t& indices, const ROOT::RVec<T>& vals, 
                                         const ROOT::RVec<TDet>& det_id_vec, int target_det);

    /**
     * @brief Projects nested Sub-Detector indices into a flat array (WITH LAYER CHECK).
     * @details Used for banks like Calorimeter that contain multiple layers.
     */
    template<typename T, typename TDet, typename TLayer>
    ROOT::RVec<T> FlattenSubDetectorInfo(const detector_matrix_t& indices, const ROOT::RVec<T>& vals, 
                                         const ROOT::RVec<TDet>& det_id_vec, const ROOT::RVec<TLayer>& layer_id_vec, 
                                         int target_det, int target_layer);

    // --- Physics Region Vector Math ---

    /**
     * @brief Merges REC and RECFT banks. Prefers Forward Tagger Based (FTB) value if valid.
     */
    template<typename T>
    ROOT::RVec<T> MergeFTB(const ROOT::RVec<T>& base, const ROOT::RVec<T>& ftb);

    /**
     * @brief Sums valid energies across ANY number of detector layers.
     * @details Uses C++17 variadic templates for high-speed SIMD folding.
     */
    template <typename... Vecs>
    rad::RVecResultType SumValid(const Vecs&... vecs);
    /**
     * @brief Resolves ambiguity by preferring Layer A > Layer B > Layer C...
     * @details Uses C++17 variadic templates. Takes ANY number of layers in order of preference.
     */
    // template <typename T, typename... Vecs>
    // ROOT::RVec<T> Fallback(const ROOT::RVec<T>& first, const Vecs&... rest);
    template <typename... Vecs>
    rad::RVecResultType Fallback(const Vecs&... vecs);
    /**
     * @brief Resolves array choice based on whether the particle track has charge.
     */
    template<typename TCharge, typename T>
    ROOT::RVec<T> SelectByCharge(const ROOT::RVec<TCharge>& charge, const ROOT::RVec<T>& chg_val, const ROOT::RVec<T>& neu_val);

    // --- Dynamic String Builders for JIT Compilation ---

    /** @brief Dynamically constructs a column name (e.g., "rec_ECAL_energy_L1") */
    std::string ColName(const std::string& prefix, const std::string& det, const std::string& item, int layer = -1);

    /** @brief Dynamically builds a generic function call string (e.g., "myFunc(col1, col2)") */
    std::string BuildFunctionString(const std::string& funcName, const std::vector<std::string>& cols);

    /** @brief Dynamically builds a layer-specific fallback/summation string based on a layer map. */
    std::string BuildLayerFunctionString(const std::string& funcName, const std::string& prefix, 
                                         const std::string& det, const std::string& item, 
                                         const std::vector<int>& layers);

} // namespace util
} // namespace clas12
} // namespace rad


// =================================================================================
// IMPLEMENTATIONS
// =================================================================================

namespace rad {
namespace clas12 {
namespace util {

    template<typename T, typename Tn>
    inline detector_matrix_t ReverseIndexN(const ROOT::VecOps::RVec<short>& vec, Tn nentries) {
        if (nentries == 0) return detector_matrix_t();
        if (nentries < vec.size()) nentries = vec.size();
        
        detector_matrix_t result(nentries);
        T entry = 0;
        for (auto idx : vec) {
            // Only add valid track indices
            if (idx >= 0 && idx < result.size()) {
                result[idx].push_back(entry);
            }
            ++entry;
        }
        return result;
    }
// 1. For detectors with layers (e.g. FTOF, ECAL)
  
#include <atomic>
#include <iostream>

     
// 1. For detectors with layers (e.g. FTOF, ECAL)
        template <typename T_val, typename T_det, typename T_lay>
        inline rad::RVecResultType FilterSubDetectorInfo(const ROOT::RVec<short>& pindex, const ROOT::RVec<T_val>& vals, const ROOT::RVec<T_det>& dets, const ROOT::RVec<T_lay>& layers, int subdet, int layer, size_t n) {
            // Initialize with standardized InvalidEntry (NaN) instead of 0.0
            rad::RVecResultType res(n, rad::consts::InvalidEntry<rad::ResultType_t>());
            
            for (size_t i = 0; i < pindex.size(); ++i) {
                if (pindex[i] >= 0 && dets[i] == subdet && (layer < 0 || layers[i] == layer)) {
                    size_t p_idx = (size_t)(pindex[i] + 2); // +2 for injected beam offset
                    
                    // Check for invalid entry rather than 0.0
                    if (p_idx < n && std::isnan(res[p_idx])) {
                        res[p_idx] = static_cast<rad::ResultType_t>(vals[i]);
                    }
                }
            }
            return res;
        }

        // 2. For detectors without layers (e.g. HTCC, FTCAL)
        template <typename T_val, typename T_det>
        inline rad::RVecResultType FilterDetectorInfo(const ROOT::RVec<short>& pindex, const ROOT::RVec<T_val>& vals, const ROOT::RVec<T_det>& dets, int subdet, size_t n) {
            rad::RVecResultType res(n, rad::consts::InvalidEntry<rad::ResultType_t>());
            
            for (size_t i = 0; i < pindex.size(); ++i) {
                if (pindex[i] >= 0 && dets[i] == subdet) {
                    size_t p_idx = (size_t)(pindex[i] + 2);
                    
                    if (p_idx < n && std::isnan(res[p_idx])) {
                        res[p_idx] = static_cast<rad::ResultType_t>(vals[i]);
                    }
                }
            }
            return res;
        }

        // 3. For simple flat arrays (e.g. CovMat, Traj)
        template <typename T_val>
        inline rad::RVecResultType FilterSimpleInfo(const ROOT::RVec<short>& pindex, const ROOT::RVec<T_val>& vals, size_t n) {
            rad::RVecResultType res(n, rad::consts::InvalidEntry<rad::ResultType_t>());
            
            for (size_t i = 0; i < pindex.size(); ++i) {
                if (pindex[i] >= 0) {
                    size_t p_idx = (size_t)(pindex[i] + 2);
                    
                    if (p_idx < n && std::isnan(res[p_idx])) {
                        res[p_idx] = static_cast<rad::ResultType_t>(vals[i]);
                    }
                }
            }
            return res;
        }
    // template<typename T, typename TDet>
    // inline ROOT::RVec<T> FlattenSubDetectorInfo(const detector_matrix_t& indices, const ROOT::RVec<T>& vals, 
    //                                             const ROOT::RVec<TDet>& det_id_vec, int target_det) {
    //     ROOT::RVec<T> out(indices.size(), rad::consts::InvalidEntry<T>());
    //     for (size_t i = 0; i < indices.size(); ++i) {
    //         if (indices[i].empty()) continue;
    //         for (auto hit_idx : indices[i]) {
    //             if (hit_idx >= 0 && hit_idx < vals.size() && hit_idx < det_id_vec.size()) {
    //                 if (det_id_vec[hit_idx] == target_det) {
    //                     out[i] = vals[hit_idx];
    //                     break; 
    //                 }
    //             }
    //         }
    //     }
    //     return out;
    // }

    // template<typename T, typename TDet, typename TLayer>
    // inline ROOT::RVec<T> FlattenSubDetectorInfo(const detector_matrix_t& indices, const ROOT::RVec<T>& vals, 
    //                                             const ROOT::RVec<TDet>& det_id_vec, const ROOT::RVec<TLayer>& layer_id_vec, 
    //                                             int target_det, int target_layer) {
    //     ROOT::RVec<T> out(indices.size(), rad::consts::InvalidEntry<T>());
    //     for (size_t i = 0; i < indices.size(); ++i) {
    //         if (indices[i].empty()) continue;
    //         for (auto hit_idx : indices[i]) {
    //             if (hit_idx >= 0 && hit_idx < vals.size() && hit_idx < det_id_vec.size() && hit_idx < layer_id_vec.size()) {
    //                 if (det_id_vec[hit_idx] == target_det) {
    //                     if (target_layer < 0 || layer_id_vec[hit_idx] == target_layer) {
    //                         out[i] = vals[hit_idx];
    //                         break; 
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     return out;
    // } 

  template<typename T>
    inline ROOT::RVec<T> MergeFTB(const ROOT::RVec<T>& base, const ROOT::RVec<T>& ftb) {
        ROOT::RVec<T> out = base;
        for (size_t i = 0; i < base.size(); ++i) {
            // FIX: Use IsInvalidEntry to safely handle NaNs
            if (i < ftb.size() && !rad::consts::IsInvalidEntry(ftb[i]) && ftb[i] != 0) {
                 out[i] = ftb[i];
            }
        }
        return out;
    }
  

  // =====================================================================
        // SAFE FALLBACK (Unified & RAD Standardized)
        // =====================================================================
 // =====================================================================
        // SAFE FALLBACK (Variadic Wrapper)
        // =====================================================================
 #include <atomic>
        #include <iostream>

        // =====================================================================
        // SAFE FALLBACK (Variadic Wrapper with DEBUG)
        // =====================================================================
        template <typename... Vecs>
        inline rad::RVecResultType Fallback(const Vecs&... vecs) {
            // Pack the variadic arguments into a vector internally
            std::vector<rad::RVecResultType> arrays = { vecs... };
            
            if (arrays.empty()) return rad::RVecResultType();
            size_t n = arrays.front().size();
            
            // Default everything to Invalid (NaN)
            rad::RVecResultType res(n, rad::consts::InvalidEntry<rad::ResultType_t>());
            
            // --- SAFE DEBUG OUTPUT ---
            // static std::atomic<int> debug_calls{0};
            // bool do_debug = false;
            // if (debug_calls < 20) { 
            //     do_debug = true;
            //     debug_calls++;
            //     std::cout << "\n[DEBUG] Variadic Fallback | N_Tracks: " << n << " | N_Arrays Mapped: " << arrays.size() << "\n";
            // }

            for (size_t i = 0; i < n; ++i) {
	      // if (do_debug) std::cout << "  -> Track " << i << ": ";
                
                bool found_valid = false;
                for (size_t j = 0; j < arrays.size(); ++j) {
                    const auto& arr = arrays[j];
                    
                    //if (do_debug) std::cout << "Arr[" << j << "]=" << arr[i] << "  ";
                    
                    // Assign the first valid hit we find
                    if (!found_valid && !rad::consts::IsInvalidEntry(arr[i])) {
                        res[i] = arr[i];
                        found_valid = true;
                        // In debug mode, we intentionally DO NOT break here 
                        // so we can see what the other fallback arrays contain!
			// if (!do_debug) break; 
                    }
                }
                //if (do_debug) std::cout << "| RESULT=" << res[i] << "\n";
            }
            return res;
        } 
        // =====================================================================
        // SAFE SUMVALID (Variadic Wrapper)
        // =====================================================================
        template <typename... Vecs>
        inline rad::RVecResultType SumValid(const Vecs&... vecs) {
            // Pack the variadic arguments into a vector internally!
            std::vector<rad::RVecResultType> arrays = { vecs... };
            
            if (arrays.empty()) return rad::RVecResultType();
            size_t n = arrays.front().size();
            
            // Sums start at 0.0 (a valid physics energy sum)
            rad::RVecResultType res(n, 0.0); 
            
            for (size_t i = 0; i < n; ++i) {
                for (const auto& arr : arrays) {
                    if (!rad::consts::IsInvalidEntry(arr[i])) {
                        res[i] += arr[i];
                    }
                }
            }
            return res;
        } 

    template<typename TCharge, typename T>
    inline ROOT::RVec<T> SelectByCharge(const ROOT::RVec<TCharge>& charge, const ROOT::RVec<T>& chg_val, const ROOT::RVec<T>& neu_val) {
        ROOT::RVec<T> res(charge.size(), rad::consts::InvalidEntry<T>());
        for(size_t i = 0; i < charge.size(); ++i) {
            if(charge[i] != 0 && i < chg_val.size() && chg_val[i] != rad::consts::InvalidEntry<T>()) res[i] = chg_val[i];
            else if (i < neu_val.size()) res[i] = neu_val[i];
        }
        return res;
    }

    // --- String Builders ---

    inline std::string ColName(const std::string& prefix, const std::string& det, const std::string& item, int layer) {
        std::string name = prefix + det + "_" + item;
        if (layer >= 0) name += "_L" + std::to_string(layer);
        return name;
    }

    inline std::string BuildFunctionString(const std::string& funcName, const std::vector<std::string>& cols) {
        std::string expr = funcName + "(";
        for (size_t i = 0; i < cols.size(); ++i) {
            expr += cols[i];
            if (i < cols.size() - 1) expr += ", ";
        }
        expr += ")";
        return expr;
    }

    inline std::string BuildLayerFunctionString(const std::string& funcName, const std::string& prefix, 
                                                const std::string& det, const std::string& item, 
                                                const std::vector<int>& layers) {
        std::vector<std::string> cols;
        for (int l : layers) cols.push_back(ColName(prefix, det, item, l));
        return BuildFunctionString(funcName, cols);
    }

} // namespace util
} // namespace clas12
} // namespace rad
