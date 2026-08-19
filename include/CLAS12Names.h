/**
 * @file clas12names.h
 * @brief Standardized bank names and detector layouts for CLAS12.
 * @details
 * Eliminates magic strings from the CLAS12 analysis framework.
 * Provides a central registry for detector-to-layer geometry mappings,
 * ensuring robust and dynamic iteration over detector components.
 */

#pragma once

#include "clas12defs.h"

#include <string>
#include <vector>
#include <map>

namespace rad {
namespace clas12 {

    // =================================================================================
    // Standardized Bank Names
    // =================================================================================
    namespace bank {
        inline const std::string Calorimeter()   { return "Calorimeter"; }
        inline const std::string Scintillator()  { return "Scintillator"; }
        inline const std::string Cherenkov()     { return "Cherenkov"; }
        inline const std::string ForwardTagger() { return "ForwardTagger"; }
        inline const std::string Track()         { return "Track"; }
        inline const std::string Traj()          { return "Traj"; }
        inline const std::string CovMat()        { return "CovMat"; }
        inline const std::string CalExtras()     { return "CalExtras"; }
        inline const std::string ScintExtras()   { return "ScintExtras"; }
    }

    // =================================================================================
    // Detector Geometry Registry
    // =================================================================================
    namespace info {
        /**
         * @brief Maps a primary Detector ID to its corresponding Layer IDs.
         * @details 
         * Provides the structural layout of the CLAS12 spectrometer. Used by the 
         * CLAS12DetectorBuilder to dynamically iterate over detector layers without 
         * hardcoding layer constants.
         * * @return const std::map<int, std::vector<int>>& Reference to the static layout map.
         */
        inline const std::map<int, std::vector<int>>& DetectorLayers() {
            // Static initialization ensures the map is only built once and safely shared
            static const std::map<int, std::vector<int>> _map = {
                {clas12::FTOF, {clas12::FTOF1A, clas12::FTOF1B, clas12::FTOF2}},
                {clas12::ECAL, {clas12::PCAL, clas12::ECIN, clas12::ECOUT}},
                {clas12::CND,  {clas12::CND1, clas12::CND2, clas12::CND3}},
                {clas12::DC,   {clas12::DC1, clas12::DC2, clas12::DC3, clas12::DC4, clas12::DC5, clas12::DC6}},
                {clas12::CVT,  {clas12::CVT1, clas12::CVT2, clas12::CVT3, clas12::CVT4, clas12::CVT5, 
                                clas12::CVT6, clas12::CVT7, clas12::CVT8, clas12::CVT9, clas12::CVT10, 
                                clas12::CVT11, clas12::CVT12}}
            };
            return _map;
        }
    }

} // namespace clas12
} // namespace rad
