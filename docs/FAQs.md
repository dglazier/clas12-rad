## Custom Track Filtering (Lambda Functions)

**How do I define custom track selections beyond just PID?**
Sometimes, identifying a particle requires checking more than just its standard PDG code. You might need to verify that the track was registered in a specific sub-detector, or that it passed a specific tracking quality cut. 

Instead of passing a simple PID filter, you can pass a custom C++ **Lambda function** directly into `SetParticleCandidates`[cite: 6]. The framework will evaluate this logic on the raw arrays before generating the combinatorial matrix[cite: 6].

**CLAS12 Example (HIPO & clas12defs)**
In CLAS12, you often need to separate tracks based on whether they hit the Forward Detector (FD) or the Central Detector (CD). 

To make this easy, the framework automatically unpacks the complex `REC::Particle` status word into a clean `rec_region` array during initialization. You can compare this array directly against the predefined geometry flags in `clas12defs.h` (e.g., `rad::clas12::FD` or `rad::clas12::CD`).

You can use these flags alongside a custom lambda to, for example, strictly identify a **positive pion** ($\pi^+$) in the Forward Detector. *(Note: You must explicitly check for `pid == 211` rather than using the absolute value, otherwise negative pions will incorrectly be added to your positive pion candidate list!)*

```cpp
#include "clas12defs.h"

// 1. Define the CLAS12 Forward Detector Lambda filter for an explicit pi+ (211)
auto FD_Pip_Filter = [](const ROOT::RVecI& pid, const ROOT::RVecI& region) {
    // Return indices where the track is a pi+ AND its region matches the FD flag
    return ROOT::VecOps::Nonzero(
        (pid == 211) && 
        (region == rad::clas12::FD)
    );
};

// 2. Inject it into the Combinatorial Engine
clas12_df.SetParticleCandidates("pip", Role_Pip, FD_Pip_Filter, {"rec_pid", "rec_region"});
```

*Note for CLAS12 Users:* The `clas12-rad` extension also supports string-based JIT-compilation via `SetParticleCandidatesExpr` for rapid prototyping[cite: 6]. You must use the exact column names as they appear in the TTree (e.g., `rec_pid` and `rec_region`):
```cpp
clas12_df.SetParticleCandidatesExpr("pip", 
    "rec_pid == 211 && rec_region == rad::clas12::FD");
```
*(While the string expression is slightly easier to read, the Lambda approach is type-safe and compiles much faster for production analyses.)*

## Inspecting Raw HIPO Banks

When developing an analysis with `clas12-rad`, you may need to verify exactly what data is stored in your input files before setting up complex detector cuts or relying on synthesized variables. The framework provides a direct diagnostic tool to peek into the raw HIPO structures.

By calling `clas12_df.InspectBanks({"MC_Lund", "REC_Particle"}, 3);`, you are instructing the framework to dump the contents of the `MC_Lund` and `REC_Particle` banks for exactly the first `3` events.

Here is an example of the output:

```text
==================================================================
=== Inspecting Banks (First 3 Events) ===

>>> DUMPING BANK: MC_Lund <<<
+-----+------------------+----------------+---------------+------------------+--------------+
| Row | MC_Lund_daughter | MC_Lund_energy | MC_Lund_index | MC_Lund_lifetime | MC_Lund_mass | 
+-----+------------------+----------------+---------------+------------------+--------------+
| 0   | 0                | 3.558740       | 1             | 0.000000         | 0.139570     | 
|     | 0                | 1.176800       | 2             | 0.000000         | 0.139570     | 
|     | 0                | 1.005000       | 3             | 0.000000         | 0.938272     | 
|     | 0                | 5.397730       | 4             | 0.000000         | 0.000511     | 
+-----+------------------+----------------+---------------+------------------+--------------+
...
+-----+----------------+-------------+------------+------------+------------+
| Row | MC_Lund_parent | MC_Lund_pid | MC_Lund_px | MC_Lund_py | MC_Lund_pz | 
+-----+----------------+-------------+------------+------------+------------+
| 0   | 0              | 211         | 0.083195   | 0.638408   | 3.497240   | 
|     | 0              | -211        | 0.264302   | -0.267596  | 1.106310   | 
|     | 0              | 2212        | -0.173562  | -0.241458  | 0.203075   | 
|     | 0              | 11          | -0.173934  | -0.129354  | 5.393380   | 
+-----+----------------+-------------+------------+------------+------------+

>>> DUMPING BANK: REC_Particle <<<

+-----+-------------------+---------------------+----------------------+------------------+-----+
| Row | REC_Particle_beta | REC_Particle_charge | REC_Particle_chi2pid | REC_Particle_pid | ... | 
+-----+-------------------+---------------------+----------------------+------------------+-----+
| 0   | 0.992860          | -1                  | -1.292901            | -211             | ... | 
|     | -99.000000        | 1                   | 9999.000000          | 0                | ... | 
|     | -99.000000        | 1                   | 9999.000000          | 0                | ... | 
|     | 0.971304          | 0                   | 9999.000000          | 22               | ... | 
|     | 0.860554          | 0                   | 9999.000000          | 2112             | ... | 
+-----+-------------------+---------------------+----------------------+------------------+-----+
...
==================================================================
```

> **How to use it:**
> *   **The Event Limit:** The second argument (in this case, `3`) specifies exactly how many events to process and print. This prevents your terminal from being flooded while giving you a sufficient sample size to inspect.
> *   **Verifying Bank Existence:** Not all HIPO files contain the same data (e.g., Monte Carlo files contain truth banks, while pure detector data will not). This diagnostic allows you to check if the specific banks you are requesting actually exist in your input files before the framework tries to parse them.
> *   **Checking Typical Values:** It allows you to see the exact structure and typical range of values stored in the raw HIPO relational database. For example, you can visually verify standard conventions, like `REC_Particle_pid` returning `-211` for a $\pi^-$, or check if variables like `REC_Particle_beta` and `REC_Particle_chi2pid` are defaulting to flag values (e.g., `-99.000000` or `9999.000000`) for neutral tracks.