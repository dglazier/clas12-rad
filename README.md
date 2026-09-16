# CLAS12-RAD: Reaction Analysis & Design for CLAS12

**A Declarative, Vectorized RDataFrame Framework for CLAS12 HIPO Analysis**

`clas12-rad` extends the core [RAD framework](https://github.com/dglazier/rad) to natively process Jefferson Lab CLAS12 HIPO data. By coupling ROOT's `RDataFrame` with custom HIPO datasources (`RHipoDS`, `RIguanaDS`), it provides high-performance, multi-threaded combinatorics, detector synthesis, and physics analysis.

The framework is **header-only** and runs directly within ROOT script execution—**no building, compilation steps, or Makefiles are required**.

---

## 📚 Documentation & Reference

* **Core RAD Repository:** [github.com/dglazier/rad](https://github.com/dglazier/rad)
  * Consult the general core documentation in `docs/`, particularly [`docs/FAQs.md`](https://github.com/dglazier/rad/blob/master/docs/FAQs.md), for foundational concepts regarding combinatorics, memory safety, and processor recipes.
* **CLAS12 Extension Repository:** [github.com/dglazier/clas12-rad](https://github.com/dglazier/clas12-rad)
  * [`docs/DetectorInfoMap.md`](https://github.com/dglazier/clas12-rad/blob/master/docs/DetectorInfoMap.md): Complete mapping reference of all synthesized CLAS12 detector layers, compound regions (FD, CD, FT), and reconstructed hit variables.
  * [`docs/FAQs.md`](https://github.com/dglazier/clas12-rad/blob/master/docs/FAQs.md): Practical recipes and frequently asked questions developed specifically for CLAS12 analysis patterns (such as status/region filtering and detector pass-throughs).

---

## ⚙️ Environment Setup & Installation

### Running on JLab ifarm
On the JLab `ifarm`, the required `hipo`, `qadb`, and `iguana` libraries, headers, and dependencies are already available. You only need to load the CLAS12 environment module, unload `clas12root`, and set your repository paths:

```csh
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12
module unload clas12root

# Clone and configure paths
git clone --recurse-submodules https://github.com/dglazier/clas12-rad.git
setenv CLAS12RAD /path/to/clas12-rad
setenv RAD ${CLAS12RAD}/rad
setenv ROOT_INCLUDE_PATH ${ROOT_INCLUDE_PATH}:${RAD}/include:${CLAS12RAD}/include
```

### Local Installation (Laptops & Off-Site Clusters)
`clas12-rad` is header-only, but requires access to the base `rad` framework. Clone with submodules to download `rad` automatically:

```bash
git clone --recurse-submodules https://github.com/dglazier/clas12-rad.git
```

Set the corresponding paths in your shell environment:

```csh
setenv CLAS12RAD /path/to/clas12-rad
setenv RAD ${CLAS12RAD}/rad
setenv ROOT_INCLUDE_PATH ${ROOT_INCLUDE_PATH}:${RAD}/include:${CLAS12RAD}/include
```

*(If you already have a standalone clone of `rad`, simply set `RAD` to your existing location instead of cloning the submodule).*

---

## 🧪 Analysis Examples Walkthrough

The `examples/` directory provides a progressive set of macros illustrating framework capabilities from basic topologies to production analysis:

* **`MinimalExample.C` (Hello World):**
  Introduces the core execution pipeline: initializing the `AnalysisManager`, configuring `CLAS12Reaction<RHipoDS>`, defining candidate PIDs from `REC::Particle`, generating combinatorics, computing a 4-vector difference (`miss`), and booking missing mass histograms.

* **`HistogramAnalysis.C` (High-Speed Histogramming & Splitting):**
  Demonstrates pure in-memory histogram accumulation without flat-tree snapshotting overhead. Uses `AddSplit` to automatically project 1D and 2D physics histograms across multiple kinematic bins (e.g., slicing distributions by discrete $Q^2$ ranges).

* **`FilterAnalysis.C` (Physics Selection & Lazy Masking):**
  Illustrates non-destructive event selection via `PhysicsSelection`. Demonstrates applying track-level minimum momentum cuts and reaction-level exclusivity cuts (e.g., Missing Mass Squared) while preserving underlying array structures via boolean masks.

* **`ResolutionAnalysis.C` (Dual-Stream MC Truth Matching):**
  Runs both Reconstructed (`Rec`) and Truth (`Truth`) streams concurrently. Leverages `SetupMatching()` to link `REC::Particle` tracks to generated particles in `MC::Lund` via `MC::GenMatch`, and uses `CrossStreamDifferences` to calculate track resolutions ($\Delta P$, $\Delta\theta$, $\Delta\phi$).

* **`PlotDetectorData.C` (Zero-Overhead Detector Synthesis):**
  Showcases the `CLAS12DetectorBuilder` class. Automatically traverses relational HIPO banks (`REC::Calorimeter`, `REC::Scintillator`, `REC::Traj`, `REC::Track`) and projects high-level region variables (FD, CD, FT) directly onto combinatorial candidates. Full variable mappings are cataloged in `DetectorInfoMap.md`.

* **`HistogramWithQadb.C` (Quality Assurance Database Interception):**
  Integrates the `QADBFilter` directly into the `RHipoDS` data source hook. Skips defect bins (e.g., `TotalOutlier`, `Misc`) early during file reading to avoid unneeded bank deserialization, and tallies valid analyzed Faraday Cup charge in nC.

* **`ProcessIguana.C` (In-RAM Pre-Processing with RIguanaDS):**
  Employs `RIguanaDS` to intercept raw HIPO banks before `RDataFrame` processing. Executes multi-threaded IGUANA algorithms on the fly (vertex filters, sector finders, fiducial cuts, momentum corrections, and inclusive kinematics) and manages non-destructive masking.

* **`ProcessDet_eppippim.C` (Full Production Benchmark):**
  A comprehensive, realistic analysis macro for exclusive $\rho$ electroproduction ($e p \to e' p' \pi^+ \pi^-$). Integrates FTB ambiguity resolution (`UseFTB()`), MC truth matching roles, region-specific lambda filters (`FD_Pip_Filter`), full detector building, stream snapshots, and performance diagnostics.

---

## 🚀 Running the Examples

Execute any example directly using ROOT. Wrap the macro and input file argument in quotes to prevent premature shell expansion of wildcards:

```bash
root -l -b -q 'MinimalExample.C("path/to/data_*.hipo")'
root -l -b -q 'HistogramAnalysis.C("path/to/data_*.hipo")'
root -l -b -q 'ProcessDet_eppippim.C'
```