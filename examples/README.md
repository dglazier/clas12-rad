## `clas12-rad` Analysis Scripts

This directory contains a suite of core analysis macros designed for reaction aware dataframes[cite: 7]. They demonstrate how to process HIPO data efficiently using declarative recipes, combinatorics, and lazy evaluation.

### 1. Script Overviews

*   **`MinimalExample.C`**: A foundational "Hello World" script. It demonstrates the strict order of operations for initialization, particle candidate definitions, combinatorial generation, and basic missing mass calculations.
*   **`PlotDetectorData.C`**: Highlights the `CLAS12DetectorBuilder` class[cite: 2]. It automatically projects relational detector banks (like Forward Detector timing and Calorimeter energy) onto combinatorial tracks with zero relational overhead[cite: 2]. *Note: All relations are specified in the github repo at `DetectorInfoMap.md`*[cite: 2].
*   **`ResolutionAnalysis.C`**: Showcases multi-stream processing by running both Reconstructed and Monte Carlo Truth streams simultaneously. It utilizes truth-matching and cross-stream differences to calculate and plot resolutions.
*   **`FilterAnalysis.C`**: Demonstrates physics selection and lazy masking. It shows how to apply minimum momentum thresholds and exclusivity cuts (like on Missing Mass Squared) before plotting and safely snapshotting the surviving signal to a flat TTree without destroying the underlying data arrays.
*   **`HistogramAnalysis.C`**: Acts as a high-speed histogramming engine that bypasses the snapshot step entirely. It utilizes lazy masking and split axes (like discrete Q2 bins) to strictly aggregate histograms in RAM.

### 2. How to Run

These scripts are designed to be executed directly from the command line using ROOT. You can pass single HIPO files or use wildcards to process entire datasets at once. 

To run a macro with a wildcard file path, wrap the macro execution in single quotes and the file path in double quotes. This ensures the exact string is passed directly into the C++ function and prevents your shell from prematurely expanding the wildcard.

```bash
root -l -b -q 'HistogramAnalysis.C("~/Jlab/clas12/data/simulation/RhoFeb24/rho-7221-9*.hipo")'
```

**Command Breakdown:**
*   **`-l`**: Skips the ROOT splash screen.
*   **`-b`**: Runs in batch mode, suppressing graphical windows to save time and prevent X11 forwarding issues on remote servers.
*   **`-q`**: Quits ROOT automatically once the macro finishes.