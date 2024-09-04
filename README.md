# Welcome to the SW4 Project on Southern Vancouver Island

The purpose of this project is to use SW4 to simulate seismic waves on Southern Vancouver Island. Our main objective is to produce different kinds of waveforms from seismic events with different parameters. We methodically create simulations that closely resemble actual seismic events, taking inputs from the Tim's catalog to improve our comprehension of seismic behavior in this area.

Here's the original repository of sw4: https://github.com/geodynamics/sw4.

## List of Files and Their Purposes:

- **ProblemSetup_SVI.ipynb**: This notebook contains the setup and configuration for SW4 simulations specific to Southern Vancouver Island. It sets up event and station locations, processes velocity model data, visualizes model boundaries, computes grid spacing, and creates input files for SW4 simulations, among other seismic wave simulation-related duties.

- **Savard_VpVs.txt**: This file is the velocity model relevant to the region with P- and S-waves velocities.

- **VCI.topo**: This file represents the topographic data for the Vancouver Island area.

- **VCI_DC.in**: This file contain input parameters or configuration settings for a specific aspect of the SW4 simulations.

- **VancouverIsland.ppmod**: This file is pfile generated for a problem which defines the materials in all point.

- **map_test.py**: A script used for generating maps related to the seismic data.

- **plot_sw4img.m**: A script used for plotting SW4 simulation outputs (special extension .sw4img).

- **results_plot.py**: A script used for plotting some results from the SW4 simulations.

- **stations.csv**: This CSV file contains data regarding the location and attributes of stations used in the seismic data analysis.

- **submit_test.txt**: This file contain instructions for submitting SW4 simulation jobs on Talapas (must remove the .txt extension).

- **sw4FileFunctions.py**: This script contains custom functions elated to file handling within the SW4 project.

### NonLinLoc Folder:

A new folder related to the NonLinLoc software has been added to the repository as of 09/04/24. The scripts within this folder are designed to assist in locating selected events from a catalog. Specifically, the folder includes:

-   An example input file compatible with the NLLoc command for event location.
-   A script to stack events, or not, and apply the P and S wave picker. The resulting seismic phase arrivals can then be formatted using another script for NLLoc.
-   A work-*in-progress* notebook for generating the necessary grid files for event location, utilizing Vel2Grid and Grid2Time.
