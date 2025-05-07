CTSM53 LUH3 Land Data Prep Workflow

This is step two in the bigger picture of generating CTSM raw PFT/CFT data for mksurfdata_esmf. 

The prior stage is to take the 1km MODIS, AVHRR, ICESAT Satellite and EarthSTAT 
cropping products processed into a consistent format including grids and land masks.
This step generates the LUH3 time series and LUH3 Descriptor files that are combined 
in the CTSM53 Land Data Tool. 

This instance is configured to produce the CMIP7 CTSM53 data.

Tools require python configured with the NCAR python environment.

On CGD machines
module load lang/python

On HPC machines
module load conda
conda activate npl

To configure the tools follow the instructions in the README.configure file

To process the CMIP7 CTSM53 LUH3 data follow the instructions in the README.process file

