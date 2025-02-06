# rrtmg-radiative-kernels

Repository that includes script for calculating RRTMG radiative kernels broadly following the metholodology introduced in Zelinka et al. (2012) and iterated on by Wall et al. (2023) and Duran et al. (2024).

The scripts are originally written in MATLAB and utilize the MATLAB/MEX interface with RRTMG. The goal is to transition these scripts into Python using a Python interface with RRTMG. Some of the Python code is not optimized as we first are attempting to reproduce the MATLAB version of the scripts. The original MATLAB script is also included in the scripts/ directory, but cannot be run without the MATLAB/MEX interface. Contact bmduran at ucsd dot edu for questions.

The included conda environment is a slight modification to the RRTMG-build-env described in the README of Brian Rose's RRTMG wrapper.


References
Brian Rose's RRTMG wrapper (as part of climlab): https://github.com/climlab/climlab-rrtmg

Duran, Brandon M., et al. "A new method for diagnosing effective radiative forcing from aerosol-cloud interactions in climate models." EGUsphere 2024 (2024): 1-41.

Wall, Casey J., Trude Storelvmo, and Anna Possner. "Global observations of aerosol indirect effects from marine liquid clouds." Atmospheric Chemistry and Physics 23.20 (2023): 13125-13141.

Zelinka, Mark D., Stephen A. Klein, and Dennis L. Hartmann. "Computing and partitioning cloud feedbacks using cloud property histograms. Part I: Cloud radiative kernels." Journal of Climate 25.11 (2012): 3715-3735.
