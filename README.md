# sparse-brdf
 
This repository contains the code for the following paper:  
T. Tongbuasirilai, J. Unger, C. Guillemot, and E. Miandji, "A Sparse Non-parametric BRDF Model", ACM Transactions on Graphics, 2022  

The code is written in C++  
The following libraries are required to compile the code:  
- Eigen 3 (this is header-only)
- boost
- matio (https://sourceforge.net/projects/matio/)
- HDF5 (required by matio)

A solution file for Visual C++ 2022, as well as a CMake file is provided.  
We have compiled the code under Windows 10 and Ubuntu-22.04. We have not tested the code under MacOS.  

The code performs the model selection, reconstruction, and interpolation of BRDFs as described in the paper.  
For model selection, we have included a sample main file named "brdf3Dcomp.cpp"  
For reconstruction, we have included a sample main file named "brdf3Drecon.cpp"  
For interpolation, we have included a sample main file named "brdf3Dinterp.cpp"  

Note that the repository does not include the model training code. This is available upon request. However, we have included a trained ensemble named "DictEnsOrth3D.mat", which is located at "project-root/bin/DataBRDF/".  
To use the code, the publicly available data sets we used in the paper should be downloaded. After applying the BRDF transformations (i.e. log-plus and cosine-weighted-log as described in Section 3.1), the result should be placed in "project-root/bin/DataBRDF/TestSet/".  
We have also provided the transformed BRDFs in MATLAB's .mat file format, which can be downloaded from: URL.  

Enquiries should be sent to  
Ehsan Miandji at ehsan.miandji@liu.se  
or  
Tanaboon Tongbuasirilai at tanaboon.tongbuasirilai@liu.se  
