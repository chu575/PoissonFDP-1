# PoissonFDP-1
# This repository includes codes for solving the designed Poisson problem using the Fictitious Domain method.

# The strategy used for addressing the immersed boundary condition is Approach 1. In this Approach, the constraint (i.e., immersed 
# boundary condition) lives not only on the immersed boundary, but also valid in the area outside the needed domain.

# Main class and main functions are listed in the source code, PoissonFDP-1.cc

# namespaces and classes regarding the FD method are listed in the head file, TestCaseSetting.h and defined in the source codes, 
# MMSassignment.cc and IBimplementation.cc.
# In addition, IBimplementation.cc only defines a weight factor which describes the area where the IB condition lives.

# If one wants to use this code to solve some other test cases, he/she needs to change TestCaseSetting.h, MMSassignment.cc and 
# IBimplementation.cc

# Historical record regarding all modifications is not listed here. This repository only stores final codes.

# The created date of this repository was on 7/28/2019
