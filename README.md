# MX_Spatial_Lib_Matlab_v1.0
A lightweight and versatile rigid body dynamics algorithm based on Roy Featherstone's spatial algorithm. Tailored for robotics application.

Matlab version of spatial dynamics. More readable, easier understand.

You can refer to the examples "MX_test1_twoLinks.m" and "MX_test2_floatBase.m" to learn how to set up spatial models and use function handles to calculate the desired matrix.

So far, it can obtain some common matrix related to kinematics and dynamics, such as absolute position and orientation of one point or body, jacobian and its derivative in the form of 6D, dynamics matrix M and F, etc. It can also calculate some centroidal matrix such as centroidal momentum matrix(A) and its derivative.

Mainly refer to Roy Featherstone's book 《Rigid Body Dynamics Algorithms》, some files about floatbase in drake's spatial library and Orin and Goswami's paper 《Centroidal dynamics of a humanoid robot》.

You can regard it as a lite matlab version of C++ library RBDL.

Written by Meng Xiang in BIT on 2021/12/2.
