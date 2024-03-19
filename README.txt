README
------
This is the README file for the MATLAB code which accompanies
the following publication:

    DiMattina, C. (2015). "Fast adaptive estimation of multidimensional 
            psychometric functions". Journal of Vision 15(9):5-5. 

If you use this code, please cite this publication. 

Version + Copyright
-------------------
Version 1.0 Beta (test version)
Copyright (C) 2015- Christopher DiMattina (cdimattina@fgcu.edu)

License Agreement
----------------- 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Detailed Description
--------------------
This code allows one to replicate the 1-D and 2-D examples from the paper. 
The scripts for each algorithm + example are:

Prior-PSI   (Fig. 4): PSI_Prior_1D.m,  PSI_Prior_2D.m
Lookup-PSI  (Fig. 8): PSI_Lookup_1D.m, PSI_Lookup_2D.m
Laplace-PSI (Fig. 9): PSI_Laplace_1D.m, PSI_Laplace_2D.m

To run each script, simply type its name at the MATLAB prompt, for example
>>PSI_Prior_1D

Functions called by these scripts are in the folder .\SupportFuns\
MatFiles they require (i.e. look-up tables) are in  .\MatFiles\

This code has been tested on MATLAB r2015a with the Optimization Toolbox

Any further questions can be e-mailed to the author (cdimattina@fgcu.edu)







   