Variable Phase Approach
The Program files are as follows:
The first 4 files form the core of the program.

1. VPA-main.f90  : This is the main program
2. VPA-Lib.f90   : This file contains many of the supporting subroutines
3. DLSODALIB.f   : This file contains the LSODA code for solving non-linear 
                   ordinary differential equations. 
4. RICBES.FOR    : The subroutine computes the spherical Bessel Functions 
                   and their derivatives.

The user must supply a subroutine to compute an interaction potential and its 
derivative. An example potential code is supplied.

5: Potential.f90   : 
If parameter KP=1 this subroutine computes a Lennard-Jones attractive potential and its derivative with parameters for a Ne-Ne collision otherwise, 
if KP is any other integer this subroutine computes an arbitrary repulsive potential as described in the write-up.

The program requires an input file. The format of this file is described in the 
text. Three sample input files are provided for the two modes in which the program 
can function. The files must be renamed as VPA_input.txt when used with the code.

7. VPA_input - Escan.TXT  :  This is an example of an input file for a scan over 
                             scattering energies.
8. VPA_input - Lscan.TXT  :  This is an example of an input file for a scan over 
                             orbital angular quantum numbers L.
9. VPA_input 10^5cm-1.TXT :  This is another example of an input file for a scan over 
                             orbital angular quantum numbers L, but at a much higher 
                             collision energy than in the above example (file 8).

The program produces several files. The only important file is LPH.txt. This 
contains four columns 
                 1. The scattering energy in cm-1; 
                 2. The VPA phase shift divided by pi.
                 3. The VPA partial scattering cross section.
                 4. The WKB phase shift divided by pi.
Four example output files corresponding to the three input files above are 
supplied. The three example outputs below correspond to the use of the attractive 
Lennard-Jones potential of file 5.

10. LPH - Escan.TXT : Output file corresponding to input file 7 above.
11. LPH - Lscan.TXT : Output file corresponding to input file 8 above.
12. LPH-a-10^5cm^-1.TXT : Output file corresponding to input file 9 above.

One further example of an output file is provided. This corresponds to the 
use of the repulsive potential (file 6) and input file 9 above.

13. LPH-r-10^5cm^-1.TXT : Output file corresponding to input file 9 above together 
                         with the repulsive potential of file 6.

