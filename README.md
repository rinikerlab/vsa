# vsa


VCD Spectra Alignment (VSA)

Lennard Böselt

Alpha version (v0.6)

The idea behind this approach is two fold:
1. Automatically identify the individual scaling factors of each normal mode. 
2. Determine the absolute stereochemistry of the given molecule.

If users of this code find any bugs, or have molecules where this given approach does not work, please share it with us.


  1. INSTALLATION GUIDE

  Requirements:
  - Tested under:
    - Centos 7.6 and Ubuntu 18.04
    - GCC compiler 4.8.5 and 5.2.0, c++11 standard

  Compilation:
  - g++ -std=c++11 *.cpp -o VSA.out
  
  2. EXAMPLE FILES


TO BE ADDED:
2. GUI
3. SI/SSO for comparison
4. ANALYSIS SCRIPTS USED IN THE ORIGINAL PUBLICATION
5. SCRIPTS TO CREATE THE CONFORMATIONAL ENSEMBLE
6. INTERFACES TO TURBOMOLE, GAUSSIAN16, DALTON
