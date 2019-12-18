An updated python script for the code is going to be uploaded at the end of january, with a simplified Input Output scheme. The old files can be found in VSA_old.

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
  
  2. OUTPUT FILES
  - The calculation provides the following output files:
  
    a) exp.txt, the normalized experimental spectrum in the range [a,b]
    
    b) +distafter.txt and -distafter.txt, the energy after optimization
    
    c) +convergence.txt and -convergence.txt, the convergence of the optimzation
    
    d) +0.txt, +1.txt, ..., +n.txt and -0.txt, -1.txt, ..., -n.txt, the aligned conformer to the experimental spectrum.
    
