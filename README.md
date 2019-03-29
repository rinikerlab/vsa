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
  - The Example archiv contains all rigide molecules and methylcyclohexanon presented in the original publication. Just create an alias to your ~/.bashrc, change to the desired molecule and run VSA input.
  
  3. OUTPUT FILES
  - The calculation provides the following output files:
  
    a) exp.txt, the normalized experimental spectrum in the range [a,b]
    
    b) +distafter.txt and -distafter.txt, the energy after optimization
    
    c) +convergence.txt and -convergence.txt, the convergence of the optimzation
    
    d) +0.txt, +1.txt, ..., +n.txt and -0.txt, -1.txt, ..., -n.txt, the aligned conformer to the experimental spectrum.
    
    
TO BE ADDED:
1. GUI
2. SI/SSO for comparison
3. ANALYSIS SCRIPTS USED IN THE ORIGINAL PUBLICATION
4. SCRIPTS TO CREATE THE CONFORMATIONAL ENSEMBLE
5. INTERFACES TO TURBOMOLE, GAUSSIAN16, DALTON
