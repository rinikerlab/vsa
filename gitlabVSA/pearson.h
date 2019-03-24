//pearson.h
//Helper file for pearson.cpp, which implements the computation of the pearson coefficient and the optimiyation
//Lennard Boeselt
#ifndef pearson_h
#define pearson_h
using namespace std;
class Pearson{
    public:
        void settheo(vector <vector <vector <double > > > &theo);       // Set the theoretical spectrum
        void setexp(vector <vector <double > > &expspectrum);           // Set the experimental spectrum
        void setdist(vector<vector<double > > &dist);                   // Set the QM distribution
        vector <vector <double > > getdist();                           // Obtain the optimized distribution
        void normalizeexp();                                            // normalizes the exp spectrum
        void normalizetheo();                                           // normalizes the theo spectrum
        void optimize();                                                // optimize the weights
        void createsuperposition(unsigned int number, double amount);   // computes the superposition of the theoretical spectrum of n conformer spectrum
        void computepearson();                                          // computes the pearson coefficient 
        void setnumberofconformers(unsigned int n);                     // sets the number of conformers
        vector <double> getconvergence();                               // obtain the convergence
        vector<vector <double > > gettheospectrum();                    // obtain the theoretical spectrum
    private:
        vector <vector<double > > dist;                                 // optimized distribution
        vector <vector<double > > expspec;                              // experimental spectrum
        vector <vector <vector< double > > > theospectra;               // theoretical spectrum of n conformers
        vector <double> convergence;                                    // convergence
        double numberofconformers;                                      // number of conformers
        double currentpearson=0;                                        // current pearson
        double previouspearson=0;                                       // previous pearson
        double expbound1=0;                                             // first boundary provided by the input file
        double expbound2=0;                                             // second boundary provided by the input file
        double zaehler=0;                                               
        double nenner1 = 0;
        int steps=0;
        vector <double> expvector;
        vector <vector <double > > theospectrum {(2200), vector<double>(2)};    // combined theoretical spectrum
};
#endif
