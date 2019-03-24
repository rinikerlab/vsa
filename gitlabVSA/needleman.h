//Helper file for needleman.cpp
// Lennard Boeselt
#ifndef needleman_h
#define needleman_h
using namespace std;
namespace needleman{
double Needleman(vector<vector<vector<double > > > &combined, int spectrumnumber, vector<vector <double> > &exppeaks, vector<vector <double> > &aligned); // Implements the Needleman algorithm
double Backtrace(vector<vector<string > >  &p_mat, vector<vector<double> > &al_mat,vector<vector<vector<double > > >&combined, vector<vector <double > >&exppeaks, vector<vector <double> > &aligned,int spectrumnumber);                                                                                                                                                // Implements the Backtrace algorithm
double Diagonal(double s1x, double emangle, double s1y, double s2x, double s2y,double normalize);                                                         // Helper function for the needleman wunsch algorithm
string Pointer(double di, double ho, double ve);                                                                                                          // Helper function for the needleman wunsch algorithm
//double funct(double s2x);                                                                                                                               
bool acceptance(double s1x, double s2x); // acceptance criterium
void shifting(vector<vector<double > > &aligned, vector<vector<double > > &shifted);                                                                      // Shifting function to align the unmatched peaks
};
#endif
