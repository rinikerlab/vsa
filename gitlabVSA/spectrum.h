// spectrum.h
// Lennard Boeselt
#ifndef spectrum_h
#define spectrum_h
using namespace std;
class VSA{
    public:
        void settheo(vector<vector<vector <double > > > &combined,int inverse);
        void setnumberofconformers(unsigned int n);
        void setexp(vector<vector <double > > &exppeaks);
        void generatespectrum();
        void setshifted(vector<vector <double > > a);
        void print();
        vector<double> align();
        vector<vector<vector <double > > > gettheospectrum();
        vector<vector<vector<double > > > getaligned();
    private:
        unsigned int numberofconformers;
        vector<vector<vector<double > > > combined;
        vector<vector <double > > exppeaks;
        vector<vector <double > > shifted;
        vector<vector<vector <double > > >theospectrum{300, vector<vector<double> > (1650, vector<double>(2))};
        vector<vector<vector <double > > >aligned;
};

#endif
