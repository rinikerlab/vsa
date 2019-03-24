// write_output.h
// Class to write out all interesting files
// Lennard Boeselt
#ifndef write_output_h
#define write_output_h
using namespace std;
class Writer{
    public:
        void setnumber(unsigned int n);
	    void setscore(vector <double> &a);
        void write(int n, string outputname);
        void setconvergence(vector <double> convergence);
        void writeconvergence(string a);
        void setexpspec(vector<vector<double > > &expspec);
        void writeexp();
        void normalizetheo();
        void normalizeexp();
        void setdist(vector<vector<double > > &a);
        void writedist(string outputname);
        void setoutput(vector<vector<vector<double > > > &theospectrum);
    private:
        vector<vector<vector<double > > > theospectrum;
        vector<vector <double > > dist;
	vector <double> scorevec;
        vector <double> convergence;
        vector<vector<double > > expspec;
        unsigned int numberofconformers;
};

#endif
