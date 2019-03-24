// read_input.h
// Helper file for read_input.cpp, read in all gaussian log files
// Lennard Boeselt
#ifndef read_input_h
#define read_input_h
using namespace std;

class Reader{
    public:
        void setnumber(unsigned int n);
        unsigned int getnumber();
        void readboltzmann(string boltzmannpath);
        void readexp(string exppath);
        void readpath(string path);
        void readlogfiles(vector <string> path);
        void readinput(string inputfile,string boltzmannpath);
        void readexppeaks(string peaks);
	void setscoredir(string a);
        vector <vector <double > > getscore();
	//get functions
	void setscore();
        vector<vector<vector <double > > > getcombined();
        vector<vector <double > > getexp();
        vector<vector <double > > getexppeaks();
        vector<vector <double > > getboltzmann();
	vector <vector <double > > getIR(string a);
    private:
	vector <vector <double > > score;
	string scoredir;
        unsigned int numberofspectra; // number of spectra
        vector<vector<double > > boltzmann; //Spectrum number and Boltzmann energy
        vector<vector<double > > exp; //wavenumber peaks vs Intensity
        vector<vector<vector <double > > > combined; //For all theoretical spectra: Freq vs rot vs EM-Angle
        string inputfile;
        vector<vector <double > > exppeaks;
        unsigned int numberofconformers;
};
#endif
