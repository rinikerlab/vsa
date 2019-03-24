#include <vector> 
#include <sstream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include "storevariables.h"
#include "read_input.h"
using namespace std;
// READER, Lennard Boeselt
void Reader::setnumber(unsigned int n){
    numberofconformers = n;
}
vector <vector <double > > Reader::getIR(string IR){        // Read in IR spectrum
	ifstream infile;
	infile.open(IR);
	vector <vector <double > > output;
	vector<double> tmp;
	while(!infile.eof()){
		tmp.clear();
		double x;
		double y;
		infile >> x >> y;
		if(infile.eof()) break;
			tmp.push_back(x);
			tmp.push_back(y);
			output.push_back(tmp);
			//}
		}
	double max = 0;
	for(unsigned int i = 0; i < output.size(); ++i){
		if(abs(output[i][1]) > max and global.getexpbound1() <= output[i][0] and output[i][0] <= global.getexpbound2()){
			max=abs(output[i][1]);
			cout << "max: " << max << "\n";
		}
	}
	for(unsigned int i=0; i < output.size(); ++i){
		output[i][1]=output[i][1]/max;
	}
	tmp.clear();
	return output;
}
void Reader::setscoredir(string a){
	this->scoredir=a;
}
void Reader::setscore(){        // set energies and numbers of log file ?
	ifstream infile;
	infile.open(scoredir);
	vector<double> tmp;
	unsigned int i = 0;
	while(!infile.eof()){
		tmp.clear();
		double number;
		double energy;
		infile >> number >> energy;
		if(infile.eof()) break;
		tmp.push_back(number);
		tmp.push_back(energy);
		score.push_back(tmp);
		++i;	
	    }
	tmp.clear();
	infile.close();
	}
vector <vector <double > > Reader::getscore(){  // output energies and numbers of log file
	return this->score;
}
unsigned int Reader::getnumber(){               // get number of theoretical spectra
    return numberofspectra;
}
void Reader::readboltzmann(string boltzmannpath){//read in qm boltzmann distribution
            //cout << "readboltzmann"<< " ";
        ifstream infile;
        infile.open(boltzmannpath);
        unsigned int i = 0;
        vector<double > tmp;
        while(!infile.eof()){
            tmp.clear();
            double number;
            double energy;
            infile >> number >> energy;
            if(infile.eof()) break;
            tmp.push_back(number);
            tmp.push_back(energy);
            boltzmann.push_back(tmp);
            i++;
            }
        tmp.clear();
        infile.close();
}
void Reader::readexp(string exppath){ //reads one to few; Bug // Should be fixed
            //cout << "readexp"<< " ";
        ifstream infile;
        string STRING;
        infile.open(exppath);
        vector<double > tmp;
        while(!infile.eof()){
            tmp.clear();
            double x,y;
            infile >> x >> y;
            if(infile.eof()) break;
            tmp.push_back(x);
            tmp.push_back(y);
            exp.push_back(tmp);
        }
        tmp.clear();
            /*for(unsigned int i = 0; i < exp.size() ; ++i){
                cout << exp[i][0] << exp[i][1] << "\n";
            }*/
        infile.close();
}
void Reader::readpath(string path){ // Read log files, rotational strength and so on and so forth
            //cout << "readpath"<< " ";
    ifstream infile;
    string STRING;
            infile.open(path);
            vector <double> freq;
            vector <double> rot;
            vector <double> em;
            while(!infile.eof()){
                getline(infile,STRING);
                // Search keywords in gaussian log file and read in information
                if((12<STRING.size() and STRING.size() <= 73) and STRING.substr(0,18).compare(" Frequencies --  ")==1){
                    string::size_type sz;
                     try{
                            double val = stod(STRING.substr(18,27),&sz);
                            freq.push_back(val);
                            val = stod(STRING.substr(41,50),&sz);
                            freq.push_back(val);
                            val = stod(STRING.substr(64,73),&sz);
                            freq.push_back(val);
                        }
                        catch(exception& e){
                            continue;
                        }
                        STRING="";
                    }
                    if((12<STRING.size() and STRING.size() <= 73) and STRING.substr(0,18).compare(" Rot. str.   --  ")==1){
                        string::size_type sz;
                        try{
                            double val = stod(STRING.substr(18,27),&sz);
                            rot.push_back(val);
                            val = stod(STRING.substr(41,50),&sz);
                            rot.push_back(val);
                            val = stod(STRING.substr(64,73),&sz);
                            rot.push_back(val);
                        }
                        catch(exception& e){
                            continue;
                        }
                        STRING="";
                    }
                    if((12<STRING.size() and STRING.size() <= 73) and STRING.substr(0,18).compare(" E-M angle   --  ")==1){
                        string::size_type sz;
                        try{
                            double val = stod(STRING.substr(18,27),&sz);
                            em.push_back(val);
                            val = stod(STRING.substr(41,50),&sz);
                            em.push_back(val);
                            val = stod(STRING.substr(64,73),&sz);
                            em.push_back(val);
                        }
                        catch(exception& e){
                            continue;
                        }
                        STRING="";
                    }
                }
                    vector <vector <double > > dummy;
                    vector <double> tmp;
                    for(unsigned int j = 0 ; j < freq.size(); ++j){
                        tmp.push_back(freq[j]);
                        tmp.push_back(rot[j]);
                        tmp.push_back(em[j]);
                        dummy.push_back(tmp);
                        tmp.clear();
                    }
                    combined.push_back(dummy);
                    dummy.clear();
                    tmp.clear();
                    rot.clear();
                    em.clear();
                    freq.clear();
                    infile.close();
            }
        void Reader::readlogfiles(vector <string> path){//read in each specified log file
            //cout << "readlogfiles"<< " ";
            for(unsigned int i = 0; i < path.size(); ++i){
                readpath(path[i]);
            }
        }
        void Reader::readinput(string inputfile,string boltzmannpath){
            //cout << "readinput"<< " ";
            readboltzmann(boltzmannpath);
           // cout << "1"<< " ";
            ifstream infile(inputfile);
            string STRING;
            unsigned int line = 0;
            while(!infile.eof()){
                getline(infile,STRING);
             //   cout << "2"<< " ";
                if(line == 0){//boltzmann.size()
                    for(unsigned int i = 0; i < numberofconformers; ++i){ //warum hier -1?!
                        string a = STRING+"VCD"+to_string((int)boltzmann[i][0])+".log";
                        readpath(a);
                        //cout << "for loop" << endl;
                        //tmp.push_back(STRING+to_string((int)boltzmann[i][0]));
                    }
                }
               // cout << "3"<< " ";
                //cout << "leave for loop "<<endl;
                if(line==1){
                    readexp(STRING);
                }
                line=line+1;
               
            }/*
            for(unsigned int i = 0; i < 1; ++i){
                for(unsigned int j = 0; j < combined[i].size(); ++j){
                    //cout <<"{";
                    for(unsigned int k = 0; k < combined[i][j].size(); ++k){
                      //  cout << combined[i][j][k]<< " ";
                        if(k!=combined[i][j].size()-1) cout << ", " ;
                    }
                    cout <<"}," << endl;
                }
            }*/
        }
        void Reader::readexppeaks(string input){
            ifstream infile(input);
            string STRING;
            vector <double> tmp;
            while(!infile.eof()){
                double a,b;
                infile >> a >> b;
                if(infile.eof()) break;
		if(a >=global.getexpbound1() and a <= global.getexpbound2()){
                	tmp.push_back(a);
                	tmp.push_back(b);
                	exppeaks.push_back(tmp);
                }
		tmp.clear();
            }
            tmp.clear();//security clear
            /*for(unsigned int i = 0; i < exppeaks.size(); ++i){
                cout << exppeaks[i][0] << " " << exppeaks[i][1] << endl;
            }*/
        
        }
        vector<vector<vector <double > > > Reader::getcombined(){
            return combined;
        }
        vector<vector <double> > Reader::getexp(){
            return exp;
        }
        vector<vector <double > > Reader::getexppeaks(){
            return exppeaks;
        }

        vector<vector <double > > Reader::getboltzmann(){
            return this->boltzmann;
        }

/*int main(){
    Reader a;
    cout << "main ";
    a.readinput("simvastatinHF.txt","filesoutput");
    return 0;
}*/
