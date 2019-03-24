// Implements the main function
// Lennard Boeselt
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "read_input.h"
#include "spectrum.h"
#include "needleman.h"
#include "storevariables.h"
#include "writer_output.h"
#include "pearson.h"
#include "storevariables.h"
storevariables global;
using namespace std;
vector <vector <double > > findpeaks(vector <vector <double  > > &IR){  // Helper function and quick hack to compute the peak position if not provided
	vector <vector <double > > peaks;
	vector <double > tmp;
	double max = 0;
	for(unsigned int i = 0; i < IR.size(); ++i){
		if(abs(IR[i][1]) > max and IR[i][0] >= global.getexpbound1() and IR[i][0] <=global.getexpbound2()){     
			max = abs(IR[i][1]);
		}
	}
	for(unsigned int i = 1; i < IR.size()-1; ++i){
		if(abs(IR[i][1]) > abs(IR[i-1][1]) and abs(IR[i][1]) > abs(IR[i+1][1]) and IR[i][0] >= global.getexpbound1() and IR[i][0] <=global.getexpbound2()){ // simply peak criterium
			tmp.push_back(IR[i][0]);
			tmp.push_back(IR[i][1]);
			peaks.push_back(tmp);
			tmp.clear();
		}
	}
	return peaks;

}

void run(){
   if(true){                                                    // IMPORTANT KEYWORD FOR THE GUI
        Reader reader;
        string inputfile1 = global.gettheoretical();            // get global variables for the computation
        string inputfile2 = global.getboltzmann();
        string inputfile3 = global.getexperimentalpeaks();
	    vector <vector <double > > a;
        string outfile1;
        string outfile2;
        if(global.getinverse()==0){                             // set names of output files
            outfile1 = "+beforeoptimization";
            outfile2 = "+afteroptimization";}
        else{
            outfile1="-beforeoptimization";
            outfile2="-afteroptimization";
        }
        int numberofconformers=global.getnumberofconformers(); // get number of conformers
        reader.setnumber(global.getnumberofconformers());
        reader.readinput(inputfile1,inputfile2);
        reader.readexppeaks(inputfile3);
        vector<vector<double > > boltzmann=reader.getboltzmann();
	    if(global.getweighting()==1){	
		    reader.setscoredir(global.getscoredir());
		    reader.setscore();
		    boltzmann=reader.getscore();
		    for(unsigned int i=0; i < boltzmann.size(); ++i){
			    cout << "score: "<<boltzmann[i][1] << "\n";}
	    }
        vector<vector <double > > exppeaks = reader.getexppeaks();
        vector<vector<vector <double > > > combined=reader.getcombined();
        vector<vector <double > > exp=reader.getexp();
        string IR = global.getIRdir();
	    if(false){ 
		    cout << "Read in peaks from IR spectrum "<<"\n";
		    vector <vector <double > > IR = reader.getIR(global.getIRdir());
	    	a=findpeaks(IR);
		    exppeaks.clear();
		    exppeaks=a;	
	    }
	    if(true){
		    vector <vector <double > > a;
		    cout << "Specified peaks are used" << "\n";
		    ifstream inFile;
		    inFile.open("newpeaks.dat");            // get specified peaks
		    double x,y;
		    while(inFile >> x >> y){
			    if(x >= global.getexpbound1() and x <=global.getexpbound2()){
				    cout << "while\n";
				    vector <double > tmp;
				    tmp.push_back(x);
				    tmp.push_back(y);
				    a.push_back(tmp);
				    tmp.clear();
			    }
		    } 	
		    inFile.close();
            exppeaks.clear();                       // Delete exppeaks to make sure that it is really clear
		    exppeaks=a;                             // get experimental peaks
	    }
        // sort peaks so that the experimental peak list is in the same order as the experimental spectrum
        sort(exppeaks.begin(),exppeaks.end(), [](const vector<double>& a, const std::vector< double >& b){ return a[0] < b[0]; } );     
        sort(exp.begin(),exp.end(), [](const vector<double>& a, const std::vector< double >& b){ return a[0] < b[0]; } );
        // cout the selected peaks
	    for(int i = 0; i < exppeaks.size(); ++i){
                 cout <<"peaks: " << exppeaks[i][0] << " " << exppeaks[i][1] << "\n";
        }
        VSA vsa;
        vsa.setnumberofconformers(numberofconformers);
        vsa.settheo(combined,global.getinverse());
        vsa.setexp(exppeaks);
        vector<vector<vector <double > > >theospectrumquick{300, vector<vector<double> > (1650, vector<double>(2))};    // Maximum 300 conformers, range until 0...1650
        for(unsigned int k = 0; k < combined.size(); ++k){
           for(unsigned int i = 0; i < combined[k].size(); ++i){
                for(unsigned int j = 0; j < 1650; ++j){
                    theospectrumquick[k][j][0]=j;
                    theospectrumquick[k][j][1]=theospectrumquick[k][j][1]+combined[k][i][1]/(1+((combined[k][i][0]-j)/3.)*((combined[k][i][0]-j)/3));
                }
           }
        }
        // CREATE CLASSES TO ALIGN SPECTRUM
        Writer unalig;
        unalig.setnumber(global.getnumberofconformers());
        unalig.setoutput(theospectrumquick);
        unalig.setexpspec(exp);
        unalig.write(2, outfile1);
        vector<double> scorevec = vsa.align();
        vector<vector<vector<double> > > aligned = vsa.getaligned();
        // OUTPUT
        Writer b;
        b.setnumber(numberofconformers);
        b.setoutput(aligned);
        b.setexpspec(exp);
        b.write(0,outfile1);
        b.writeexp();
        // OUTPUT BOLTZMANN
        Pearson pearson2;
        pearson2.setnumberofconformers(numberofconformers);
        pearson2.settheo(aligned);
        pearson2.setdist(boltzmann);
        pearson2.setexp(exp);
        pearson2.createsuperposition(0,0);
        pearson2.normalizetheo();
        vector<vector <double> > tmp=pearson2.gettheospectrum();
        vector <vector <vector <double > > > tmp2;
        tmp2.push_back(tmp);
        b.setoutput(tmp2);
        b.write(1,outfile1);
        tmp.clear();    // CLEAN UP
        tmp2.clear();   // CLEAN UP
        // CREATE OPTIMIZED SPECTRUM
        Pearson pearson;
        pearson.setnumberofconformers(numberofconformers);
        pearson.settheo(aligned);
        pearson.setdist(boltzmann);
        pearson.setexp(exp);
        pearson.optimize();       
        vector <double> convergence = pearson.getconvergence();
        vector <vector <double > > tmp3;
        tmp3.clear();
        tmp.clear();
        tmp2.clear();
        // OUTPUT OPTIMIZED SPECTRUM
        tmp3=pearson.getdist();
        tmp=pearson.gettheospectrum();
        tmp2.push_back(tmp);
        Writer c;
        c.setconvergence(convergence);
        c.setscore(scorevec);
        // WRITE OUT OTHER FILES
        if(global.getinverse()==1){
            c.writeconvergence("-");
        }
        else{
            c.writeconvergence("+");
        }
        c.setdist(tmp3);
        c.writedist("distafter");
        c.setoutput(tmp2);
        c.write(1,outfile2);
	cout << "scores: " << "\n";
	for(unsigned int i = 0; i < scorevec.size(); ++i){
		cout << i<< " "<< scorevec[i] << "\n";
	}
	if(!IR.empty()){
	 std::ofstream outfilewrite("peaks.dat");
		for(unsigned int i=0; i < exppeaks.size(); ++i){
			 outfilewrite << exppeaks[i][0] << " " << exppeaks[i][1] << "\n";
		}
	 outfilewrite.close();
	}
   }
}

int main(int argc, char *argv[], char *envp[]){
 ifstream infile(argv[1]);                                          // Read input00 file
 string a;
 unsigned int i = 0;
 while(infile >> a){                                                // readin file and set global parameters. This here is very poorly done and can be easily improved, the file does not realize comments yet
    if(i==0){global.setoutputdir(a);}
        else if(i==1){global.setboltzmann(a);}                      // SET BOLTZMANN DISTRIBUTION
        else if(i==2){global.settheoretical(a);}                    // SET THEORETICAL SPECTRUM
        else if(i==3){global.setexperimentalpeaks(a);}              // SET EXPERIMENTAL PEAKS
        else if(i==4){global.setinverse(std::stoi(a));}             // COMPUTE INVERSE SPECTRUM?
        else if(i==5){global.setnumberofconformers(std::stoi(a));}  // SET NUMBER OF CONFORMERS
        else if(i==6){global.setexpbound1(std::stoi(a));}           // SET EXPERIMENTAL BOUNDARIES
        else if(i==7){global.setexpbound2(std::stoi(a));}           
        else if(i==8){global.settheobound1(std::stoi(a));}          // SET THEORETICAL BOUNDARIES (NOT USED HERE, IRRELEVANT!)
        else if(i==9){global.settheobound2(std::stoi(a));}
        else if(i==10){global.setconvcrit(std::stod(a));}   
        else if(i==11){global.setstepsize(std::stod(a));}           
        else if(i==12){global.setmaxstep(std::stod(a));}
        else if(i==13){global.setweighting(std::stoi(a));}
        else if(i==14){global.setscoredir(a);}
        else if(i==15){global.setIRdir(a);}
        else if(i==16){global.setequalformat(std::stod(a));}
        else if(i==17){global.setupwardsdownwards(std::stoi(a));}
        else if(i==18){global.setglobalcutoff(std::stod(a));}
        ++i;
    }
    run();                      // run the algorithm
    return 0;
}
