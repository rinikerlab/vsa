#include <vector>
#include <sstream>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "storevariables.h"
#include "writer_output.h"
//#include "QDebug"
//File to write out all interesting stuff
//Lennard Boeselt
using namespace std;

void Writer::setnumber(unsigned int n){
    this->numberofconformers=n;
}
void Writer::write(int n, string outputname){
    if(n==1){
        stringstream ss;
        ss << global.getoutputdir() <<"/" << outputname << ".txt";
        string s = ss.str();
        ofstream fout;
        fout.flush();
        fout.open(s,ios::out);
        for(unsigned int j = 0; j < 1650; ++j){
            fout << theospectrum[0][j][0] << " " << theospectrum[0][j][1] <<"\n";
                }
        fout.close();
        
	stringstream aa;
	aa << global.getoutputdir() << "/" << outputname << "score.txt";
	string a = aa.str();
	fout.flush();
	fout.open(a,ios::out);
	for(unsigned int i = 0; i < scorevec.size() ; ++i){
		fout << i << " " << scorevec[i] << "\n";
	}
	fout.close();
	}
    else{
    normalizetheo();
    for(unsigned int i = 0; i < numberofconformers; ++i){
        stringstream ss;
        //qDebug() << "vor dem If" << "\n";
        if(n==2){

          //  qDebug() << "Im If" << "\n";
            ss << global.getoutputdir() <<"/"<< i << ".txt";
        }

        else if(global.getinverse()==0 and n!=2){
            ss << global.getoutputdir() <<"/+" << i << ".txt";
        }
        else if(global.getinverse()==1 and n!=2){
            ss << global.getoutputdir() <<"/-" << i << ".txt";
        }
        string s = ss.str();
        ofstream fout;
        fout.flush();
        fout.open(s,ios::out);
        for(unsigned int j = 0; j < 1650; ++j){
            fout << theospectrum[i][j][0] << " " << theospectrum[i][j][1] <<"\n";
            //qDebug() << theospectrum[i][j][0] << " " << theospectrum[i][j][1] <<"\n";
           /* fout << endl;
            fout << theospectrum[i][j][1];
            fout << endl;
            cout << "write: "<<theospectrum[i][j][0] << " " << theospectrum[i][j][1] <<"\n";*/
        }
        fout.close();
    }
    }
}
void Writer::setoutput(vector<vector<vector <double> > > &theospectrum){
    this->theospectrum = theospectrum;
}
void Writer::setexpspec(vector<vector<double > > &expspec){
    this->expspec=expspec;
}
void Writer::normalizeexp(){
    double max=0;
    //cout << "1'" << "\n";
        for(unsigned int i = 0; i < expspec.size(); ++i){
           if(expspec[i][0] >= global.getexpbound1() and expspec[i][0] <=global.getexpbound2()){
                if(abs(expspec[i][1]) > max){
                        max=abs(expspec[i][1]);
  			cout << "MAX:"  << " ";
                        cout << max << "\n";
           }
          }
        }
        for(unsigned int i = 0; i < expspec.size(); ++i){
        expspec[i][1]=expspec[i][1]/max;
    }

}
void Writer::setdist(vector<vector<double > > &a){
    this->dist= a;
}
void Writer::writedist(string outputname){
    stringstream ss;
    if(global.getinverse()==0){
    ss << global.getoutputdir() <<"/+" << outputname << ".txt";
    }
    else{
        ss << global.getoutputdir() <<"/-" << outputname << ".txt";
    }
    string s = ss.str();
    ofstream fout;
    fout.flush();
    fout.open(s,ios::out);
    for(unsigned int j = 0; j < this->dist.size(); ++j){
        fout << this->dist[j][0] << " " << dist[j][1] <<"\n";
        }
    fout.close();
}
void Writer::setconvergence(vector<double> convergence){
    this->convergence=convergence;
}

void Writer::writeconvergence(string a){
    stringstream ss;
    ss << global.getoutputdir() << "/" << a << "convergence.txt";
    string s = ss.str();
    ofstream fout;
    fout.flush();
    fout.open(s,ios::out);
    for(unsigned int i = 0; i < this->convergence.size(); ++i){
        fout << this->convergence[i]<<"\n";
    }
    fout.close();
}

void Writer::writeexp(){
    stringstream ss;
    ss << global.getoutputdir() <<"/" << "exp" << ".txt";
    normalizeexp();
    string s = ss.str();
    ofstream fout;
    fout.flush();
    fout.open(s,ios::out);
    for(unsigned int i = 0; i < expspec.size(); ++i){
        if(expspec[i][0]>=global.getexpbound1() and expspec[i][0] <=global.getexpbound2()){
            fout << expspec[i][0] << " "<<expspec[i][1] <<"\n";
        }
    }
    fout.close();

}

void Writer::setscore(vector<double > &a){
	this->scorevec=a;
}

void Writer::normalizetheo(){
    //cout << "1'" << "\n";
    for(unsigned int i = 0; i < numberofconformers; ++i){
    double max=0;
        for(int j = global.gettheobound1(); j < global.gettheobound2(); ++j){
            if(abs(theospectrum[i][j][1]) > max){
             max=abs(theospectrum[i][j][1]);
             cout << "MAX_THEO: "<<max << "\n";
            }
        }
        for(unsigned int j = 0; j < theospectrum[i].size(); ++j){
                theospectrum[i][j][1]=theospectrum[i][j][1]/max;
        }
   }
}
