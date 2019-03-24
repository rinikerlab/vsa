#include <vector>
#include <cmath>
#include <iostream>
#include "spectrum.h"
#include "needleman.h"
using namespace std;
// Lennard Boeselt

void VSA::settheo(vector<vector<vector <double > > > &combined, int inverse){
          this->combined=combined;
          if(inverse==1){
            for(unsigned int i = 0; i < numberofconformers; ++i){
                for(unsigned int j = 0; j < combined[i].size(); ++j){
                    //cout << "inverse"<<"\n";
                    this->combined[i][j][1] = -this->combined[i][j][1];
                    this->combined[i][j][2] = 180-this->combined[i][j][2];
                }
            }
          }
    }
void VSA::setexp(vector<vector <double > > &exppeaks){
            this->exppeaks=exppeaks;
}
vector<vector<vector <double > > >  VSA::gettheospectrum(){
    return this->theospectrum;
}
void VSA::setnumberofconformers(unsigned int n){
    this->numberofconformers = n;
}
        void VSA::generatespectrum(){
            for(unsigned int i = 0; i < this->theospectrum.size(); ++i){
                for(unsigned int j = 0; j < this->theospectrum[i].size(); ++j){
                    this->theospectrum[i][j][0]=0;
                    this->theospectrum[i][j][1]=0;
                }
            }
        for(unsigned int k = 0; k < aligned.size(); ++k){
           for(unsigned int i = 0; i < aligned[k].size(); ++i){
                for(unsigned int j = 0; j < 1650; ++j){
                    this->theospectrum[k][j][0]=j;
                    this->theospectrum[k][j][1]=this->theospectrum[k][j][1]+j/(2.236*pow(10,-39))*aligned[k][i][1]*(3/M_PI)/((aligned[k][i][0]-j)*(aligned[k][i][0]-j)+9);
                    //this->theospectrum[k][j][1]=this->theospectrum[k][j][1]+aligned[k][i][1]/(1+((aligned[k][i][0]-j)/3.)*((aligned[k][i][0]-j)/3));
                }
           }
        }
        }
        void VSA::setshifted(vector<vector <double > > a){
            this->shifted=a;
        }
        vector <double> VSA::align(){
            vector<vector <double > > tmp;
            vector<vector <double > > tmp2;
	    vector<double> scorevec;
            for(unsigned int i = 0; i < numberofconformers; ++i){//combined.size(); ++i){
                double score = needleman::Needleman(combined,i,exppeaks,tmp);
                needleman::shifting(tmp,tmp2);
		scorevec.push_back(score);
		/*for(unsigned int j = 0; j < tmp2.size(); ++j){
                    for(unsigned int k = 0; k < tmp2[i].size(); ++k){
                    cout << tmp2[j][k] << " ";
                    }
                    cout << endl;
                }*/
                aligned.push_back(tmp2);
                tmp2.clear();
                tmp.clear();
                //cout << score << endl;
            }
            tmp2.clear();
            tmp.clear();
            generatespectrum();
            return scorevec;
	}
        vector<vector<vector<double > > > VSA::getaligned(){
            return theospectrum;
        }
        void VSA::print(){
            for(unsigned int k = 0; k < theospectrum.size(); ++k){
            for(unsigned int i=0; i < theospectrum[k].size(); ++i){
                cout << theospectrum[k][i][0] << " " << theospectrum[k][i][1] << endl;
                }
            }
        }
