/*
 Needleman-Wunsch Implementation
 Lennard Boeselt
*/


#include <vector>
#include <typeinfo>
#include <iostream>
#include "storevariables.h"
#include <stdlib.h>
#include <algorithm>
#include <string>
#include <cmath>
#include "needleman.h"
#include <fstream>
using namespace std;

// Aligns the provided peaks of the experimental spectrum with the peaks of the theoretical spectrum
double needleman::Needleman(vector<vector<vector<double > > > &combined, int spectrumnumber, vector<vector <double> > &exppeaks, vector<vector <double> > &aligned){
    unsigned int n = combined[spectrumnumber].size()+1;             // number of theoretical peaks
    unsigned int m = exppeaks.size()+1;                             // number of experimental peaks
    vector<vector<double > > al_mat(m,vector<double > (n) );        // score for each matrix entry
    vector<vector<string > > p_mat(m,vector<string > (n) );         // Matrix for Backtracking
    for(unsigned int i = 0; i < m; ++i){                            // Assign the first row and the first column with no-matches
        al_mat[i][0] = 0;
        p_mat[i][0] = 'V';
    }
    for(unsigned int i = 0; i < n; ++i){
        al_mat[0][i] = 0;
        p_mat[0][i] = 'H';
    }
    p_mat[0][0]="S";                                                    // Start point
    double normalize = 0;                                               // init normalization value
    for(unsigned int i = 0; i < combined[spectrumnumber].size(); ++i){
	    if(abs(combined[spectrumnumber][i][1]) > normalize){
		    normalize = abs(combined[spectrumnumber][i][1]);            // find maximum value
	    }
    } 
    for(unsigned int i = 1; i < m; i++){                                // iterate through experimental peaks
        for(unsigned int j = 1; j < n; j++){                            // iterate through theoretical peaks
           double di=Diagonal(combined[spectrumnumber][j-1][0],combined[spectrumnumber][j-1][1],combined[spectrumnumber][j-1][2],exppeaks[i-1][0],exppeaks[i-1][1],normalize);                                  // Compute matched value
	       di=al_mat[i-1][j-1]+di;
           double ho=al_mat[i][j-1]+0;                                  // Compute unmatch values
           double ve=al_mat[i-1][j]+0;
           al_mat[i][j] = max(di,max(ho,ve));                           // take maximum
           p_mat[i][j] = Pointer(di,ho,ve);
        }
    }
    double checker = Backtrace(p_mat,al_mat, combined,exppeaks, aligned,spectrumnumber);    // perform backtrace
    double returnvalue=al_mat[m-1][n-1];                                                    // return maximum score
    cout << "al_mat: "<<al_mat[m-1][n-1] << "\n";
    al_mat.clear();                                                                         // clear up
    p_mat.clear();
    return returnvalue;
}

double needleman::Backtrace(vector<vector <string > >& p_mat, vector<vector <double > > &al_mat, vector<vector<vector<double > > >&combined, vector<vector <double > >&exppeaks, vector<vector <double> > &aligned,int spectrumnumber){
        // Implements the Backtrace algorithm
        unsigned int n = combined[spectrumnumber].size();
        unsigned int m = exppeaks.size();
        vector <double> tmp; 
	    double score = 0;
	    while(true){                            //iterates through p_mat, until one hits the "S" buttom
                if(p_mat[m][n]=="D"){
                        tmp.push_back(combined[spectrumnumber][n-1][0]);
                        tmp.push_back(combined[spectrumnumber][n-1][1]);
                        tmp.push_back(exppeaks[m-1][0]);
		            	tmp.push_back(exppeaks[m-1][1]);
                        n--;
                        m--;
                  }
                else if(p_mat[m][n]=="V"){
                        tmp.push_back(0);
                        tmp.push_back(0);
                        tmp.push_back(exppeaks[m-1][0]);
			            tmp.push_back(exppeaks[m-1][1]);
                        m--;
                }
                else if(p_mat[m][n] == "H"){
                        tmp.push_back(combined[spectrumnumber][n-1][0]);
                        tmp.push_back(combined[spectrumnumber][n-1][1]);
                        tmp.push_back(0);
			            tmp.push_back(0);
                        n--;
                }       
                else{
                        break;
                }
                aligned.push_back(tmp);
                tmp.clear();
        }
	return 0;
}

double needleman::Diagonal(double s1x, double s1y, double emangle, double s2x, double s2y,double normalize){
    // Computes the score
    // firstterm is the wavenumber domain
    // secondterm is the intensity domain
    if(global.getalgorithm()==0){
        double emangle_save = emangle*M_PI/180.;
        double firstterm = -1;
        double secondterm = 0;
        double a=0.028682557873482607;
        double b=0.02621550108552149;
        double c=0.9692913138559619;	
	    //firstterm=0;
	    double mean=0.99;
	    double sigma =0.00021946117399612649*8;
	    if(acceptance(s1x,s2x)){
	        firstterm=1./(sqrt(2*M_PI)*sigma)*exp(-(s2x/s1x-mean)*(s2x/s1x-mean)/(2*sigma)); //sqrt(2*M_PI*sigma)
        }
        secondterm=0;
        secondterm=abs(s2y/M_PI*s1x/(2.236*pow(10,-39))*s1y);
	    return firstterm*secondterm;
    }
}
bool needleman::acceptance(double s1x, double s2x){
    bool accept = false;
    std::cout << global.getupwardsdownwards();
    if(global.getupwardsdownwards()==0 and s2x/s1x <= 1.00 and abs(s2x-s1x) < global.getglobalcutoff()){
        accept = true;
    } 
    if(global.getupwardsdownwards()==1 and abs(s2x-s1x) < global.getglobalcutoff()){
        accept = true;
    } 
    return accept;
}
string needleman::Pointer(double di, double ho, double ve){
    // Helper function for the needleman wunsch algorithm
    double pointer = max(di,max(ho,ve));
    if(di == pointer){ 
       return "D";
    }
    else if(ho==pointer){
        return "H";
    }
    else{
        return "V";
    }

}
void needleman::shifting(vector<vector<double > > &aligned, vector<vector<double > > &shifted){
    // Function, which shifts the unmatched peaks
    vector<double> tmp;                                 //aligned 0 -> theox, aligned 1 -> theox, aligned 2 -> expx, aligned 3 -> expy
    vector < vector < double > > shiftlog;
    for(unsigned int i = 0; i < aligned.size(); ++i){
        if(aligned[i][2]!=0 and aligned[i][0]!=0){        // == if unmatched
            vector <double > tmp2;
            tmp.push_back(aligned[i][2]);
	        tmp.push_back(aligned[i][1]);
	        tmp2.push_back(aligned[i][0]);                  //theoretical peak
	        tmp2.push_back(aligned[i][2]/aligned[i][0]);    //exp divided by theo, to get the shift factor
	        shifted.push_back(tmp);
	        shiftlog.push_back(tmp2);
	        tmp2.clear();
        }
        else if(aligned[i][2]==0 and aligned[i][0]!=0){
            double min=99999;
            double shiftvalue=1;
	        double factor=1;
 	        vector <double> tmp2; 
            for(unsigned int j = 0; j < aligned.size(); ++j){
                if(aligned[j][2]!=0 and aligned[j][0]!=0){
                    if(abs(min) > abs(aligned[j][0]-aligned[i][0])){//minimaler abstand zwischen theox1 and theox2
                        min=aligned[j][0]-aligned[i][0]; //substract theox1 and theox unassigned
                        shiftvalue=aligned[j][2]/aligned[j][0]; //shiftvalue
			            factor = abs(aligned[j][3]/aligned[j][1]); // to multiply with theory
                    }
                }
            }
            tmp.push_back(aligned[i][0]*shiftvalue); //modification: multiplication
            tmp.push_back(aligned[i][1]); //factor);
            shifted.push_back(tmp);
	        tmp2.push_back(aligned[i][0]); //theoretical peak
	        tmp2.push_back(shiftvalue); //exp divided by theo
	        shiftlog.push_back(tmp2);
	        tmp2.clear();           //clean up
        }
        tmp.clear();                //clean up
     }
     tmp.clear();                   //clean up
}

