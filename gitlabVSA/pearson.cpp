// Algorithm to compute the pearson coefficient and create the superpositions of the spectra
// this module has to many exercises and should probably be simplified
// Lennard Boeselt
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <string>
#include "pearson.h"
#include "storevariables.h"
using namespace std;

void Pearson::settheo(vector <vector <vector <double > > > &theo){
    this->theospectra=theo;
}

void Pearson::setexp(vector <vector <double > > &expspectrum){
    this->expspec=expspectrum;
}
void Pearson::optimize(){                                                   // Gradient ascent algorithm to optimize the weights
    vector <double> tmp;
    vector<vector <vector< double> > >save_dist;
    vector <vector <double > >dev;
    save_dist.push_back(this->dist);
    bool breaker=false;
    createsuperposition(0,0.0);
    computepearson();
    cout << this->currentpearson << "\n";
    convergence.push_back(this->currentpearson);
    double stepsize=global.getstepsize();
    for(unsigned int k = 1; k < global.getmaxstep()+1; ++k){
        tmp.clear();
        for(unsigned int i = 0; i < this->numberofconformers; ++i){
            createsuperposition(i,0.001);                                 // Compute the gradients numerically, gradient 1
            computepearson();
            double a = this->currentpearson;
            createsuperposition(i,-0.001);
            computepearson();
            double b = this->currentpearson;
            tmp.push_back((abs(a)-abs(b))/0.002);                         // Compute the gradients numerically, gradient 2
        }
        for(unsigned int i = 0; i < dist.size(); ++i){
            dist[i][1]=dist[i][1]+stepsize*tmp[i];                        // change the energies according to the computed gradient
        }
        dev.push_back(tmp);
        save_dist.push_back(dist);                                        // save the energies
        createsuperposition(0,0);
        computepearson();
        convergence.push_back(this->currentpearson);                      // push back the convergence
        cout << this->currentpearson << "\n";
        for(unsigned int i = 0; i < dev[k].size(); ++i){
    	    if(abs(dev[k-1][i]) < global.getconvcrit()){
                breaker=true;
            }
            else{
                breaker=false;
                break;
            }
        }
        if(breaker==true){ break;}
    }
}

vector <double> Pearson::getconvergence(){
    return this->convergence;                                               // get the convergence
}

void Pearson::setdist(vector <vector <double > > &dist){
    this->dist=dist;                                                        // set the distribution
}
void Pearson::normalizeexp(){                                               // Normalize the experimental spectrum
    double max=0;
    for(unsigned int i = expbound1; i < expbound2; ++i){
        if(abs(expspec[i][1]) > max){
            max=abs(expspec[i][1]);                                     
        }
    }
    for(unsigned int i = 0; i < expspec.size(); ++i){
        expspec[i][1]=expspec[i][1]/max;
    }
}

void Pearson::createsuperposition(unsigned int number, double amount){      // create superposition based on the energies
    //cout << "1''" << "\n";
    double sum = 0;
    if(global.getweighting()==0){
        for(unsigned int i = 0; i < this->numberofconformers; ++i){
            if(number==i){
                sum = sum + exp(-(dist[i][1]+amount)/2.479);
             }
            else{
                sum = sum + exp(-dist[i][1]/2.479);
            }
        }
    }
    else if(global.getweighting()==1){
        for(unsigned int i = 0; i < this->numberofconformers; ++i){
		    if(number==i){
		    	sum=sum+abs(dist[i][1]+amount);                             //abs to make sure that dist > 0
		    }
		    else{
			    sum=sum+abs(dist[i][1]);                                    //to make sure that dist > 0
		    }
    	}
    }
    for(unsigned int i = 0; i < theospectrum.size(); ++i){
        theospectrum[i][1] = 0;                                             //here probably .clear() also an option, but I want to make sure that it is working correctly
    }
    for(unsigned int i = 0; i < this->numberofconformers; ++i){
   	    if(global.getweighting()==0){
       		 for(unsigned int j = 0; j < 1650; ++j){
            		theospectrum[j][0]=j;
        		if(number==i){
            			theospectrum[j][1]=theospectrum[j][1]+theospectra[i][j][1]*exp(-(dist[i][1]+amount)/2.479)/sum;
         		 }
        		else{
            			theospectrum[j][1]=theospectrum[j][1]+theospectra[i][j][1]*exp(-dist[i][1]/2.479)/sum;
        		}   
        	}
	    }
	if(global.getweighting()==1){
       		 for(unsigned int j = 0; j < 1650; ++j){
            		theospectrum[j][0]=j;
        		if(number==i){
            			theospectrum[j][1]=theospectrum[j][1]+theospectra[i][j][1]*abs(dist[i][1]+amount)/sum;

         		}
        		else{
            			theospectrum[j][1]=theospectrum[j][1]+theospectra[i][j][1]*abs(dist[i][1])/sum;
        		}   
		    }
    	}
	}	
}
vector<vector <double > > Pearson::gettheospectrum(){
    return this->theospectrum;                                              // return theoretical spectrum
}

void Pearson::normalizetheo(){
    double max=0;
    for(int i=global.getexpbound1()-50; i < global.getexpbound2()+50; ++i){
        if(theospectrum[i][0] >= global.gettheobound1() and theospectrum[i][0] <= global.gettheobound2()){
            if(max < abs(theospectrum[i][1])){
                max=abs(theospectrum[i][1]);
            }
        }
    }
    for(unsigned int i = 0; i <theospectrum.size(); ++i){
        theospectrum[i][1]=theospectrum[i][1]/max;
    }
}
void Pearson::computepearson(){                                             // compute pearson coefficient of given graph
    this->previouspearson=currentpearson;
    if(steps==0){
        for(unsigned int i = 0; i < expspec.size(); ++i){
            if( global.getexpbound1()-50 <= expspec[i][0] and expspec[i][0] <= global.getexpbound1()){
                this->expbound1=i;
            }
            if(global.getexpbound2() <= expspec[i][0] and expspec[i][0] <= global.getexpbound2()+50){
                this->expbound2=i;
            }
        }
        normalizeexp(); 
    }
    normalizetheo();
    int counter = 0;
    vector <double> tmpy;//outside if so that i dont need to reallocate
    if(steps==0){
            //cout <<"Here i am" << "\n";
        steps++;
        for(unsigned int i = this->expbound1; i < this->expbound2 ; ++i){
            if((expspec[i][0] >= global.getexpbound1()+6*counter) and (expspec[i][0] <= global.getexpbound1()+6*(counter+1))){
                tmpy.push_back(expspec[i][1]);
             //       cout << "expspec: " << expspec[i][0] << " " << expspec[i][1] << "\n";
            }
            else if((tmpy.size() > 0) and (expspec[i][0] > global.getexpbound1()+6*(counter+1))){
                counter++;
                double a_j=0;
                for(unsigned int j = 0; j < tmpy.size(); ++j){
                     a_j=a_j + tmpy[j];
                }
                this->expvector.push_back(a_j/tmpy.size());
                cout << "a_j: "<< a_j/tmpy.size() << " size: "<< tmpy.size() <<  "\n";
                tmpy.clear();
                tmpy.push_back(expspec[i][1]);
                a_j=0;
                }
            if(global.getexpbound1()+6*(counter+1) > global.getexpbound2()){break;}
            }
        this->nenner1=0;
        for(unsigned int i = 0; i < expvector.size(); ++i){
            nenner1=nenner1+expvector[i]*expvector[i];
        }
        this->nenner1=sqrt(nenner1);
        //cout << "6" << "\n";
        }

        //cout << "nenner1"<< nenner1 << "\n";
        vector <double> theovector;
        tmpy.clear(); //security clear
        theovector.clear();
        counter = 0;
        for(unsigned int i = 0; i < theospectrum.size(); ++i){
            if(theospectrum[i][0]>= global.getexpbound1()+6*counter and theospectrum[i][0] <= global.getexpbound1()+6*(counter+1)){
                //cout << "1. IF" << "\n";
                tmpy.push_back(theospectrum[i][1]);
            }
            else if(tmpy.size() >0 and theospectrum[i][0] >= global.getexpbound1()){
                //cout << "2. ELSE IF" << "\n";
                counter++;
                double a_j=0;
                for(unsigned int j = 0; j < tmpy.size(); ++j){
                    a_j+=tmpy[j];
                }
                theovector.push_back(a_j/tmpy.size());
                tmpy.clear();
                tmpy.push_back(theospectrum[i][1]);
            }
            if(global.getexpbound1()+6*(counter+1) > global.getexpbound2()){break;}
        }
        double nenner2 = 0;
        for(unsigned int i = 0; i < theovector.size(); ++i){
            nenner2=nenner2+theovector[i]*theovector[i];
        }
        //cout << "7" << "\n";
        nenner2=sqrt(nenner2); 
        //cout << "nenner2"<< nenner2 << "\n";
        double zaehler=0;
        //cout << theovector.size() << "\n";
        //cout << expvector.size() << "\n";
        for(unsigned int i = 0; i < theovector.size(); ++i){
            zaehler=zaehler+this->expvector[i]*theovector[i];
        }
        //cout << "8" << "\n";
        this->currentpearson=zaehler/(this->nenner1*nenner2);
        //cout << "9" << "\n";
    }
    void Pearson::setnumberofconformers(unsigned int n){
        this->numberofconformers=n;
    }
    vector<vector<double > > Pearson::getdist(){
        return this->dist;
    }
