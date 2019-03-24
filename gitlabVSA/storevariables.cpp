//Storage for specified variables
//Lennard Boeselt
#include "storevariables.h"

void storevariables::setexperimental(std::string a){
    this->experimentalpath=a ;
}
void storevariables::settheoretical(std::string a){
    this->theoreticalpath=a;
}
void storevariables::setscoredir(std::string a){
	this->scoredir=a;
}
void storevariables::setexperimentalpeaks(std::string a){
    this->experimentalpeakspath=a;
}
void storevariables::setoutputdir(std::string a){
    this->outputdir=a;
}
void storevariables::setnumberofconformers(int a){
    this->numberofconformers=a;
}
void storevariables::setboltzmann(std::string a){
    this->boltzmann=a;
}
void storevariables::setmaxstep(unsigned int a){
	this->maxstep=a;
}
void storevariables::setIRdir(std::string a){
	this->IRdir=a;
}
std::string storevariables::getIRdir(){
	return this->IRdir;
}
void storevariables::setexpbound1(int a){
    this->expbound1=a;
}
void storevariables::setexpbound2(int a){
    this->expbound2=a;
}
void storevariables::setalgorithm(unsigned int a){
    this->algorithm=a;
}
void storevariables::setavalue(double a){
    this->avalue = a;
}
void storevariables::setupwardsdownwards(int a){
    this->upwardsdownwards = a;
}        
void storevariables::setglobalcutoff(double a){        
}
double storevariables::getavalue(){
    return this->avalue;
}
unsigned int storevariables::getalgorithm(){
    return this->algorithm;
}

void storevariables::settheobound1(int a){
    this->theobound1=a;
}

void storevariables::settheobound2(int a){
    this->theobound2=a;
}

std::string storevariables::getexperimental(){
    return this->experimentalpath ;
}

std::string storevariables::getoutputdir(){
    return this->outputdir ;
}
std::string storevariables::gettheoretical(){
    return this->theoreticalpath;
}
std::string storevariables::getexperimentalpeaks(){
    return this->experimentalpeakspath;
}
double storevariables::getstepsize(){
    return this->stepsize;
}
void storevariables::setconvcrit(double a){
    this->convcrit=a;
}
double storevariables::getconvcrit(){
    return this->convcrit;
}
int storevariables::getequalformat(){
	return this->equalformat;
}
void storevariables::setequalformat(int a){
	this->equalformat=a;
}
void storevariables::setstepsize(double a){
    this->stepsize=a;
}

int storevariables::getexpbound1(){
    return this->expbound1;
}

int storevariables::getexpbound2(){
    return this->expbound2;
}

int storevariables::gettheobound1(){
    return this->theobound1;
}

int storevariables::gettheobound2(){
    return this->theobound2;
}
unsigned int storevariables::getinverse(){
    return this->inverse;
}
void storevariables::setinverse(unsigned int a){
    this->inverse=a;
}

std::string storevariables::getboltzmann(){
    return this->boltzmann;
}
unsigned int storevariables::getmaxstep(){
 return this->maxstep;
}
void storevariables::setweighting(int a){
 this->weighting=a;
}
int storevariables::getweighting(){
 return this->weighting;
}
std::string storevariables::getscoredir(){
	return this->scoredir;
}
int storevariables::getnumberofconformers(){
    return this->numberofconformers;
}
int storevariables::getupwardsdownwards(){
    return this->upwardsdownwards;
}        
double storevariables::getglobalcutoff(){        
    return this->cutoff;
}

