// save all variables 
// Lennard Boeselt
#ifndef STOREVARIABLES_H
#define STOREVARIABLES_H
#include <string>
class storevariables
{
public:
    //storevariables();
    //setter
    void setexperimental(std::string a);
    void settheoretical(std::string a);
    void setexperimentalpeaks(std::string a);
    void setboltzmann(std::string a);
    void setnumberofconformers(int a);
    void setmaxstep(unsigned int a);
    void setexpbound1(int a);
    void setexpbound2(int a);
    void settheobound1(int a);
    void settheobound2(int a);
    void setinverse(unsigned int a);
    void setoutputdir(std::string a);
    void setstepsize(double a);
    void setconvcrit(double a);
    void setalgorithm(unsigned int a);
    void setIRdir(std::string a);
    void setavalue(double a);
    void setscoredir(std::string a);
    std::string getscoredir();
    void setweighting(int a);
    void setupwardsdownwards(int a);
    void setglobalcutoff(double a);
    //getter
    std::string getexperimental();
    std::string gettheoretical();
    std::string getexperimentalpeaks();
    std::string getboltzmann();
    std::string getoutputdir();
    std::string getIRdir();
    unsigned int getinverse();
    unsigned int getalgorithm();
    double getstepsize();
    double getavalue();
    int getweighting();
    double getconvcrit();
    int getexpbound1();
    int getexpbound2();
    int gettheobound1();
    int gettheobound2();
    void setequalformat(int a);
    int getequalformat();
    int getnumberofconformers();
    unsigned int getmaxstep();
    int getupwardsdownwards();
    double getglobalcutoff();
private:
    std::string experimentalpath;
    std::string theoreticalpath;
    std::string experimentalpeakspath;
    std::string boltzmann;
    std::string outputdir;
    std::string scoredir;
    std::string IRdir = "";
    int weighting = 0;
    unsigned int algorithm=0;
    double convcrit = 0.001;
    unsigned int maxstep=1000;
    double stepsize = 10;
    unsigned int inverse = 0;
    int numberofconformers=1;
    int expbound1=1150;
    int expbound2=1450;
    int theobound1=1150;
    int theobound2=1450;
    int equalformat = 1;
    double avalue = 1;
    double cutoff = 40;
    int upwardsdownwards = 0;
};
extern storevariables global;
#endif // STOREVARIABLES_H
