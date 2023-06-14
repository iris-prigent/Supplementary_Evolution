#include "header.h"
#include "mt.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <future>
#include <thread>
#include <algorithm>
#include <iterator>

using namespace std;

MTRand eng;

    
int main(){
    parameters par;
    par.cd=0.1;
    par.n=10;
    par.nd=2000;
    par.sigmu=0.005;
    par.pu=0.01;
    par.iniz=0.09;
    par.inid=0.36;
    par.tau=50;
    par.r=5;
    par.gen=100000;
    par.sample=500;
    par.f0=2;
    par.sigma=0.03;
    par.chi=0.08;
    generation(par);
}
