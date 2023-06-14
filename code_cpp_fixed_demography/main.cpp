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
    // Here we define the parameter set to be used 
    par.cd=0.1; // Cost of dispersal
    par.n=10; // Deme size
    par.nd=1000; // Number of demes
    par.tau=5; // Time available for consumption (tau_1 in the main text)
    par.r=5; // renewal rate per generation (gamma in the main text)
    par.pu=0.01; // Probability of a mutation
    par.sigmu=0.005; // Standard error in the phenotypic effect of mutations
    par.iniz=0.01; // Initial value for the attack rate
    par.inid=0.1; // Initial value for the dispersal
    par.gen=100000; // Number of generations
    par.sample=500; // Frequency at which the whole population is sampled
    generation(par);
}
