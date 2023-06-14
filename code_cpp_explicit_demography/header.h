#include <vector>
#include <iostream>
#include <future>
#include <thread>
#include <iterator>

using namespace std;

// global variables

//#define filePar "par"

struct parameters
{
	// Population parameters
    int nd; // # of demes
	double cd; // cost of dispersal

    double iniz; // initial consumption
	double inid; // initial dispersal
	int n; // initial deme size

    // intra- and inter-generational effects
    double tau; // consumption period
    double r; // renewal rate

	//demographic parameters
	double f0; // baseline fecundity
	double sigma; // scale the strengh of selection 
	double chi; // density dependent survival
 
	// Mutation parameters
	double sigmu; // Standard error in the phenotypic effect of mutations
	double pu; // Probability of mutation

    int gen,sample; // # of generations and freq of sampling

};

struct pheno
{
	//Definition of an phenotype as two traits: dispersal and consumption 
	double d;
	double z;
};

struct deme
{
	//Definition of a deme as a vector of individuals with a shared environment
	vector <pheno> indiv;
	double envt;
};

struct population
{
	// Definition of a population as a vector of demes
	vector <deme> idem;

};

double gaussdev();
double poisdev(const double xm);
double gammln(const double xx);

population initial(parameters para);
double benefits(parameters para, deme dem0);
double mutation(parameters para, double pheno);
population iteration(parameters para, population pop);
void generation(parameters para);