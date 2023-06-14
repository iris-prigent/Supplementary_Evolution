#include <vector>
#include <iostream>
#include <future>
#include <thread>
#include <iterator>

using namespace std;

// global variables

struct parameters
{

	int n; // Deme size
    int nd; // Number of demes
	double cd; // Cost of dispersal

    double iniz; // Initial value for the attack rate
	double inid; // Initial value for the dispersal

    double tau; // Time available for consumption (tau_1 in the main text)
    double r; // renewal rate per generation (gamma in the main text)
 
	double sigmu; // Standard error in the phenotypic effect of mutations
	double pu; // Probability of a mutation

    int gen,sample; // Number of generations and frequency of sampling

};

struct sample
{
	// This vector allows us to sample the population at specific generations (in addition to every s generation)
	vector <int> sampled;
};

struct pheno
{
	//Definition of a phenotype as two traits: dispersal and consumption 
	double d;
	double z;
};

struct deme
{
	//Definition of a deme as a collection of individuals with a shared environment
	vector <pheno> indiv;
	double envt;
};

struct population
{
	// Definition of a population as a collection of demes
	vector <deme> idem;

};

// Functions

double gaussdev();
population initial(parameters para);
pheno mutation(parameters para, pheno phen);
population iteration(parameters para, population pop);
sample sampling(parameters para);
void generation(parameters para);