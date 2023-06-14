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
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
extern MTRand eng;

//Generates the initial population, a monomorphic population at the ecological equilibrium
//The parameter para.n is now the patch size at the first generation
population initial(parameters para){
    population pop0;
    deme dem0;
    pheno phen0;
    //Create an individual whose phenotypic values are the initial attack rate and dispersal defined in the parameter set 
    phen0.d=para.inid;
    phen0.z=para.iniz;
    //Duplicate this individual for each spot in each deme.
    for (int j = 0; j < para.n; j++){
        dem0.indiv.push_back(phen0);
    } 
    //In each deme, the resource density is at ecological equilibrium
    dem0.envt=(exp(para.r)-exp(para.iniz*para.tau))/(exp(para.r)-1);
    for (int i = 0; i < para.nd; i++){
        pop0.idem.push_back(dem0);
    }
    return pop0;
}


//Mutations
pheno mutation(parameters para, pheno phen){//Given an individual with phenotype pheno...
    if (eng.rand()<para.pu){
         // ...it mutates on both its traits with probability para.pu
        phen.d+=gaussdev()*para.sigmu;
        phen.z+=gaussdev()*para.sigmu;
        //We truncate the distribution of resulting phenotypes such that the attack rate remains positive and dispersal remains between 0 and 1
        if (phen.d<0){
            phen.d=0;
        }
        if (phen.d>1){
            phen.d=1;
        }
        if (phen.z<0){
            phen.z=0;
        }
    }
    return phen;
}

population iteration(parameters para, population pop){
    vector <pheno> ofdeme[para.nd];// Declare an array for the fecundity of the population
    population popof; // Declare the offspring population       
    for (int i = 0; i < para.nd; i++){ //For each patch...
    deme dem0;
        double sumconso = 0;
        for (int j=0; j < pop.idem[i].indiv.size(); j++){//...compute the total consumption in the patch
            sumconso+=pop.idem[i].indiv[j].z;
        }
        int ni=pop.idem[i].indiv.size();
        dem0.envt=pop.idem[i].envt*(exp(para.r))/(pop.idem[i].envt*(exp(para.r)-1)+exp(sumconso*para.tau/max(1,ni))); 
        // Determines the state of the environment at the beginning of the next generation
        // NB: We scale para.tau such that the consumption time is inversely proportional to the number of individual in the patch. This means that the more individuals in a patch, the less they can accumulate resources
        popof.idem.push_back(dem0);
        for (int j = 0; j < pop.idem[i].indiv.size(); j++){//For each individual...
            double fec=para.f0+1/(pow(para.sigma,2))*(pop.idem[i].indiv[j].z/sumconso)*pop.idem[i].envt*(1-exp(-sumconso*para.tau/max(1,ni)))*(1-pop.idem[i].indiv[j].z);
            //..compute its fecundity...
            int nbabies=poisdev(fec);
            //...and sample its number of offspring in a poisson distribution with mean value its fecundity
            for (int baby=0;baby<nbabies;baby++){//For each offspring...
                if (eng.rand()<pop.idem[i].indiv[j].d){//...determine if it disperse
                    if (eng.rand()>para.cd){// For each dispersing offspring, determine if it survives...
                        int x2=i;
                        while (x2==i){
                        x2=eng.rand(para.nd); //...and if it does attribute it to a randomly sampled patch that is different from its natal patch
                        }
                        ofdeme[x2].push_back(pop.idem[i].indiv[j]); 
                    }
                }
                else{//if the offspring does not disperse, attribute it to its natal patch
                ofdeme[i].push_back(pop.idem[i].indiv[j]);
                }
            }
        }
    }
  
    for (int i = 0; i < para.nd; i++){ // For each patch...

        for (int j = 0; j<ofdeme[i].size() ;j++) {//...consider every offspring that has been attributed to it...
            if(eng.rand()<(1/(1+para.chi*ofdeme[i].size()))){//... and determine if it survive density-dependent competition
                popof.idem[i].indiv.push_back(mutation(para,ofdeme[i][j]));
            }
        }
    }
    return popof;
}


void generation(parameters para){
// Create a file to keep track of the mean and variance in the attack rate, dispersal and resource density in the population
   char nomFichier[200];
	stringstream nom;
	nom << "C:\\Users\\PATH\\variables.txt";
	nom >> nomFichier;
   ofstream myfile(nomFichier);
//Initiate the population
    population pop=initial(para);
//Compute and write in file the size of the population, as well as the mean and variance in the attack rate, dispersal and resource density for the initial population
    double meanz=para.iniz, varz=0, meand=para.inid, vard=0, meane=pop.idem[1].envt, vare=0;
    myfile <<0<<" "<<meanz<<" "<<varz<<" "<<meand<<" "<<vard<<" "<<meane<<" "<<vare<<" "<<para.nd*para.n<<endl;
    cout<<"mean in z="<<meanz<<" var in z="<<varz<<" mean in e="<<meane<<" var in e="<<vare<<endl;   
    for (int it = 0; it < para.gen; it++){   
        // iterate the life-cycle for para.gen generations
        cout<<endl<<it;
        pop=iteration(para,pop);
       meanz=0,varz=0, meand=0,vard=0, meane=0,vare=0;
       //Compute and write in file the population size, the mean and variance in the attack rate, dispersal and resource density for the population
       int nindiv=0;
        for (int i = 0; i < para.nd; i++){
            for (int j = 0; j < pop.idem[i].indiv.size() ; j++){
                nindiv+=1;
                meanz+=pop.idem[i].indiv[j].z;
                varz+=pow(pop.idem[i].indiv[j].z,2);
                meand+=pop.idem[i].indiv[j].d;
                vard+=pow(pop.idem[i].indiv[j].d,2);
            }
            meane+=pop.idem[i].envt/para.nd;
            vare+=pow(pop.idem[i].envt,2)/para.nd;
        }
        meanz/=nindiv;
        meand/=nindiv;
        vard=vard/nindiv-pow(meand,2);
        varz=varz/nindiv-pow(meanz,2);
        vare-=pow(meane,2);
        myfile <<it+1<<" "<<meanz<<" "<<varz<<" "<<meand<<" "<<vard<<" "<<meane<<" "<<vare<<" "<<nindiv<<endl;
        if ((it+1)%para.sample==0)
                 {//Each para.sample generation, sample the entire population and write in files the dispersal and attack rate of every individual, the patch it lives in, and resource density in each patch
                 char nomFichierz[100];
	             stringstream nomz;
                 nomz << "C:\\Users\\PATH\\sampling_z_generation_"<<it+1<<".txt";
	             nomz >> nomFichierz;
                 ofstream myfilez(nomFichierz);

                char nomFichierd[100];
	             stringstream nomd;
                 nomd << "C:\\Users\\PATH\\sampling_d_generation_"<<it+1<<".txt";
	             nomd >> nomFichierd;
                 ofstream myfiled(nomFichierd);
 
                 char nomFichieri[100];
	             stringstream nomi;
                 nomi << "C:\\Users\\PATH\\sampling_i_generation_"<<it+1<<".txt";
	             nomi >> nomFichieri;
                 ofstream myfilei(nomFichieri);
                
                 char nomFichiere[200];
	             stringstream nome;
                 nome << "C:\\Users\\PATH\\sampling_e_generation_"<<it+1<<".txt";
	             nome >> nomFichiere;
                 ofstream myfilee(nomFichiere);

                 for (int i = 0; i < para.nd; i++){
                 for (int j = 0; j < pop.idem[i].indiv.size() ; j++){
                    myfilez<<pop.idem[i].indiv[j].z<<endl;
                    myfiled<<pop.idem[i].indiv[j].d<<endl;
                    myfilei<<i<<endl;
               
                 }
                 myfilee <<pop.idem[i].envt<<endl;
             }
         }
        
    }  

}