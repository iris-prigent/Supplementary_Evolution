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
population initial(parameters para){
    population pop0;
    deme dem0;
    pheno phen0;
    //Create an individual whose phenotypic values are the initial attack rate and dispersal defined in the parameter set 
    phen0.d=para.inid;
    phen0.z=para.iniz;
    //Copy this individual in each spot in each deme
    for (int j = 0; j < para.n; j++){
        dem0.indiv.push_back(phen0);
    } 
    //In each deme, the resource abundance is at ecological equilibrium
    dem0.envt=(exp(para.r)-exp(para.n*para.iniz*para.tau))/(exp(para.r)-1);
    for (int i = 0; i < para.nd; i++){
        pop0.idem.push_back(dem0);
    }
    return pop0;
}



//Mutations (the function necessary to do step (iii) of the procedure described in Appendix F4)
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
    double fec[para.nd][para.n]; // Declare an array for the fecundity of the population
    population popof = pop; // Declare the offspring population

    for (int i = 0; i < para.nd; i++){ //For each deme...
        double sumconso = 0;
        for (int j=0; j < para.n; j++){//...compute the total consumption in the deme
            sumconso+=pop.idem[i].indiv[j].z;
        }
        for (int j = 0; j < para.n; j++){// Compute the fecundity of each individual in the deme
        // Step (i) of the procedure described in Appendix F4 (eq. F-15)
            fec[i][j]=(pop.idem[i].indiv[j].z/sumconso)*pop.idem[i].envt*(1-exp(-sumconso*para.tau))*(1-pop.idem[i].indiv[j].z);
        }
        // Compute the resource abundance at the beginning of the next generation
        // Step (iv) of the procedure described in Appendix F4 (we compute the state of the resource )
        popof.idem[i].envt=pop.idem[i].envt*(exp(para.r))/(pop.idem[i].envt*(exp(para.r)-1)+exp(sumconso*para.tau)); // determines the state of the environment at the beginning of the next generation
    }

    //We form the offspring generation (step (ii) of the procedure described in Appendix F4 )
    vector <double> p; //Declare a vector for the relative probability for each individual to have sucessful offspring in the first deme
    // Philopatric elements (for individuals in the first deme)
    p.push_back(fec[0][0]*(1-pop.idem[0].indiv[0].d)); 
    //Add to the vector the weighted fecundity of the first individual
    for (int j = 1; j < para.n; j++){
    //Add to the vector the (cumulated) weighted fecundity of other individuals
        p.push_back(p.back()+fec[0][j]*(1-pop.idem[0].indiv[j].d));}

    //Parapatric elements (for individuals outstide the first deme)
    for (int i=1; i < para.nd; i++){
        for (int j = 0; j < para.n; j++){
            p.push_back(p.back()+fec[i][j]*(1-para.cd)*pop.idem[i].indiv[j].d/(para.nd-1));}
    //Add the (cumulated) weighted fecundity of remaining individuals
    }

    for (int i = 0; i < para.nd; i++){ // Filling the breeding spots of the deme considered
        for (int j = 0; j < para.n; j++){
            // Sampling from the parental generation the next individual in each breeding spot 
            int k=-1;
            int index_ind=-1;
            int index_deme=0;
			double r=eng.rand( p.back() );
            do{
                k++;
                if (index_ind==(para.n-1)){
                    index_deme++;
                    index_ind=0;}
                else{
                    index_ind++;}
            }
            while (p[k]<=r);
        //Once the parent has been sampled, it duplicates in the offspring population with some errors (mutations; step (iii) of the procedure described in Appendix F4)
        popof.idem[i].indiv[j]=mutation(para, pop.idem[index_deme].indiv[index_ind]); 
        }
        //For the next deme considered, update the vector of cumulated weighted fecundity 
        if (i<(para.nd-1))
        {
            for (int i2=i; i2<para.nd; i2++){//for each deme starting with the current one...
                for (int j=0; j<para.n;j++){
                    if (i2==i+1)//...update the vector so that the fecundity of the individuals in the next deme are weighted by the appropriate value
                    { 
                        p[i2*para.n+j]=p[i2*para.n+j-1]+fec[i2][j]*(1-pop.idem[i2].indiv[j].d);
                    }
                    else{//...and update the (cumulated) weighted fecundity of every other elements in the population
                        p[i2*para.n+j]=p[i2*para.n+j-1]+fec[i2][j]*pop.idem[i2].indiv[j].d*(1-para.cd)/(para.nd-1);
                    }
                }
            }
        }
    }
    return popof;
}


void generation(parameters para){
// Create a file to keep track of the mean and variance in the attack rate, dispersal and resource abundance in the population
   char nomFichier[200];
	stringstream nom;
	nom <<"C:\\Users\\PATH\\variables.txt";
	nom >> nomFichier;
   ofstream myfile(nomFichier);
//Initiate the population
    population pop=initial(para);
//Compute and write in file the mean and variance in the attack rate, dispersal and resource abundance for the initial population
    double meanz=para.iniz, varz=0, meand=para.inid, vard=0, meane=pop.idem[1].envt, vare=0;
    myfile <<0<<" "<<meanz<<" "<<varz<<" "<<meand<<" "<<vard<<" "<<meane<<" "<<vare<<endl;
    cout<<"mean in z="<<meanz<<" var in z="<<varz<<" mean in e="<<meane<<" var in e="<<vare<<endl;   
    for (int it = 0; it < para.gen; it++){   
        // iterate the life-cycle for para.gen generations
        cout<<it<<endl;
        pop=iteration(para,pop);
       meanz=0,varz=0, meand=0,vard=0, meane=0,vare=0;
       //Compute and write in file the mean and variance in the attack rate, dispersal and resource abundance for the population each generation
        for (int i = 0; i < para.nd; i++){
            for (int j = 0; j < para.n ; j++){
                meanz+=pop.idem[i].indiv[j].z/(para.nd*para.n);
                varz+=pow(pop.idem[i].indiv[j].z,2)/(para.nd*para.n);
                meand+=pop.idem[i].indiv[j].d/(para.nd*para.n);
                vard+=pow(pop.idem[i].indiv[j].d,2)/(para.nd*para.n);
            }
            meane+=pop.idem[i].envt/para.nd;
            vare+=pow(pop.idem[i].envt,2)/para.nd;
        }
        vard-=pow(meand,2);
        varz-=pow(meanz,2);
        vare-=pow(meane,2);
        myfile <<it+1<<" "<<meanz<<" "<<varz<<" "<<meand<<" "<<vard<<" "<<meane<<" "<<vare<<endl;
        if ((it+1)%para.sample==0)
                {//Each para.sample generation, sample the entire population and write in files the dispersal, attack rate of every individual, and resource abundance in each deme
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
                
                char nomFichiere[200];
	            stringstream nome;
                nome <<"C:\\Users\\PATH\\sampling_e_generation_"<<it+1<<".txt";
	            nome >> nomFichiere;
                ofstream myfilee(nomFichiere);

                for (int i = 0; i < para.nd; i++){
                for (int j = 0; j < para.n ; j++){
                   myfilez<<pop.idem[i].indiv[j].z<<endl;
                   myfiled<<pop.idem[i].indiv[j].d<<endl;
               
                }
                myfilee <<pop.idem[i].envt<<endl;
            }
        }
        
    }  

}