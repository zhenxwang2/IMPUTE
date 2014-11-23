#include "VcfFileReader.h"
#include "StringBasics.h"
#include "StringHash.h"
#include "MemoryAllocators.h"

#include <stdio.h>
#include <limits.h>
#include <math.h>

#include <iostream>     // std::cout, std::endl
#include <iomanip>

using namespace std

// Load Genotypes
void Initialize(char ** HapArray, string file_name){


}

int CountStudyHap(){}
int CountRefHap(){}
int CountMarkers(){}

// Run HMM 
void RunLeftHMM(char * StudyHap, double ** Lmat){
    
    InitialFirstVector(StudyHap, Lmat[0]);
    for (int j=1; j<Ms; j++) {
        Transpose(StudyHap[j], Lmat[j-1], Lmat[j], theta);
        CondGL(StudyHap[j], glvec, epsilon);
        for (int i=0; i<Nr; i++) {
            Lmat[j][i] = Lmat[j][i] * glvec[i];
        }
    }

}

RunRightHMM(char * StudyHap, double ** Rmat){}

ConbineHMM(double ** Lmat, double ** Rmat, double ** V){}


// HMM Functions

void InitialFirstVector(int StudyHap, double * Sstart){
    int Nr = Sstart->size();
    for (int i = 0; i < Nr; i++) {
        Sstart[i] = 1.0/Nr;
    }
    return;
}

void  Transpose(int StudyHap, double * Sfrom, double * Sto, double &theta)
{
// For each study haplotype: StudyHap
// Calculate the current state probabilities transitioned from a previous position with probability vector Sfrom
// State space is consisted of N_r reference haplotypes
// Sfrom, Sto is a vector of size N_r, Sfrom[i] is the probability of being in the state of reference haplotype S_i, i= 1, ..., N_r, on previous marker position
// theta is transition rate, i.e., recombination rate, assume it is the same for all markers
    int Nr = Sfrom->size(); // total number of reference haplotypes
    double p_sum; // sum of vector Sfrom

	for(int i = 0; i < Nr; i++){
        p_sum = 0.0;
        for (int j=0; j < Nr; j++) {
                p_sum += Sfrom[j];
        }
        
        Sto[i] = (1.0 - theta) * Sfrom[i] + theat / (double)N_r * p_sum;
	}
    
    return;
}

void CondGL(int StudyHap, double * GV, double &espilon){

    int Nr = GV->size();
    char *study_allel = GetStudyHap(StudyHap);

    for (int i=0; i<Nr; i++) {
        char *ref_allel = GetRefHap(RefHap[i]);
        if (strcmp(study_allel, ref_allel) == 0) {
            GV[i] = 1 - epsilon;
        }
        else GV[i] = epsilon;
    }
    return;
}



// Imputation

void Impute(char * StudyHap, double ** V){}


