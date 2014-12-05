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
void RunLeftHMM(int StudyHap, double ** Lmat){
    
    InitialFirstVector(StudyHap, Lmat[0]);
    if (GetStudyHap(StudyHap, 0) != NULL) {
        CondGL(StudyHap[0], Lmat[0], epsilon);
    }
    
    for (int j=1; j<Ms; j++) {
        if (GetStudyHap(StudyHap, j) != NULL) {
            Transpose(StudyHap, Lmat[j-1], Lmat[j], theta);
            CondGL(StudyHap, j, Lmat[j], epsilon, freq);
        }
    }

}

RunRightHMM(int StudyHap, double ** Rmat){}

ConbineHMM(double ** Lmat, double ** Rmat, double ** V){}


// HMM Functions
void  Transpose(int StudyHap, double * Sfrom, double * Sto, double &theta)
{
// For each study haplotype: StudyHap
// Calculate the current state probabilities transitioned from a previous position with probability vector Sfrom
// State space is consisted of N_r reference haplotypes
// Sfrom, Sto is a vector of size N_r, Sfrom[i] is the probability of being in the state of reference haplotype S_i, i= 1, ..., N_r, on previous marker position
// theta is transition rate, i.e., recombination rate, assume it is the same for all markers
    int Nr = Sfrom->size(); // total number of reference haplotypes
    double p_sum; // sum of vector Sfrom
    
    if (theta == 0) {
        for (int i=0; i < Nr; i++) {
            Sto[i] = Sfrom[i];
        }
    }

    else{
        p_sum = 0.0;
        for (int j=0; j < Nr; j++) {
            p_sum += Sfrom[j];
        }
        p_sum *= theta / (double)N_r;
        
        double q = 1.0 - theta;
        
        if (p_sum < 1e-10) {
            p_sum *= 1e15;
            q *= 1e15;
        }
        
        for(int i = 0; i < Nr; i++){
            Sto[i] = q * Sfrom[i] + p_sum;
        }
    }
    
    return;
    
}

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
    
    if (theta == 0) {
        for (int i=0; i < Nr; i++) {
            Sto[i] = Sfrom[i];
        }
    }

    else{
        p_sum = 0.0;
        for (int j=0; j < Nr; j++) {
            p_sum += Sfrom[j];
        }
        p_sum *= theta / (double)N_r;
        
        double q = 1.0 - theta;
        
        if (p_sum < 1e-10) {
            p_sum *= 1e15;
            q *= 1e15;
        }
        
        for(int i = 0; i < Nr; i++){
            Sto[i] = q * Sfrom[i] + p_sum;
        }
    }
    
    return;
    
}

void CondGL(int StudyHap, int position, double * GV, double &epsilon, double &freq){

    int Nr = GV->size();
    char *study_allel = GetStudyHap(StudyHap, position);

    double prand = epsilon * freq;
    double pmatch = (1.0 - epsilon) + prand;
    
    for (int i=0; i<Nr; i++) {
        char *ref_allel = GetRefHap(RefHap[i], position);
        if (strcmp(study_allel, ref_allel) == 0) {
            GV[i] *= pmatch;
        }
        else GV[i] *= prand;
    }
    return;
}



// Imputation

void Impute(char * StudyHap, double ** V){}


