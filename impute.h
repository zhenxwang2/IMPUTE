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

// Run HMM 

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

// Impute


