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
    N_r = Sfrom->size(); // total number of reference haplotypes
    double p_sum; // sum of vector Sfrom

	for(int i = 0; i < N_r; i++){
        p_sum = 0.0;
        for (int j=0; j < N_r; j++) {
                p_sum += Sfrom[j];
        }
        
        Sto[i] = (1.0 - theta) * Sfrom[i] + theat / (double)N_r * p_sum;
	}
    
    return;
}
// Impute


