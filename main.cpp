#include "impute.h"


using namespace std;


int main(int argc, char ** argv){

    string file_ref = "reference.vcf";
    string file_study = "study.vcf";
    
    // Load Genotypes
    int Ns; // number of study haplotypes
    int Nr; // number of reference haplotypes
    int Ms; // number of study markers
    int Mr; // number of reference markers
    //int iter; // # of mcmc iterations
    
    char ** HapStudy[Ns][Ms]; // array to save study haplotypes data
    char ** HapRef[Mr][Nr]; // array to save reference haplotype data
    Initialize(HapStudy, file_study); // load in data
    Initialize(HapRef, file_ref);
    
    double ** V[Nr][Ms]; // used to do imputation
    
    for (int i = 0; i< Ms; i++) {
        char * hap_study = HapStudy[i];
        RunHMM(hap_study, V, iter);
        Impute(hap_study, V);
    }
    //print(hap_study);
    
  return 1;
}
