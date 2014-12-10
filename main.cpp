#include "LoadVCF.h"
#include "CalcGL.h"
#include "HMM.h"
#include "Impute.h"



//marker/haplotype index starting from 0

//declared globally
char ** samplematrix;
char ** referencematrix;
    

int main(int argc, char ** argv){
  
  if(argc < 4){cout << "command format: ./impute sample.vcf reference.vcf output.vcf" << endl;}
  
  char *samplevcf = argv[1];
  char *referencevcf = argv[2];
  string outputvcf = argv[3];
    
  
    int number_markers = countMarkers(argv[1]);
    int num_sample_haplotypes = countHaplotypes(argv[1]);
    
    //load sample vcf
    //allocate memory
    samplematrix = AllocateCharMatrix(num_sample_haplotypes, number_markers);
    //load genotypes
    
    loadStudy(samplematrix, num_sample_haplotypes, number_markers, samplevcf);
    int number_reference_markers = countMarkers(argv[2]);
    int num_reference_haplotypes = countHaplotypes(referencevcf);
    
    referencematrix = AllocateCharMatrix(number_reference_markers, num_reference_haplotypes);
    loadReference(referencematrix, num_reference_haplotypes, number_reference_markers, referencevcf);
      
    //allocate memory for HMM
    // Define double array to save GenotypeLikelyhood and Frequences
    double ** matrix= AllocateDoubleMatrix(num_reference_haplotypes,number_markers);
    double ** freqs= AllocateDoubleMatrix(5, number_markers);
    
    
    //walk forward to produce the prob matrix
    //walk reverse to calculate backward prob matrix, and overwrite matrix to combined probs.
    int numIterations; // number of iterations of HMM
    
    for (int i=0; i<numIterations; i++) {
        double ** probmatrix = AllocateCharMatrix(num_sample_haplotypes, number_markers)
        RunLeftHmm(samplematrix, referencematrix, freqs )
        RunRightHmmCombine(samplematrix,referencematrix, freqs)
        
        
        // Impute
    }

  return 1;
}

