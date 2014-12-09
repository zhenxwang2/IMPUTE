#include "impute.h"
#include "CalcGL.h"

//marker/haplotype index starting from 0
//char GetStudyHaplotype(int haplotype, int marker)

//char GetReferenceHaplotype(int haplotype, int marker);
//
char ** samplematrix;
char ** referencematrix;

int main(int argc, char ** argv){
  
  if(argc < 4){cout << "command format: ./impute sample.vcf reference.vcf output.vcf" << endl;}
  
  char *samplevcf = argv[1];
  char *referencevcf = argv[2];
  string outputvcf = argv[3];
  
  int number_markers = countMarkers(samplevcf);
  int num_sample_haplotypes = countHaplotypes(samplevcf);
  
    // Load study sample haplotype data
  samplematrix = AllocateCharMatrix(num_sample_haplotypes, number_markers);
  loadStudyGenotypes(samplematrix, num_sample_haplotypes, number_markers, samplevcf);
  
  //to be added
  // loadReferenceGenotypes as variable haplotypes
  int num_ref_markers = countMarkers(referencevcf);
  int num_reference_haplotypes = countHaplotypes(referencevcf);
  referencematrix = AllocateCharMatrix(num_reference_haplotypes, num_ref_markers);
  loadReferenceGenotypes(referenceMatrix, num_reference_haplotypes, num_ref_markers, referencevcf);
  
    int num_ref_haplotypes = countHaplotypes(referencevcf); //number of reference haplotypes
    
    //allocate memory for HMM
    // Define double array to save GenotypeLikelyhood and Frequences
    double ** matrix= AllocateDoubleMatrix(num_ref_haplotypes,number_markers);
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

