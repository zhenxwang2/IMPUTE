#include "impute.h"
#include "transpose.h"
#include "conditionondata.h"

//marker/haplotype index starting from 0
//char GetStudyHaplotype(int haplotype, int marker)

//char GetReferenceHaplotype(int haplotype, int marker);
//

int main(int argc, char ** argv){
  
  if(argc < 4){cout << "command format: impute.cpp sample.vcf reference.vcf output.vcf" << endl;}
  
  char *samplevcf = argv[1];
  char *referencevcf = argv[2];
  string outputvcf = argv[3]
  
  int number_markers = countMarkers(samplevcf);
  int num_sample_haplotypes = countStudyHaplotypes(samplevcf);
  
  char ** samplematrix = AllocateCharMatrix(num_sample_haplotypes, number_markers);
  loadStudyGenotypes(samplematrix, num_sample_haplotypes, number_markers, samplevcf);
  
  //to be added
  // loadReferenceGenotypes
  double ** matrix=
  double ** freqs=
  double ** probmatrix = AllocateCharMatrix(num_sample_haplotypes, number_markers)
   RunLeftHmm(samplematrix, referencematrix, freqs )
   RunRightHmmCombine(samplematrix,referencematrix, freqs)
  //allocate memory for HMM
  //walk forward to produce the prob matrix
  //walk reverse to calculate backward prob matrix, and overwrite matrix to combined probs.


  return 1;
}



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
