void Haplotyper::ConditionOnData(float * matrix, int marker, char genotype)
   {
   // We treat missing genotypes as uninformative about the mosaic's
   // underlying state. If we were to allow for deletions and the like,
   // that may no longer be true.
   if (genotype == GENOTYPE_MISSING)
      return;
 
   for (int i = 0; i < states; i++)
      {
      double factors[2];
 
      factors[0] = Penetrance(marker, haplotypes[i][marker], genotype);
      factors[1] = Penetrance(marker, haplotypes[i][marker] + 1, genotype);
 
      for (int j = 0; j <= i; j++, matrix++)
         *matrix *= factors[haplotypes[j][marker]];
      }
   }
