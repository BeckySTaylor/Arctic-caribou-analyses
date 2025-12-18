##Script to run the DoS calculator script (which does not need changing). Enter the path and names of the non-synonymous and synonymous VCF files, respectively.
##The give the name of a file containing a comma-separated list of individuals to be included ('-p', example included as another file), the outgroup individuals for '-o', the name of the output file, and the minumum number for Pn + Dn. 

./dos_calculator.sh \
  -n ../Apr2025_96indv_Arctic_NM1_NALout_autosomes_Reindeer_5missing_nonSyn_biallelic.recode.vcf.gz \
  -s ../Apr2025_96indv_Arctic_NM1_NALout_autosomes_Reindeer_5missing_Syn_biallelic.recode.vcf.gz \
  -p All_Inds.txt \
  -o 22832,27689,27694,45932,45933 \
  -O dos_results_All_grouped_Nalout_5missing.tsv \
  -m 4
