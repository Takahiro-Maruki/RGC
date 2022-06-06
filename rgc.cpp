// Updated on 05/24/22

/*

Program rgc.cpp to refine GATK genotype calls by inferring the ancestral genotype and filtering
individual genotypes at each SNP site.  SNP sites are excluded if at least two MA lines have
mutations.  The inferred ancestral genotype can be a reference homozygote and the inferred
individual genotype can be an alternative homozygote.  Only the MA line with mutation is allowed
to contain >= 2 mutant nucleotide reads.  Major and minor nucleotide reads are used for calculating
the likelihood of the homozygous and heteroygous genotypes, respectively.  The alternative allele
depth is examined as long as the genotype call is available, regardless of the depth.  Only sites
with genotype calls in all MA lines with mean depth at least six are printed out.  The two-tailed
binomial test is carried out using the cumulative binomial function to examine the goodness of fit
between the binomial expectation and obserbed nucleotide read count.

Input: Table of hard-filtered GATK SNP genotype calls.

*/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
using namespace std;

// Function declaration
double binomial(int n, int x, double p);

// Function returning a binomial probability
double binomial(int n, int x, double p) {
	// Local variable declarations
	double ln_p;	// natural logarithm of the binomial probability
	double prob;	// binomial probability
	int ag, bg, cg;	// counters

	ln_p = 0.0;
	for (ag=n; ag>=1; ag--) {
		ln_p = ln_p + log(ag);
	}
	for (bg=x; bg>=1; bg--) {
                ln_p = ln_p - log(bg);
        }
	for ( cg=(n-x); cg>=1; cg-- ) {
                ln_p = ln_p - log(cg);
        }
	ln_p = ln_p + x*log(p) + (n - x)*log( (double)1.0 - p );
	prob = exp(ln_p);
	return prob;
}

int main(int argc, char *argv[])
{
	// Default values of the options
	const char* in_file = {"In_rgc.txt"};
	const char* out_file = {"Out_rgc.txt"};
	int min_dp = 1;	// minimum required total coverage
	int max_dp = 2000000000; // maximum allowed total coverage
	int mc = 6;	// minimum required individual coverage
	double cv_het = 0.025;	// critical value for heteozygous cumulative binomial probability
	double cv_hom = 0.05;  // critical value for the homozygous cumulative binomial probability
	int print_help = 0;
	int argz = 1; // argument counter

	// Read specified settings
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-in") == 0) {
			in_file = argv[++argz];
		} else if (strcmp(argv[argz], "-min_dp") == 0) {
                        sscanf(argv[++argz], "%d", &min_dp);
                } else if (strcmp(argv[argz], "-max_dp") == 0) {
                        sscanf(argv[++argz], "%d", &max_dp);
		} else if (strcmp(argv[argz], "-mc") == 0) {
			sscanf(argv[++argz], "%d", &mc);
		} else if (strcmp(argv[argz], "-cv_het") == 0) {
			sscanf(argv[++argz], "%lf", &cv_het);
		} else if (strcmp(argv[argz], "-cv_hom") == 0) {
                        sscanf(argv[++argz], "%lf", &cv_hom);
		} else if (strcmp(argv[argz], "-out") == 0) {
			out_file = argv[++argz];
		} else {
			fprintf(stderr, "unknown option %s\n", argv[argz]);
			print_help = 1;
			break;
		}
		argz++;
	}
	if (print_help) { // print error/usage message ?
		fprintf(stderr, "USAGE: %s {<options>}\n", argv[0]);
		fprintf(stderr, "	options:\n");
		fprintf(stderr, "	-h: print the usage message\n");
		fprintf(stderr, "	-in <s>: specify the input file\n");
		fprintf(stderr, "       -min_dp <d>: specify the minimum DP\n");
                fprintf(stderr, "       -max_dp <d>: specify the maximum DP\n");
		fprintf(stderr, "       -mc <d>: specify the minimum required coverage\n");
		fprintf(stderr, "       -cv_het <f>: specify the critical value for the heterozygous cumulative binomial probability\n");
		fprintf(stderr, "       -cv_hom <f>: specify the critical value for the homozygous cumulative binomial probability\n");
		fprintf(stderr, "       -out <s>: specify the name of the output file\n");
		exit(1);
	}

	FILE *outstream;

	// Open the output file
	outstream = fopen(out_file, "w");
	if (outstream == NULL ) { // Exit on failure
		fprintf(stderr, "Cannot open %s for writing.\n", out_file);
		exit(1);
	}

	string line; // String buffer for the input file

	ifstream inputFile(in_file); // Try to open the input file
	if ( !inputFile.is_open() ) { // Exit on failure
		fprintf(stderr, "Cannot open %s for reading.\n", in_file);
		exit(1);
	}

	vector <string> ma_line;	// MA line ID

	// Clear the vector
	ma_line.clear();

	// Read the header
	string h_chrom, h_pos, h_ref, h_alt, h_qual, h_filter, h_ac, h_an, h_dp;
	int dotp;	// dot position
	string id;
	int lg;         // line counter
	getline(inputFile, line);
	istringstream ss(line);
	ss >> h_chrom >> h_pos >> h_ref >> h_alt >> h_qual >> h_filter >> h_ac >> h_an >> h_dp;
	string str; // Temporarily stores line-specific labels
	while (true) {
		ss >> str;
		dotp = str.find(".");
		id = str.substr(0, dotp);
		if ( find(ma_line.begin(), ma_line.end(), id) == ma_line.end() ) {
			ma_line.push_back(id);
		}
		if ( ss.eof() ) {
			break;
		}
	}
	int num_lines = ma_line.size();
	printf("%d lines to be analyzed in %s\n", num_lines, in_file);

	// Print out the header of output file
	fprintf(outstream, "CHROM\tPOS\tREF\tALT\tANC\tQUAL\tAC\tAN\tDP\t");
	// printf("CHROM\tPOS\tREF\tALT\tANC\tQUAL\tAC\tAN\tDP\t");
	for (lg=0; lg<num_lines-1; lg++) {
		fprintf(outstream, "%s\t", ma_line.at(lg).c_str());
		// printf("%s\t", ma_line.at(lg).c_str());
	}
	fprintf(outstream, "%s\n", ma_line.at(lg).c_str());
	// printf("%s\n", ma_line.at(lg).c_str());
	// For testing
	// fprintf(outstream, "CHROM\tPOS\tAC\n");
	// printf("CHROM\tPOS\tAC\n");

	// Read the input data
	string chrom, ref, alt, filter;
	int pos, ac, an, dp;
	int num_alt;	// number of alternative alleles
	string t_gt, ad, gq;
	int cov;
	float qual;
	int tot_cov_ref, tot_cov_alt;	// total reference and alternative coverage
	int cp;		// comma position
	string s_ad1, s_ad2;
	int * ad1;	// pointer to an array containing allele depth 1
	int * ad2;      // pointer to an array containing allele depth 2
	int dg;	// depth counter
	string a1, a2;	// alleles in a genotype
	string * gt;	// pointer to ar array containing genotype calls
	double prob_anc_het, prob_anc_hom;	// cumulative probabilities of heterozygous and homozygous ancestral genotypes
	string anc_geno;	// inferred ancestral genotype
	double prob_ind_het, prob_ind_hom;      // cumulative probabilities of heterozygous and homozygous individual genotypes
	int f_ac, f_an;	// filtered allele count and allele number
	int f_tot_cov_ref, f_tot_cov_alt;   // filtered total reference and alternative coverage
	int num_mut_lines;	// number of MA lines having mutations
	int num_altr_alt;       // number of MA lines with at least two alternative nucleotide reads
	int num_genos;		// number of MA lines with genotype calls
	int mean_dp;	// mean depth over MA lines
	double base_prob_ind_hom;	// baseline probability of homozygous individual genotype

	gt = (string *)malloc( (num_lines+1)*sizeof(string) );
        ad1 = (int *)malloc( (num_lines+1)*sizeof(int) );
        ad2 = (int *)malloc( (num_lines+1)*sizeof(int) );

	while ( getline(inputFile, line) ) {
		istringstream ss(line);
		ss >> chrom >> pos >> ref >> alt >> qual >> filter >> ac >> an >> dp;
		num_alt = alt.length();
		if (num_alt == 1 && filter == "PASS" && dp >= min_dp && dp <= max_dp) {
			tot_cov_ref = 0;
			tot_cov_alt = 0;
			num_altr_alt = 0;
			num_genos = 0;
			for (lg = 1; lg <= num_lines; lg++) {	// Read individual data
				ss >> t_gt >> ad >> cov >> gq;
				// For debugging
				/*
				if (chrom == "scaffold_11" && pos == 2737668) {
					printf("%s: %d\n", ma_line.at(lg-1).c_str(), cov);
				}
				*/
				if (cov >= mc) {
					if (t_gt != "./.") {	// genotype call available
						gt[lg] = t_gt;
						num_genos = num_genos + 1;
						cp = ad.find(",");
						s_ad1 = ad.substr(0, cp);
						s_ad2 = ad.substr(cp+1);
						ad1[lg] = atoi(s_ad1.c_str());
						ad2[lg] = atoi(s_ad2.c_str());
						tot_cov_ref = tot_cov_ref + ad1[lg];
						tot_cov_alt = tot_cov_alt + ad2[lg];
						if (ad2[lg] >= 2) {
							num_altr_alt = num_altr_alt + 1;
						}
						// printf("Depth of allele %s: %d\tDepth of allele %s: %d\n", a1.c_str(), ad1, a2.c_str(), ad2);
					} else {
						gt[lg] = "NA";
						ad1[lg] = 0;
						ad2[lg] = 0;
					}
				} else {
					if (t_gt != "./.") {    // genotype call available
						gt[lg] = "NA";
						num_genos = num_genos + 1;
						cp = ad.find(",");
                                        	s_ad1 = ad.substr(0, cp);
                                        	s_ad2 = ad.substr(cp+1);
                                        	ad1[lg] = atoi(s_ad1.c_str());
                                        	ad2[lg] = atoi(s_ad2.c_str());
                                        	tot_cov_ref = tot_cov_ref + ad1[lg];
                                        	tot_cov_alt = tot_cov_alt + ad2[lg];
						if (ad2[lg] >= 2) {
                                                        num_altr_alt = num_altr_alt + 1;
                                                }
					} else {
						gt[lg] = "NA";
                                                ad1[lg] = 0;
                                                ad2[lg] = 0;
					}
				}
			}
			mean_dp = (double)dp/(double)num_genos;
			if (num_genos == num_lines && mean_dp > (double)6.0) {	// genotype calls in all MA lines with mean depth at least six
				// Infer the ancetral genotype
				if (tot_cov_alt <= tot_cov_ref) {	// alternative allele is the minor allele
					prob_anc_het = binomial(tot_cov_ref+tot_cov_alt, tot_cov_alt, 0.5);
				} else {
                                        prob_anc_het = binomial(tot_cov_ref+tot_cov_alt, tot_cov_ref, 0.5);
				}
				if (tot_cov_ref >= tot_cov_alt) {	// reference allele is the major allele
                                	prob_anc_hom = binomial(tot_cov_ref+tot_cov_alt, tot_cov_ref, 0.99);
				} else {
					prob_anc_hom = binomial(tot_cov_ref+tot_cov_alt, tot_cov_alt, 0.99);
				}
				if (prob_anc_het > prob_anc_hom) {	// ancestral genotype likely to be heterozygous
					anc_geno = ref + alt;
					// Filter individual genotypes
					f_ac = 0;
					f_an = 0;
					f_tot_cov_ref = tot_cov_ref;
					f_tot_cov_alt = tot_cov_alt;
					num_mut_lines = 0;
					for (lg = 1; lg <= num_lines; lg++) {
						// For debugging
						/*
						if (chrom == "scaffold_11" && pos == 2737668) {
							printf("%s: %s\n", ma_line.at(lg-1).c_str(), gt[lg].c_str());
						}
						*/
						if (gt[lg] != "NA") {	// genotype call available
							f_an = f_an + 2;
							prob_ind_het = 0.0;
							if (ad2[lg] < ad1[lg]) {	// alternative nucleotide is less frequent
								for (dg=ad2[lg]; dg>=0; dg--) {
									prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], dg, 0.5);
								}
								prob_ind_het = (double)2.0*prob_ind_het;
							} else if (ad2[lg] == ad1[lg]) {	// nucleotide reads equally abundant
								for (dg=ad2[lg]-1; dg>=0; dg--) {
									prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], dg, 0.5);
								}
								prob_ind_het = (double)2.0*prob_ind_het;
								prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], ad2[lg], 0.5);
							} else {	// reference nucleotide is less frequent
								for (dg=ad1[lg]; dg>=0; dg--) {
									prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], dg, 0.5);
								}
								prob_ind_het = (double)2.0*prob_ind_het;
							}
							if (prob_ind_het < cv_het) {	// unlikely to be heterozygous
								if (ad1[lg] == 0) {
									gt[lg] = alt + alt;
									f_ac = f_ac + 2;
									f_tot_cov_alt = f_tot_cov_alt - ad2[lg];
									num_mut_lines = num_mut_lines + 1;
								} else if (ad2[lg] == 0) {
									gt[lg] = ref + ref;
									f_tot_cov_ref = f_tot_cov_ref - ad1[lg];
									num_mut_lines = num_mut_lines + 1;
								} else {
									gt[lg] = anc_geno;
									f_ac = f_ac + 1;
								}
							} else {	// likely to be heterozygous
								gt[lg] = anc_geno;
								f_ac = f_ac + 1;
							}
						}
					}
					// Examine whether the ancestral genotype keeps appearing heterozygous after removing reads from MA lines containing mutation
					if (f_tot_cov_alt <= f_tot_cov_ref) {	// alternative allele is the minor allele
						prob_anc_het = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_alt, 0.5);
					} else {
						prob_anc_het = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_ref, 0.5);
					}
					if (f_tot_cov_ref >= f_tot_cov_alt) {	// reference allele is the major allele
						prob_anc_hom = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_ref, 0.99);
					} else {
						prob_anc_hom = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_alt, 0.99);
					}
					if (prob_anc_het > prob_anc_hom && f_ac > 0 && num_mut_lines == 1) {	// Keep the SNP site only if the ancestral genotype keeps appearing heterozygous, the site remains polymorphic, and only one MA line has mutations
						fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%f\t%d\t%d\t%d\t", chrom.c_str(), pos, ref.c_str(), alt.c_str(), anc_geno.c_str(), qual, f_ac, f_an, dp);
						// printf("%s\t%d\t%s\t%s\t%s\t%f\t%d\t%d\t%d\t", chrom.c_str(), pos, ref.c_str(), alt.c_str(), anc_geno.c_str(), qual, f_ac, f_an, dp);
						for (lg = 1; lg < num_lines; lg++) {
							fprintf(outstream, "%s\t", gt[lg].c_str());
							//printf(gt[lg].c_str());
						}
						fprintf(outstream, "%s\n", gt[lg].c_str());
						// printf("%s\n", gt[lg].c_str());
					}
				} else {	// ancestral genotype likely to be major homozygous
					if (tot_cov_ref >= tot_cov_alt) {       // reference allele is the major allele
						anc_geno = ref + ref;
					} else {
						anc_geno = alt + alt;
					}
					// Filter individual genotypes
					f_ac = 0;
                        	        f_an = 0;
                        	        f_tot_cov_ref = tot_cov_ref;
                        	        f_tot_cov_alt = tot_cov_alt;
					num_mut_lines = 0;
                        	        for (lg = 1; lg <= num_lines; lg++) {
						if (gt[lg] != "NA") {   // genotype call available
							f_an = f_an + 2;
							prob_ind_hom = 0.0;
							if (tot_cov_ref >= tot_cov_alt) {	// reference allele is the major allele
								base_prob_ind_hom = binomial(ad1[lg]+ad2[lg], ad1[lg], 0.99);
								prob_ind_hom = base_prob_ind_hom;
								for (dg=ad1[lg]+ad2[lg]; dg>=0; dg--) {
									if (binomial(ad1[lg]+ad2[lg], dg, 0.99) < base_prob_ind_hom) {
										prob_ind_hom = prob_ind_hom + binomial(ad1[lg]+ad2[lg], dg, 0.99);
									}
								}
							} else {
								base_prob_ind_hom = binomial(ad1[lg]+ad2[lg], ad2[lg], 0.99);
								for (dg=ad1[lg]+ad2[lg]; dg>=0; dg--) {
									if (binomial(ad1[lg]+ad2[lg], dg, 0.99) < base_prob_ind_hom) {
										prob_ind_hom = prob_ind_hom + binomial(ad1[lg]+ad2[lg], dg, 0.99);
									}
								}
							}
							if (prob_ind_hom < cv_hom) {    // unlikely to be major homozygous
								prob_ind_het = 0.0;
								if (ad2[lg] < ad1[lg]) {	// alternative nucleotide is less abundant
									for (dg=ad2[lg]; dg>=0; dg--) {
										prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], dg, 0.5);
									}
									prob_ind_het = (double)2.0*prob_ind_het;
								} else if (ad2[lg] == ad1[lg]) {	// nucleotide reads equally abundant
									for (dg=ad2[lg]-1; dg>=0; dg--) {
										prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], dg, 0.5);
									}
									prob_ind_het = (double)2.0*prob_ind_het;
									prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], ad2[lg], 0.5);
								} else {	// reference nucleotide is less abundant
									for (dg=ad1[lg]; dg>=0; dg--) {
										prob_ind_het = prob_ind_het + binomial(ad1[lg]+ad2[lg], dg, 0.5);
									}
									prob_ind_het = (double)2.0*prob_ind_het;
								}
								if (prob_ind_het > cv_het && ad1[lg] >= 2 && ad2[lg] >= 2) {	// likely to be heterozygous
									gt[lg] = ref + alt;
									f_ac = f_ac + 1;
									f_tot_cov_ref = f_tot_cov_ref - ad1[lg];
									f_tot_cov_alt = f_tot_cov_alt - ad2[lg];
									num_mut_lines = num_mut_lines + 1;
								} else {
									if (ad1[lg] == 0) {
										gt[lg] = alt + alt;
										f_ac = f_ac + 2;
                                                        	       	 	f_tot_cov_alt = f_tot_cov_alt - ad2[lg];
                                                        	        	num_mut_lines = num_mut_lines + 1;
									} else {
										gt[lg] = anc_geno;
									}
								}
							} else {        // likely to be homozygous
								gt[lg] = anc_geno;
							}
						}
					}
					// Examine whether the ancestral genotype keeps appearing homozygous after removing reads from MA lines containing mutation
					if (f_tot_cov_alt <= f_tot_cov_ref) {	// alternative allele is the minor allele
                                        	prob_anc_het = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_alt, 0.5);
					} else {
						prob_anc_het = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_ref, 0.5);
					}
					if (f_tot_cov_ref >= f_tot_cov_alt) {	// reference allele is the major allele
                                        	prob_anc_hom = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_ref, 0.99);
					} else {
						prob_anc_hom = binomial(f_tot_cov_ref+f_tot_cov_alt, f_tot_cov_alt, 0.99);
					}
					if (prob_anc_het <= prob_anc_hom && f_ac > 0 && num_mut_lines == 1 && num_altr_alt == 1) {      // Keep the SNP site only if the ancestral genotype keeps appearing homozygous, the site remains polymorphic, and only one MA line has mutations with at least two mutant nucleotide reads
						fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%f\t%d\t%d\t%d\t", chrom.c_str(), pos, ref.c_str(), alt.c_str(), anc_geno.c_str(), qual, f_ac, f_an, dp);
						// printf("%s\t%d\t%s\t%s\t%s\t%f\t%d\t%d\t%d\t", chrom.c_str(), pos, ref.c_str(), alt.c_str(), anc_geno.c_str(), qual, f_ac, f_an, dp);
						for (lg = 1; lg < num_lines; lg++) {
							fprintf(outstream, "%s\t", gt[lg].c_str());
							//printf(gt[lg].c_str());
						}
						fprintf(outstream, "%s\n", gt[lg].c_str());
						// printf("%s\n", gt[lg].c_str());
					}
				}
				// For testing
				// fprintf(outstream, "%s\t%d\t%s\n", chrom.c_str(), pos, anc_geno.c_str());
				// printf("%s\t%d\t%s\n", chrom.c_str(), pos, anc_geno.c_str());
			}
		}
	}

	return 0;
}
