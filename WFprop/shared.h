/* 

Dr. Christopher J.R. Illingworth, Department of Genetics, University of Cambridge, 2017

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>


using namespace std;

struct haplo {
	string seq1;
	string seq2;
	unsigned int size;
	double fitness;
	int colour;
};

struct nuc {
	int A;
	int C;
	int G;
	int T;
};

struct amap {
	char n1;
	char n2;
};

struct hap_loc {
	unsigned int hap;
	unsigned int loc;
	int value;
};

struct sample {
	int time;
	float freq;
};

struct mutant {
	int time;
	int clone;
	double selection;
	int c_index;
};

struct clone {
	vector<double> individuals;
	vector<double> fitness;
};

struct inpair {
	string seq1;
	string seq2;
	int clone1;
	int clone2;
};

struct dat_store {
	int position;
	vector<double> freqs;
};

struct trajectory {
	int start;
	int finish;
	int position;
	int open;
	int length;
	float max_freq;
	int traj_id;
	vector<float> freq;
	vector<float> sig;
	vector<unsigned int> overlaps;
	int uv_index;
	int n_over;
	double logL;
	int d_timer;
	int f_timer;
	int put_end;
	int put_len;
};

struct run_params {
	
	int N; // Population size
	int F; // number of founders
	int L; // Size of genomes
	int tSim; // Number of generations
	int seed; // seed for random number genrator
	int samp_freq; // sample every samp_freq generations
	int depth; // read depth
	int sel_loc;// locus under selection
	int s_out; //Number of individuals from which sequences are reported
	double sigma; // Selection magnitude; at the moment of one locus only
	double rho; // recombination rate
	double mu; // mutation rate
	string file; // .arp file name
};

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

void GetParameters (run_params& p, int argc, const char **argv);
void ReadData (run_params p, int pop_size, int n_seqs, int s_print, ifstream& arp_file, vector<int>& diff_pos, vector<string>& seqs, vector<haplo>& population, gsl_rng *rgen);
void GetNucleotides (vector<nuc>& nucs, vector<haplo> population);
void GetAlleles (run_params p,vector<amap>& alleles, vector<nuc> nucs);



double FindFit(int sel_loc, char allele, string seq1, string seq2, vector<double> fitness, double h);

void   MakeMutatedHaplotype (int sel_loc, char allele, double h, vector<double> fitness, vector<hap_loc>& mutants, vector<amap> alleles, vector<haplo>& population, gsl_rng *rgen);

void   CalculateRecombination (int L, double rho, int gen_length, int sel_loc, char allele, double h, vector<double> fitness, vector<int> diff_pos, vector<inpair>& parents, vector<haplo>& population, gsl_rng *rgen);

void   ConstructPopulation (int sel_loc, char allele, double h, vector<double> fitness, vector<inpair> parents, vector<haplo>& population);

void   AlleleFreqs (vector<amap> alleles, vector<haplo> population, vector<double>& freqs);


void   GetSelectedAllele (int pop_size, int L, int& sel_loc, char& allele, vector<nuc> nucs, gsl_rng *rgen);
void   GetSelectedAllele2 (int pop_size, int L, int& sel_loc, char& allele, vector<nuc> nucs, gsl_rng *rgen);

void   InitialiseHaploPopulation (int L, int N, vector<haplo>& population, vector<haplo>& pop_temp, vector<haplo>& pop_sample);
void   FindHapMutationLoci (int L, double mu, vector<haplo>& population, vector<hap_loc>& mutants, gsl_rng *rgen);
void   CombineHapSame (int L, vector<haplo>& population);
void   DeleteEmptyHap (vector<haplo>& population, vector<haplo>& pop_temp);
void   CopyHap (vector<haplo>& population, vector<haplo>& pop_temp);
double GlobalFit (vector<haplo>& population);
void   HapSelProbs(double glob_fit, vector<double>& p_select, vector<haplo>& population);
void   HapNextGeneration (int Nt, gsl_rng *rgen, vector<double>& p_select, vector<haplo>& population);
void   DiagPrintPopulation2 (int t_k, int L, ofstream& pop_track, vector<haplo>& population);


void   FreqCalc (int N, int L, gsl_vector *freqs, vector<haplo>& population);
void   PrintFreqs (int L, ofstream& frq_track, gsl_vector *freqs);
void   FreqReset (int L, double sigma, gsl_vector *freqs, gsl_matrix *fitness, vector<haplo>& population,gsl_rng *rgen);

void   PrintInh (int L, ofstream& sig_track, gsl_matrix *fitness);
void   PrintEff (int L, ofstream& sig_track, gsl_matrix *fitness,vector<haplo>& population);
void   ChooseParents (int N, vector<inpair>& parents, vector<haplo> population, gsl_rng *rgen);

