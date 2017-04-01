/*

Dr. Christopher J.R. Illingworth, Department of Genetics, University of Cambridge, 2017

*/

#include "shared.h"
#include <iostream>
#include <string>

void GetParameters (run_params& p, int argc, const char **argv) {
	string p_switch;
	int x=1;
	
	p.N=1000;                 // Population size
	p.L=5000;                 // Genome length
	p.tSim=300;               // Number of generations
	p.sigma=0;                // Magnitude of selection
	p.samp_freq=10;           // Sample frequency (generations)
	p.F=20;                   // Number of founder lines
	p.depth=100;              // Sample depth
	p.nseqs=100;              // Number of sequences in .arp file
	p.rho=1e-8;               // Recombination rate
	p.mu=3e-9;                // Mutation rate per base
	p.sel_loc=-1;             // Position of locus under selection
	p.s_out=0;                // Output sequences to file
	p.file="n2000_1_1_1.arp"; // Input file from fastsimcoal
	p.seed=(int) time(NULL);
	
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		cout << p_switch << "\n";
		if (p_switch.compare("--N")==0) {
			x++;
			p.N=atoi(argv[x]);
		}
		else if (p_switch.compare("--F")==0) {
			x++;
			p.F=atoi(argv[x]);
		}
		else if (p_switch.compare("--L")==0) {
			x++;
			p.L=atoi(argv[x]);
		}
		else if (p_switch.compare("--tSim")==0) {
			x++;
			p.tSim=atoi(argv[x]);
		}
		else if (p_switch.compare("--depth")==0) {
			x++;
			p.depth=atoi(argv[x]);
		}
		else if (p_switch.compare("--nseqs")==0) {
			x++;
			p.nseqs=atoi(argv[x]);
		}
		else if (p_switch.compare("--s_out")==0) {
			x++;
			p.s_out=atoi(argv[x]);
		}
		else if (p_switch.compare("--seed")==0) {
			x++;
			p.seed=atoi(argv[x]);
		}
		else if (p_switch.compare("--samp_freq")==0) {
			x++;
			p.samp_freq=atoi(argv[x]);
		}
		else if (p_switch.compare("--sigma")==0) {
			x++;
			p.sigma=atof(argv[x]);
		}
		else if (p_switch.compare("--sel_loc")==0) {
			x++;
			p.sel_loc=atoi(argv[x]);
		}
		else if (p_switch.compare("--rho")==0) {
			x++;
			p.rho=atof(argv[x]);
		}
		else if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		}
		else if (p_switch.compare("--file")==0) {
			x++;
			p.file=argv[x];
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void ReadData (run_params p, int pop_size, int n_seqs, int s_print, ifstream& arp_file, vector<int>& diff_pos, vector<string>& seqs, vector<haplo>& population, gsl_rng *rgen) {
	cout << "ReadData\n";
	int d;
	int diffs;
	vector<string> all_seqs;  //Read all of the fastsimcoal output to here, then select random individuals
	int g=n_seqs/2;
	int a[g];
	for (int i=0;i<g;i++) {
		a[i]=i;
	}
	gsl_ran_shuffle(rgen,a,g,sizeof(int)); //Randomly order previous generation

	
	while (!arp_file.eof()) {
		string ref;
		arp_file >> ref;
	
		if (ref.compare("#Number")==0) {
			for (int i=1;i<=9;i++) {
				arp_file >> ref;
			}
			arp_file >> diffs;
			
			for (int i=1;i<=6;i++) {
				arp_file >> ref;
			}
			arp_file >> d;
			
			for (int i=1;i<=diffs;i++) {
				arp_file >> ref;
				if (i==1) {
					diff_pos.push_back(d);
				} else {
					diff_pos.push_back(atoi(ref.c_str()));
				}
				
			}
			for (int i=1;i<6;i++) {
				arp_file >> ref;
			}
			for (int i=0;i<n_seqs;i++) {
			
				arp_file >> ref;
				arp_file >> ref;
				arp_file >> ref;
				all_seqs.push_back(ref);
			}
			
			ofstream seq_file;

			
			seq_file.open("Sequence_individuals.dat");
			
			unsigned int tot=0;
			
			for (int i=0;i<p.F;i++) {  //Generate a random selection of sequences
			
				haplo hap;
				hap.seq1=all_seqs[2*a[i]];
				hap.seq2=all_seqs[2*a[i]+1];
				if (i<s_print) {
			
					seq_file << i << " 1 " << hap.seq1 << "\n";
					seq_file << i << " 2 " << hap.seq2 << "\n";
				}
				hap.size=floor(pop_size/p.F);
				tot=tot+hap.size;
				hap.fitness=0;
				population.push_back(hap);
			}
			
			cout << "Total " << tot << "\n";
			
		
			break;
		}
	}
}

void GetNucleotides (vector<nuc>& nucs, vector<haplo> population) {
	cout << "Number of differentiating loci " <<population[0].seq1.size() << "\n";
	for (int i=0;i<population[0].seq1.size();i++) {
		nuc n;
		n.A=0;
		n.C=0;
		n.G=0;
		n.T=0;
		nucs.push_back(n);
	}
	//Collect bases at each locus
	for (int i=0;i<population.size();i++) {
		for (int j=0;j<population[i].seq1.size();j++) {
			
			if (population[i].seq1[j]=='A') {
				nucs[j].A=nucs[j].A+population[i].size;
			}
			if (population[i].seq1[j]=='C') {
				nucs[j].C=nucs[j].C+population[i].size;
			}
			if (population[i].seq1[j]=='G') {
				nucs[j].G=nucs[j].G+population[i].size;
			}
			if (population[i].seq1[j]=='T') {
				nucs[j].T=nucs[j].T+population[i].size;
			}
			if (population[i].seq2[j]=='A') {
				nucs[j].A=nucs[j].A+population[i].size;
			}
			if (population[i].seq2[j]=='C') {
				nucs[j].C=nucs[j].C+population[i].size;
			}
			if (population[i].seq2[j]=='G') {
				nucs[j].G=nucs[j].G+population[i].size;
			}
			if (population[i].seq2[j]=='T') {
				nucs[j].T=nucs[j].T+population[i].size;
			}
		}
	}
}

void GetAlleles (run_params p,vector<amap>& alleles, vector<nuc> nucs) {
	ofstream all_file;

	
	all_file.open("Alleles.dat");
	
	for (unsigned int i=0;i<nucs.size();i++) {
		amap a;
		int max=0;
		if (nucs[i].A>max) {
			a.n1='A';
			max=nucs[i].A;
		}
		if (nucs[i].C>max) {
			a.n1='C';
			max=nucs[i].C;
		}
		if (nucs[i].G>max) {
			a.n1='G';
			max=nucs[i].G;
		}
		if (nucs[i].T>max) {
			a.n1='T';
			max=nucs[i].T;
		}
		max=0;
		int def=0;
		if (a.n1!='A') {
			if (nucs[i].A>max) {
				a.n2='A';
				def=1;
				max=nucs[i].A;
			}
		}
		if (a.n1!='C') {
			if (nucs[i].C>max) {
				a.n2='C';
				def=1;
				max=nucs[i].C;
			}
		}
		if (a.n1!='G') {
			if (nucs[i].G>max) {
				a.n2='G';
				def=1;
				max=nucs[i].G;
			}
		}
		if (a.n1!='T') {
			if (nucs[i].T>max) {
				a.n2='T';
				def=1;
				max=nucs[i].T;
			}
		}

		
		
		if (def==0) {
			a.n2='A';
		}
		all_file << "i= " << i << " " << nucs[i].A << " " << nucs[i].C << " " << nucs[i].G << " " << nucs[i].T << " " << a.n1 << " " << a.n2 << "\n";
	
		alleles.push_back(a);
	}
}

void GetSelectedAllele (int pop_size, int L, int& sel_loc, char& allele, vector<nuc> nucs, gsl_rng *rgen) {
	int loc_pos=0;
	while (sel_loc==0) {
		int loc_id=floor(gsl_rng_uniform(rgen)*L);
		vector<int> nc;
		nc.push_back(nucs[loc_id].A);
		nc.push_back(nucs[loc_id].C);
		nc.push_back(nucs[loc_id].G);
		nc.push_back(nucs[loc_id].T);
		
		int i1=0;
		int t1=0;
		for (int i=0;i<nc.size();i++) {
			if (nc[i]>t1) {
				t1=nc[i];
				i1=i;
			}
		}
		t1=0;
		for (int i=0;i<nc.size();i++) {
			if (nc[i]>t1&&i!=i1) {
				t1=nc[i];
				loc_pos=i;
			}
		}
		if (t1>=pop_size/5) { //Polymorphism exists at frequency 10% or greater: pop_size is number of diploids
			
			sel_loc=loc_id;
			
		}
	}
	if (loc_pos==0) {
		allele='A';
	} else if (loc_pos==1) {
		allele='C';
	} else if (loc_pos==2) {
		allele='G';
	} else if (loc_pos==3) {
		allele='T';
	}
}

void GetSelectedAllele2 (int pop_size, int L, int& sel_loc, char& allele, vector<nuc> nucs, gsl_rng *rgen) {
	int loc_pos=0;
	int loc_id=sel_loc;
	vector<int> nc;
	nc.push_back(nucs[loc_id].A);
	nc.push_back(nucs[loc_id].C);
	nc.push_back(nucs[loc_id].G);
	nc.push_back(nucs[loc_id].T);
	
	int i1=0;
	int t1=0;
	for (int i=0;i<nc.size();i++) {
		if (nc[i]>t1) {
			t1=nc[i];
			i1=i;
		}
	}
	t1=0;
	for (int i=0;i<nc.size();i++) {
		if (nc[i]>t1&&i!=i1) {
			t1=nc[i];
			loc_pos=i;
		}
	}
	if (loc_pos==0) {
		allele='A';
	} else if (loc_pos==1) {
		allele='C';
	} else if (loc_pos==2) {
		allele='G';
	} else if (loc_pos==3) {
		allele='T';
	}
}




double FindFit(int sel_loc, char allele, string seq1, string seq2, vector<double> fitness, double h) {
	double fit=0;
	
	if (seq1[sel_loc]==allele) { //Single locus
		if (seq2[sel_loc]==allele) {
			fit=1+fitness[sel_loc];
		} else {
			fit=1+(h*fitness[sel_loc]);
		}
	} else {
		if (seq2[sel_loc]==allele) {
			fit=1+(h*fitness[sel_loc]);
		} else {
			fit=1;
		}
	}
	return fit;
}


void FindHapMutationLoci (int L, double mu, vector<haplo>& population, vector<hap_loc>& mutants, gsl_rng *rgen) {
	//Generates vector mutants of new haplotypes
	for (unsigned int i=0;i<population.size();i++) {
		int nMut=gsl_ran_poisson(rgen,population[i].size*2*L*mu);// Number of mutations in this haplotype - factor of 2 for diploid
		
		hap_loc new_mut;
		for (int j=1;j<=nMut;j++) {
			unsigned int ran=floor(gsl_rng_uniform(rgen)*L);
			int same=0;
			for (unsigned int k=0;k<mutants.size();k++){
				if ((mutants[k].hap==i)&&(mutants[k].loc==ran)){
					mutants[k].value++;
					same=1;
					break;
				}
			}
			if (same==0) {
				new_mut.hap=i;
				new_mut.loc=ran;
				new_mut.value=1;
				mutants.push_back(new_mut);
			}
		}
	}
}

void MakeMutatedHaplotype (int sel_loc, char allele, double h, vector<double> fitness, vector<hap_loc>& mutants, vector<amap> alleles, vector<haplo>& population, gsl_rng *rgen) {
	for (unsigned int i=0;i<mutants.size();i++){
		haplo hap = population[mutants[i].hap];
		int ran=floor(gsl_rng_uniform(rgen)*2);
		if (ran==0) {
			
			if (hap.seq1[mutants[i].loc]==alleles[i].n1){
				hap.seq1[mutants[i].loc]=alleles[i].n2;
			} else {
				hap.seq1[mutants[i].loc]=alleles[i].n1;
			}
			
		} else {
			
			if (hap.seq2[mutants[i].loc]==alleles[i].n1){
				hap.seq2[mutants[i].loc]=alleles[i].n2;
			} else {
				hap.seq2[mutants[i].loc]=alleles[i].n1;
			}
			
		}
		hap.size=mutants[i].value;
		hap.fitness=FindFit(sel_loc,allele,hap.seq1,hap.seq2,fitness,h);
		population.push_back(hap);
		population[mutants[i].hap].size=population[mutants[i].hap].size-mutants[i].value;
	}
}

void CombineHapSame (int L, vector<haplo>& population){
	for (unsigned int i=0;i<population.size()-1;i++) {
		for (unsigned int j=i+1;j<population.size();j++) {
			int diff=1;
			if (population[i].seq1==population[j].seq1&&population[i].seq2==population[j].seq2) {
				diff=0;
			}
			if (population[i].seq1==population[j].seq2&&population[i].seq2==population[j].seq1) {
				diff=0;
			}
			if (diff==0) {
				population[i].size=population[i].size+population[j].size;
				population[j].size=0;
			}
		}
	}
}

void DeleteEmptyHap (vector<haplo>& population, vector<haplo>& pop_temp) {
	pop_temp.clear();
	
	for (unsigned int i=0;i<population.size();i++) {
		if (population[i].size>0) {
			pop_temp.push_back(population[i]);
			
		}
	}
}

void CopyHap (vector<haplo>& population, vector<haplo>& pop_temp) {
	pop_temp.clear();
	for (unsigned int i=0;i<population.size();i++) {
		pop_temp.push_back(population[i]);
	}
}

double GlobalFit (vector<haplo>& population){
	double glob_fit=0;
	for (unsigned int i=0;i<population.size();i++) {
		glob_fit=glob_fit+(population[i].size*population[i].fitness);
	}
	return(glob_fit);
}

void HapSelProbs(double glob_fit, vector<double>& p_select, vector<haplo>& population){
	for (unsigned int i=0;i<population.size();i++) {
		double ps = (population[i].size*population[i].fitness)/glob_fit;
		p_select.push_back(ps);
		
	}
}

void HapNextGeneration (int Nt, gsl_rng *rgen, vector<double>& p_select, vector<haplo>& population) {
	//Probabilities for multinomial are probabilities of selecting from each haplotype, given fitness
	double prob_select_mirror[population.size()];
	unsigned int next_gen[population.size()];
	for (unsigned int i=0;i<population.size();i++) {
		prob_select_mirror[i]=p_select[i];
	}
	gsl_ran_multinomial(rgen,population.size(),Nt,prob_select_mirror,next_gen);

	for (unsigned int i=0;i<population.size();i++) {

		population[i].size=next_gen[i];
	}
}




void ChooseParents (int N, vector<inpair>& parents, vector<haplo> population, gsl_rng *rgen) {

	int a[N];
	for (int i=0;i<N;i++) {
		a[i]=i;
	}
	gsl_ran_shuffle(rgen,a,N,sizeof(int)); //Randomly order previous generation
	
	for (int t=0;t<N;t++) {
		int k=0;
		
		for (unsigned int i=0;i<population.size();i++) {
		
			for (unsigned int j=0;j<population[i].size;j++) {
		
				if (a[k]==t) {
		
					parents[t].seq1=population[i].seq1;
					if (t+1<N) {
						parents[t+1].seq2=population[i].seq2;
					} else {
						parents[0].seq2=population[i].seq2;
					}
				}
				if (k+1<N) {
					if (a[k+1]==t) {
		
						parents[t].seq1=population[i].seq2;
						if (t+1<N) {
							parents[t+1].seq2=population[i].seq1;
						} else {
							parents[0].seq2=population[i].seq1;
						}
					}
				} else {
					if (a[0]==t) {
		
						parents[t].seq1=population[i].seq2;
						if (t+1<N) {
							parents[t+1].seq2=population[i].seq1;
						} else {
							parents[0].seq2=population[i].seq1;
						}
					}
				}
				k++;
			}
		}
	
	}

}

void CalculateRecombination (int L, double rho, int gen_length, int sel_loc, char allele, double h, vector<double> fitness, vector<int> diff_pos, vector<inpair>& parents, vector<haplo>& population, gsl_rng *rgen) {

	int loc;
	double rec=rho*gen_length;
    double rc=0;
	string s1;
	string s2;
	int rec_count=0;

	for (unsigned int i=0;i<parents.size();i++) {
		unsigned int nr=gsl_ran_poisson(rgen,rec);
		rec_count=rec_count+nr;
		
		if (nr>0) {
			
			for (unsigned int r=0;r<nr;r++) {
				rc=gsl_rng_uniform(rgen);
				loc=floor(rc*gen_length);
		
				//Search through locations
				int loc_index=0;
				for (int j=0;j<diff_pos.size();j++) {
					if (diff_pos[j]>loc) {
						loc_index=j;
						break;
					}
				}
		
				s1=parents[i].seq1;
				s2=parents[i].seq2;
				
				parents[i].seq1.replace(parents[i].seq1.begin()+loc_index,parents[i].seq1.end(),s2.begin()+loc_index,s2.end());
				parents[i].seq2.replace(parents[i].seq2.begin()+loc_index,parents[i].seq2.end(),s1.begin()+loc_index,s1.end());
			}
		}
	}

}

void ConstructPopulation (int sel_loc, char allele, double h, vector<double> fitness, vector<inpair> parents, vector<haplo>& population) {

	population.clear();
	
	for (int i=0;i<parents.size();i++) {
	
		haplo hap1;
		hap1.seq1=parents[i].seq1;
		hap1.seq2=parents[i].seq2;
		hap1.size=1;
		hap1.fitness=FindFit(sel_loc,allele,hap1.seq1,hap1.seq2,fitness,h);
		population.push_back(hap1);
	}
}

void PrintParents (vector<inpair> parents) {
	for (unsigned int i=0;i<parents.size();i++) {
		cout << "Parents " << i << "\n";
		cout << parents[i].seq1 << "\n";
		cout << parents[i].seq2 << "\n";
		cout << parents[i].clone1 << " " << parents[i].clone2 << "\n";
	}
}

void AlleleFreqs (vector<amap> alleles, vector<haplo> population, vector<double>& freqs) {
	int total=0;
	
	for (int i=0;i<population.size();i++) {
		total=total+2*population[i].size;
		for (int j=0;j<population[i].seq1.size();j++) {
			if (i==0) {
				freqs.push_back(0);
			}
			if (alleles[j].n1=='A') {
				if (population[i].seq1.compare(j,1,"A")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
				if (population[i].seq2.compare(j,1,"A")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
			}
			if (alleles[j].n1=='C') {
				if (population[i].seq1.compare(j,1,"C")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
				if (population[i].seq2.compare(j,1,"C")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
			}
			if (alleles[j].n1=='G') {
				if (population[i].seq1.compare(j,1,"G")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
				if (population[i].seq2.compare(j,1,"G")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
			}
			if (alleles[j].n1=='T') {
				if (population[i].seq1.compare(j,1,"T")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
				if (population[i].seq2.compare(j,1,"T")==0) {
					freqs[j]=freqs[j]+population[i].size;
				}
			}
		}
	}
}



