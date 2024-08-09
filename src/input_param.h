#ifndef PARAM_H   
#define PARAM_H
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<unordered_map>
#include<set>
#include<map>
#include "utilities.h"
using namespace std;


class input_param{
	public:
		double pop1_upper,pop1_lower,pop2_upper,pop2_lower,af,sample_frac,mean_llh;
		int eff_sample,min_sample,min_sample_tar,min_sample_ref,thread_num,analysis_mode,pop_num,min_depth,max_admix_pop,max_homo_pop,seed,min_rep_size,rep_num,experiment_times;
		string ann_flag,ann_file,var_list_file,pop_file,vcf_file,output_file,depth_file,likelihood_file;
		set<string> pop1,pop2,pop_all;
		bool verbose,ml_mode,lh_mode,StrPrint,isGS,isEXP,no_splicing,kmeans,flip,dist_mode,exhaust_disc_mode,exhaust_valid_mode,output_per_sample_gt,bipop_mode,unipop_mode,gene_mode,unique_valid_mode,combine_all,STR_mode,INDEL_mode;

		void reset_freq_param(){
			min_depth=0;
			pop1_upper=1;
			pop1_lower=0.5;
			pop2_upper=0.5;
			pop2_lower=0;
			min_sample=3;
			eff_sample=-1;
			min_sample_tar=-1;
			min_sample_ref=-1;
			min_rep_size=1;
			max_admix_pop=0;
			max_homo_pop=1;
			seed=123;
			sample_frac=1;
			experiment_times=10;
			af=0.0;
			mean_llh=1.0;
			pop_num=0;
			rep_num=1;
			ann_flag="";
			var_list_file="";
			depth_file="";
			flip=false;
			verbose=false;
			no_splicing=false;
			dist_mode=false;
			gene_mode=false;
			ml_mode=false;
			exhaust_disc_mode=false;
			exhaust_valid_mode=false;
			output_per_sample_gt=false;
			unique_valid_mode=false;
			bipop_mode=false;
			STR_mode=false;
			INDEL_mode=false;
			combine_all=false;
			unipop_mode=false;
			lh_mode=false;
			kmeans=false;
			isEXP=false;
			isGS=false;
			StrPrint=false;
		}

		void print_input_parameters(){
			if(!no_splicing){
				cout<<"Genes with splicing variants allowed."<<endl;
			}else{
				cout<<"Only genes without splicing variants are allowed."<<endl;
			}

			if(ann_flag.length()>0){
				cout<<"Qualifying variants with annotation flag: "<<ann_flag<<endl;
			}else{
				cout<<"No flag in annotation selected."<<endl;
			}

			if(!(dist_mode||exhaust_disc_mode)){
				if(exhaust_valid_mode){
					cout<<"Variant frequency list for "<<pop1.size()<<" selected populations will be outputed"<<endl;
					cout<<"Target frequency range for the validation population: ["<<pop1_lower<<","<<pop1_upper<<"]"<<endl;	
					cout<<"Minimum number of samples for inclusion: "<<min_sample<<endl;				
				}else{
					cout<<"Variant frequency list for "<<pop1.size()+pop2.size()<<" selected populations will be outputed"<<endl;
					cout<<"Target frequency range for population 1: ["<<pop1_lower<<","<<pop1_upper<<"]"<<endl;
					dp.print(pop1);
					cout<<"Target frequency range for population 2: ["<<pop2_lower<<","<<pop2_upper<<"]"<<endl;
					dp.print(pop2);
					cout<<"Minimum number of samples for inclusion: "<<min_sample<<endl;					
				}
			}else{
				if(exhaust_disc_mode){
					cout<<"Target frequency range for any possible population: ["<<pop1_lower<<","<<pop1_upper<<"]"<<endl;
					cout<<"Minimum number of samples for inclusion: "<<min_sample<<endl;
				}

			}
			cout<<endl;
		}

		void load_population_label(string a,vector<string> b){
			pop1.clear();
			pop2.clear();
			pop1.insert(a);
			for (int i = 0; i < b.size(); ++i){
				pop2.insert(b[i]);
			}
		}

		void load_population_label(string a,string b){
			pop1.clear();
			pop2.clear();
			pop1.insert(a);
			pop2.insert(b);
		}

		void load_population_label(vector<string> a,string b){
			pop1.clear();
			pop2.clear();
			pop2.insert(b);
			for (int i = 0; i < a.size(); ++i){
				pop1.insert(a[i]);
			}
		}

		void load_population_label(vector<string> a,vector<string> b){
			pop1.clear();
			pop2.clear();
			for (int i = 0; i < a.size(); ++i){
				pop1.insert(a[i]);
			}
			for (int i = 0; i < b.size(); ++i){
				pop2.insert(b[i]);
			}
		}

		void load_population_label(vector<string> a,set<string> b){
			pop1.clear();
			pop2.clear();
			for (int i = 0; i < a.size(); ++i){
				pop1.insert(a[i]);
			}
			for(set<string>::iterator i=b.begin();i!=b.end();++i){
				if(!pop1.count(*i)){
					pop2.insert(*i);
				}
			}
		}

		void load_population_label(string a,unordered_map<string,set<string> > b){
			pop1.clear();
			pop2.clear();
			pop1.insert(a);
			for (unordered_map<string,set<string> >::iterator i = b.begin(); i!=b.end(); ++i){
				if(i->first!=a){
					pop2.insert(i->first);
				}
			}
		}

		void load_population_label(unordered_map<string,set<string> > b){
			pop1.clear();
			pop2.clear();
			pop_all.clear();
			for (unordered_map<string,set<string> >::iterator i = b.begin(); i!=b.end(); ++i){
				pop1.insert(i->first);
				pop_all.insert(i->first);
			}
		}

		void print_population_info(){
			if(!(dist_mode||exhaust_disc_mode)){
				cout<<"Label(s) to be included in population 1: "<<endl;
				dp.print(pop1);
				cout<<"Label(s) to be included in population 2: "<<endl;
				dp.print(pop2);
				cout<<endl;
			}
		}

		void read_parameters(int argc, char const *argv[]);

		void launch_analysis_module();

	private:
		data_printer dp;
		map<string,string> file_dict,param_dict,analysis_dict;
};

void input_param::read_parameters(int argc, char const *argv[]){
	int idx=1;
	string cur_str;
	while(idx<argc){
		cur_str=argv[idx];
		if(cur_str=="-vl"){
			idx++;
			var_list_file=argv[idx];
			file_dict["Varaint list file: "]=argv[idx];
			idx++;
		}else if(cur_str=="-p"){
			idx++;
			pop_file=argv[idx];
			file_dict["Sample label file: "]=argv[idx];
			idx++;
		}else if(cur_str=="-i"){
			idx++;
			vcf_file=argv[idx];
			file_dict["Input variant file: "]=argv[idx];
			idx++;
		}else if(cur_str=="-o"){
			idx++;
			output_file=argv[idx];
			file_dict["Output variant list file: "]=argv[idx];
			idx++;
		}else if(cur_str=="-o-likelihood"){
			idx++;
			likelihood_file=argv[idx];
			file_dict["Output likelihood table to: "]=argv[idx];
			idx++;
		}else if(cur_str=="-t"){
			idx++;
			thread_num=stoi(argv[idx]);
			param_dict["Number of threads: "]=argv[idx];
			idx++;
		}else if(cur_str=="-a"){
			idx++;
			ann_file=argv[idx];
			file_dict["Annotation database file: "]=argv[idx];
			idx++;
		}else if(cur_str=="-d"){
			idx++;
			depth_file=argv[idx];
			file_dict["Sample depth file: "]=argv[idx];
			idx++;
		}else if(cur_str=="disc"){
			exhaust_disc_mode=true;
			analysis_dict["Analysis mode 1: "]="Exhaustive discovery mode";
			idx++;
		}else if(cur_str=="valid"){
			exhaust_valid_mode=true;
			analysis_dict["Analysis mode 2: "]="Exhaustive validation mode";
			idx++;
		}else if(cur_str=="d-v"){
			exhaust_disc_mode=true;
			exhaust_valid_mode=true;
			analysis_dict["Analysis mode 3: "]="Exhaustive discovery-validation mode";
			idx++;
		}else if(cur_str=="Dist"){
			dist_mode=true;
			analysis_dict["Analysis mode 4: "]="Output variant frequency for all labeling group";
			idx++;
		}else if(cur_str=="SigFreq"){
			unipop_mode=true;
			analysis_dict["Analysis mode 5: "]="Unique pattern search for all group";
			idx++;
		}else if(cur_str=="freq"){
			unipop_mode=true;
			analysis_dict["Analysis mode 6: "]="Flexible frequency analysis mode";
			idx++;
		}else if(cur_str=="gene"){
			gene_mode=true;
			analysis_dict["Analysis mode 7: "]="Variant aggregation analysis by gene";
			idx++;
		}else if(cur_str=="uvalid"){
			unique_valid_mode=true;
			analysis_dict["Analysis mode 8: "]="Validate unique variants";
			idx++;
		}else if(cur_str=="STRdisc"){
			STR_mode=true;
			analysis_dict["Analysis mode 9: "]="Discover Short Tandem Repeats(STR)";
			idx++;
		}else if(cur_str=="SigLh"){
			lh_mode=true;
			analysis_dict["Analysis mode 10: "]="Likelihood based profiling";
			idx++;
		}else if(cur_str=="SigMl"){
			ml_mode=true;
			analysis_dict["Analysis mode 11: "]="Selflearning likelihood based profiling";
			idx++;
		}else if(cur_str=="--tar-lower"){
			idx++;
			pop1_lower=stod(argv[idx]);
			param_dict["Target group lower signature variant frequency boundary: "]=to_string(pop1_lower);			
			idx++;
		}else if(cur_str=="--tar-upper"){
			idx++;
			pop1_upper=stod(argv[idx]);
			param_dict["Target group upper signature variant frequency boundary: "]=to_string(pop1_upper);			
			idx++;
		}else if(cur_str=="--ref-lower"){
			idx++;
			pop2_lower=stod(argv[idx]);
			param_dict["Reference group lower non-signature variant frequency boundary: "]=to_string(pop2_lower);			
			idx++;
		}else if(cur_str=="--ref-upper"){
			idx++;
			pop2_upper=stod(argv[idx]);
			param_dict["Reference group upper non-signature variant frequency boundary: "]=to_string(pop2_upper);			
			idx++;
		}else if(cur_str=="--group-num"){
			idx++;
			pop_num=stoi(argv[idx]);
			param_dict["Number of groups included for shared pattern discovery: "]=to_string(pop_num);			
			idx++;
		}else if(cur_str=="--min-total-sample"){
			idx++;
			eff_sample=stoi(argv[idx]);
			param_dict["Minimum number of total effective samples for locus inclusion: "]=to_string(eff_sample);			
			idx++;
		}else if(cur_str=="--min-sample"){
			idx++;
			min_sample=stoi(argv[idx]);
			param_dict["Minimum number of samples allowed for population inclusion: "]=to_string(min_sample);			
			idx++;
		}else if(cur_str=="--min-tar"){
			idx++;
			min_sample_tar=stoi(argv[idx]);
			param_dict["Minimum number of samples allowed for target population inclusion: "]=to_string(min_sample_tar);			
			idx++;
		}else if(cur_str=="--min-ref"){
			idx++;
			min_sample_ref=stoi(argv[idx]);
			param_dict["Minimum number of samples allowed for reference population inclusion: "]=to_string(min_sample_ref);			
			idx++;
		}else if(cur_str=="--min-depth"){
			idx++;
			min_depth=stoi(argv[idx]);
			param_dict["Minimum read depth for a variant to be included in the analysis: "]=to_string(min_depth);			
			idx++;
		}else if(cur_str=="--flag"){
			idx++;
			ann_flag=argv[idx];
			param_dict["Target annotation flag: "]=argv[idx];
			idx++;
		}else if(cur_str=="--subsample-seed"){
			idx++;
			seed=stoi(argv[idx]);
			param_dict["Random seed used for dataset partition: "]=argv[idx];
			idx++;
		}else if(cur_str=="--subsample-frac"){
			idx++;
			sample_frac=stod(argv[idx]);
			param_dict["Sampling fraction: "]=argv[idx];
			idx++;
		}else if(cur_str=="--subsample-run"){
			idx++;
			experiment_times=stoi(argv[idx]);
			param_dict["Repetitions of subsampling runs: "]=argv[idx];
			idx++;
		}else if(cur_str=="--lh-thres"){
			idx++;
			mean_llh=stod(argv[idx]);
			param_dict["Minimum average likelihood for signature qualification: "]=argv[idx];
			idx++;
		}else if(cur_str=="--min-rep-size"){
			idx++;
			min_rep_size=stoi(argv[idx]);
			param_dict["Minimum size of the repetitive units in STR: "]=to_string(min_rep_size);			
			idx++;
		}else if(cur_str=="--str-kmeans"){
			kmeans=true;
			param_dict["Choosing STR classification boundary using k-means"]="";
			idx++;
		}else if(cur_str=="--str-long"){
			isEXP=true;
			param_dict["Targetting STR with long expansions"]="";
			idx++;
		}else if(cur_str=="--GS"){
			isGS=true;
			param_dict["Targetting at GS"]="";
			idx++;
		}else if(cur_str=="--str-printall"){
			StrPrint=true;
			param_dict["Printing the header for all STR loci"]="";
			idx++;
		}else if(cur_str=="--repnum"){
			idx++;
			rep_num=stoi(argv[idx]);
			param_dict["Cutoff point for number of the repetitive units in STR: "]=to_string(rep_num);			
			idx++;
		}else if(cur_str=="--max-admix"){
			idx++;
			max_admix_pop=stoi(argv[idx]);
			param_dict["Maximum number of admixed population allowed in the analysis: "]=to_string(max_admix_pop);			
			idx++;
		}else if(cur_str=="--max-homo-pop"){
			idx++;
			max_homo_pop=stoi(argv[idx]);
			param_dict["Maximum number of homogeneous population allowed in the analysis: "]=to_string(max_homo_pop);			
			idx++;
		}else if(cur_str=="--af"){
			idx++;
			af=stod(argv[idx]);
			param_dict["Minumum minor allele frequency allowed for a site to be included: "]=to_string(af);			
			idx++;
		}else if(cur_str=="--combine-all"){
			combine_all=true;
			param_dict["Combine all populations during unique variant discovery "]="";			
			idx++;
		}else if(cur_str=="--no-splicing"){
			idx++;
			no_splicing=true;
			param_dict["No spicing mode"]="";			
		}else if(cur_str=="--flip"){
			idx++;
			flip=true;
			param_dict["Flip the ref/alt allele"]="";			
		}else if(cur_str=="--verbose"){
			idx++;
			verbose=true;
			param_dict["Run in verbose mode"]="";			
		}else if(cur_str=="--INDEL"){
			idx++;
			INDEL_mode=true;
			param_dict["INDEL mode activated"]="";			
		}else{
			cout<<"Invalid argument: "<<cur_str<<endl;
			cout<<"Tool terminated."<<endl;
			exit(-1);
		}
	}

	cout<<"Analysis mode: "<<endl;
	for (map<string,string>::iterator i = analysis_dict.begin(); i !=analysis_dict.end(); ++i){
		cout<<i->second<<endl;
	}

	cout<<"File parameters: "<<endl;
	for (map<string,string>::iterator i = file_dict.begin(); i !=file_dict.end(); ++i){
		cout<<"--"<<i->first<<i->second<<endl;
	}

	cout<<"Analysis parameters: "<<endl;
	if(!analysis_dict.count("Analysis mode 4: ")){
		for (map<string,string>::iterator i = param_dict.begin(); i !=param_dict.end(); ++i){
			cout<<"--"<<i->first<<i->second<<endl;
		}	
	}else{
		if(af!=0){
			cout<<"--MAF: "<<af<<endl;
		}
	}
	cout<<endl;

}

void input_param::launch_analysis_module(){
	if(analysis_dict.count("Analysis mode 1: ")){
		analysis_mode=1;
	}else if(analysis_dict.count("Analysis mode 2: ")){
		analysis_mode=2;
	}else if(analysis_dict.count("Analysis mode 3: ")){
		analysis_mode=3;
	}else if(analysis_dict.count("Analysis mode 4: ")){
		analysis_mode=4;
	}else if(analysis_dict.count("Analysis mode 5: ")){
		analysis_mode=5;
	}else if(analysis_dict.count("Analysis mode 6: ")){
		analysis_mode=6;
	}else if(analysis_dict.count("Analysis mode 7: ")){
		analysis_mode=7;
	}else if(analysis_dict.count("Analysis mode 8: ")){
		analysis_mode=8;
	}else if(analysis_dict.count("Analysis mode 9: ")){
		analysis_mode=9;
	}else if(analysis_dict.count("Analysis mode 10: ")){
		analysis_mode=10;
	}else if(analysis_dict.count("Analysis mode 11: ")){
		analysis_mode=11;
	}
}


#endif








