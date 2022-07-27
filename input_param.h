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
		double pop1_upper,pop1_lower,pop2_upper,pop2_lower;
		int min_sample,thread_num,analysis_mode;
		string ann_flag,ann_file,var_list_file,pop_file,vcf_file,output_file;
		set<string> pop1,pop2;
		bool no_splicing,dist_mode,exhaust_disc_mode,exhaust_valid_mode,output_per_sample_gt;

		void reset_freq_param(){
			pop1_upper=1;
			pop1_lower=0;
			pop2_upper=1;
			pop2_lower=0;
			min_sample=3;
			ann_flag="";
			no_splicing=false;
			dist_mode=false;
			exhaust_disc_mode=false;
			exhaust_valid_mode=false;
			output_per_sample_gt=false;
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
			for (unordered_map<string,set<string> >::iterator i = b.begin(); i!=b.end(); ++i){
				pop1.insert(i->first);
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
		map<string,string> file_dict,param_dict;
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
		}else if(cur_str=="--flag"){
			idx++;
			ann_flag=argv[idx];
			param_dict["Target annotation flag: "]=argv[idx];
			idx++;
		}else if(cur_str=="--disc"){
			exhaust_disc_mode=true;
			param_dict["Analysis mode 1: "]="Exhaustive discovery mode";
			idx++;
		}else if(cur_str=="--valid"){
			exhaust_valid_mode=true;
			param_dict["Analysis mode 2: "]="Exhaustive validation mode";
			idx++;
		}else if(cur_str=="--d-v"){
			exhaust_disc_mode=true;
			exhaust_valid_mode=true;
			param_dict["Analysis mode 3: "]="Exhaustive discovery-validation mode";
			idx++;
		}else if(cur_str=="--dist"){
			dist_mode=true;
			param_dict["Analysis mode 4: "]="Output variant frequency for all labeling group";
			idx++;
		}else if(cur_str=="--mins"){
			idx++;
			min_sample=stoi(argv[idx]);
			param_dict["Minimum number of samples allowed for population inclusion "]=to_string(min_sample);			
			idx++;
		}else if(cur_str=="--p1l"){
			idx++;
			pop1_lower=stod(argv[idx]);
			param_dict["Pop1 lower frequency boundary: "]=to_string(pop1_lower);			
			idx++;
		}else if(cur_str=="--p1u"){
			idx++;
			pop1_upper=stod(argv[idx]);
			param_dict["Pop1 upper frequency boundary: "]=to_string(pop1_upper);			
			idx++;
		}else if(cur_str=="--no-splicing"){
			idx++;
			no_splicing=true;
			param_dict["No spicing mode"]="";			
		}else{
			cout<<"Invalid argument: "<<cur_str<<endl;
			cout<<"Tool terminated."<<endl;
			exit(-1);
		}
	}

	cout<<"File parameters: "<<endl;
	for (map<string,string>::iterator i = file_dict.begin(); i !=file_dict.end(); ++i){
		cout<<"--"<<i->first<<i->second<<endl;
	}

	cout<<"Analysis parameters: "<<endl;
	for (map<string,string>::iterator i = param_dict.begin(); i !=param_dict.end(); ++i){
		cout<<"--"<<i->first<<i->second<<endl;
	}	
	cout<<endl;

}

void input_param::launch_analysis_module(){
	if(param_dict.count("Analysis mode 1: ")){
		analysis_mode=1;
	}else if(param_dict.count("Analysis mode 2: ")){
		analysis_mode=2;
	}else if(param_dict.count("Analysis mode 3: ")){
		analysis_mode=3;
	}else if(param_dict.count("Analysis mode 4: ")){
		analysis_mode=4;
	}
}


#endif








