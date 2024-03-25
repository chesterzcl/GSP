#ifndef MAM_H   
#define MAM_H
#include "utilities.h"
#include "pop_data.h"
#include "var_list.h"
#include "ann_data.h"
#include "input_param.h"

using namespace std;

class main_analysis_module{
	public:
		vector<string> input_vec,output_vec,output_gt_vec;

		void output_file_header(ofstream& output,var_list& var,input_param param,pop_data pop,vector<string>& pop_vec){
			string temp_str;
			output<<"#Chromosome"<<'\t'<<"Position"<<'\t'<<"Gene"<<'\t'<<"Variant Type"<<'\t'<<"Mutation Position"<<'\t'<<"Amino Acid Change";
			if(param.dist_mode||param.exhaust_disc_mode||param.unipop_mode||param.bipop_mode||param.gene_mode||param.STR_mode||param.lh_mode){
				for (unordered_map<string,set<string> >::iterator i =  pop.pop_dict.begin(); i != pop.pop_dict.end(); ++i){
					if(i->second.size()>=param.min_sample){
						output<<'\t'<<i->first;
						pop_vec.push_back(i->first);
					}
				}
				output<<endl;
			}else{
				for (set<string>::iterator i = param.pop1.begin(); i != param.pop1.end(); ++i){
					output<<'\t'<<*i;
					pop_vec.push_back(*i);
				}
				for (set<string>::iterator i = param.pop2.begin(); i != param.pop2.end(); ++i){
					output<<'\t'<<*i;
					pop_vec.push_back(*i);
				}
				output<<endl;
			}
			//Output genotype header
			if(!param.dist_mode){
				if(param.output_per_sample_gt){
					temp_str+="#CHROMOSOME\tPOSITION\tGENE\tMUTATION TYPE\tMUTATION POSITION\tAMINO ACID CHANGE";
					for (int i = 0; i < pop_vec.size(); ++i){
						for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j){
							temp_str+='\t'+pop.col_sample_dict[*j];
						}
					}
					var.var_gt_mat.push_back(temp_str);
				}		
			}
		}

		void output_var_freq_data(ofstream& output,var_list var,vector<string> line_vec,vector<string> pop_vec,unordered_map<string,string>  freq_dict,int row_cur){
			output<<var.var_mat[row_cur][0]<<'\t';
			output<<var.var_mat[row_cur][1]<<'\t';
			output<<read_char_delim_str(line_vec[7],'|')[3]<<'\t';
			output<<read_char_delim_str(line_vec[7],'|')[1]<<'\t';
			output<<read_char_delim_str(line_vec[7],'|')[13]<<'\t';
			output<<read_char_delim_str(line_vec[7],'|')[10];
			for (int i = 0; i < pop_vec.size(); ++i){
				output<<'\t'<<freq_dict[pop_vec[i]];
			}
			output<<endl;			
		}

		void store_sample_gt_data(var_list& var,int row_cur,vector<string> line_vec,string& temp_str){
			temp_str=var.var_mat[row_cur][0]+'\t'+var.var_mat[row_cur][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10]+temp_str;
			var.var_gt_mat.push_back(temp_str);			
		}

		void store_var_pop_mat(var_list& var,vector<string> var_pop_vec,vector<string> temp_vec,int row_cur){
			temp_vec.push_back(var.var_mat[row_cur][0]);
			temp_vec.push_back(var.var_mat[row_cur][1]);
			for (int i = 0; i < var_pop_vec.size(); ++i){
				temp_vec.push_back(var_pop_vec[i]);
			}
			var.var_pop_mat.push_back(temp_vec);
		}

		void pop_frequency_analysis(string input_address,string output_address,pop_data pop,input_param param,var_list& var,ann_data ann);

		void calculate_var_pop_variant_frequency(pop_data pop,vector<string> pop_vec,input_param param,var_list& var,int row_cur);

		void variants_discovery_validation(string disc_var,pop_data pop1,string disc_op,string valid_var,pop_data pop2,string valid_op,ann_data ann,input_param param1,input_param param2,var_list& var);

		void single_pop_freq_analysis(string var_input_disc,string var_input_val,string& var_output,pop_data pop1,pop_data pop2,ann_data ann,input_param param1,input_param param2,var_list& var);

		void bi_pop_freq_analysis(string var_input_disc,string var_input_val,string var_output,pop_data pop1,pop_data pop2,ann_data ann,input_param param1,input_param param2,var_list& var);

		

	private:
		ifstream input;
		ofstream output;
		string line;
		vector<string> line_vec;
		set<int> empty_set;
};

void main_analysis_module::calculate_var_pop_variant_frequency(pop_data pop,vector<string> pop_vec,input_param param,var_list &var,int row_cur){
	unordered_map<string,string>  freq_dict;
	int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0;
	int total_num,case_num,total_num2,case_num2;
	double freq;
	vector<string> var_pop_vec,temp_vec;
	string temp_str="";
	for (int i = 0; i < pop_vec.size(); ++i){
		total_num=0;
		case_num=0;
		for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j){
			if(line_vec[*j][0]=='0'||line_vec[*j][0]=='1'){
				total_num+=1;
				case_num+=check_genotype(line_vec[*j]);
				temp_str+='\t'+to_string(check_genotype(line_vec[*j]));
			}else{
				if(line_vec[*j][0]=='.'){
					temp_str+="\t2";
				}else{
					temp_str+="\t3";
				}
			}
		}

		if(total_num>=param.min_sample){
			freq=(double)case_num/(double)total_num;
			freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
			if (param.exhaust_disc_mode){
				if (freq>=param.pop1_lower&&freq<=param.pop1_upper){
					var_pop_vec.push_back(pop_vec[i]);
					ind4*=0;
				}
			}else{
				if (!param.dist_mode){
					if(param.pop1.count(pop_vec[i])){
						valid_ind1+=1;
						if (!(freq>=param.pop1_lower&&freq<=param.pop1_upper)){
							ind1*=0;
						}else{
							for (int k = 0; k < var.var_mat[row_cur].size(); ++k){
								if(var.var_mat[row_cur][k]==pop_vec[i]){
									ind5*=0;
								}								
							}
						}
					}else{
						if(param.pop2.count(pop_vec[i])){
							valid_ind2+=1;
							if(!(freq>=param.pop2_lower&&freq<=param.pop2_upper)){
								ind2*=0;
							}
						}
					}					
				}
			}
		}else{
			if (total_num==0){
				freq_dict[pop_vec[i]]="0/0";
			}else{
				freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
			}
		}
		if (total_num==case_num){
			ind3*=1;
		}else{
			ind3*=0;
		}
	}
	if(param.dist_mode){
		output<<var.var_mat[row_cur][0]<<'\t'<<var.var_mat[row_cur][1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[3]<<'\t'<<read_char_delim_str(line_vec[7],'|')[1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[13]<<'\t'<<read_char_delim_str(line_vec[7],'|')[10];
		for (int i = 0; i < pop_vec.size(); ++i){
			output<<'\t'<<freq_dict[pop_vec[i]];
		}
		output<<endl;
	}else{
		if(param.exhaust_disc_mode){
			if (ind4==0&&ind3==0){
				temp_vec.push_back(var.var_mat[row_cur][0]);
				temp_vec.push_back(var.var_mat[row_cur][1]);
				for (int i = 0; i < var_pop_vec.size(); ++i){
					temp_vec.push_back(var_pop_vec[i]);
				}
				var.var_pop_mat.push_back(temp_vec);
				//Store per-individual gt info
				if(param.output_per_sample_gt){
					temp_str=var.var_mat[row_cur][0]+'\t'+var.var_mat[row_cur][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10]+temp_str;
					var.var_gt_mat.push_back(temp_str);
				}

				output<<var.var_mat[row_cur][0]<<'\t'<<var.var_mat[row_cur][1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[3]<<'\t'<<read_char_delim_str(line_vec[7],'|')[1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[13]<<'\t'<<read_char_delim_str(line_vec[7],'|')[10];
				for (int i = 0; i < pop_vec.size(); ++i){
					output<<'\t'<<freq_dict[pop_vec[i]];
				}
				output<<endl;				
			}
		}else{
			if(param.exhaust_valid_mode){
				if(valid_ind1!=0&&ind3==0){
					if(ind5==0){
						//Store per-individual gt info
						if(param.output_per_sample_gt){
							temp_str=var.var_mat[row_cur][0]+'\t'+var.var_mat[row_cur][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10]+temp_str;
							var.var_gt_mat.push_back(temp_str);
						}

						output<<var.var_mat[row_cur][0]<<'\t'<<var.var_mat[row_cur][1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[3]<<'\t'<<read_char_delim_str(line_vec[7],'|')[1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[13]<<'\t'<<read_char_delim_str(line_vec[7],'|')[10];
						for (int i = 0; i < pop_vec.size(); ++i){
							output<<'\t'<<freq_dict[pop_vec[i]];
						}
						output<<endl;
					}
				}
			}else{
				if(valid_ind1*valid_ind2!=0&&valid_ind1==param.pop1.size()&&ind3==0){
					if(ind1*ind2!=0){
						temp_vec.push_back(var.var_mat[row_cur][0]);
						temp_vec.push_back(var.var_mat[row_cur][1]);
						for (int i = 0; i < var_pop_vec.size(); ++i){
							temp_vec.push_back(var_pop_vec[i]);
						}
						var.var_pop_mat.push_back(temp_vec);

						if(param.output_per_sample_gt){
							temp_str=var.var_mat[row_cur][0]+'\t'+var.var_mat[row_cur][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10]+temp_str;
							var.var_gt_mat.push_back(temp_str);
						}

						output<<var.var_mat[row_cur][0]<<'\t'<<var.var_mat[row_cur][1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[3]<<'\t'<<read_char_delim_str(line_vec[7],'|')[1]<<'\t'<<read_char_delim_str(line_vec[7],'|')[13]<<'\t'<<read_char_delim_str(line_vec[7],'|')[10];
						for (int i = 0; i < pop_vec.size(); ++i){
							output<<'\t'<<freq_dict[pop_vec[i]];
						}
						output<<endl;
				  	}
				}	
			}		
		}
	}
}

void main_analysis_module::pop_frequency_analysis(string input_address,string output_address,pop_data pop,input_param param,var_list& var,ann_data ann){
	input.open(input_address);
	check_file_open_status(input,input_address);
	output.open(output_address);
	check_file_open_status(output,output_address);
	int row_cur=0,row_cter=1;
	vector<string> pop_vec;
	string temp_str;
	output<<"#CHROMOSOME"<<'\t'<<"POSITION"<<'\t'<<"GENE"<<'\t'<<"MUTATION TYPE"<<'\t'<<"MUTATION POSITION"<<'\t'<<"AMINO ACID CHANGE";
	if(param.dist_mode||param.exhaust_disc_mode){
		for (unordered_map<string,set<string> >::iterator i =  pop.pop_dict.begin(); i != pop.pop_dict.end(); ++i){
			output<<'\t'<<i->first;
			pop_vec.push_back(i->first);
		}
		output<<endl;
	}else{
		for (set<string>::iterator i = param.pop1.begin(); i != param.pop1.end(); ++i){
			output<<'\t'<<*i;
			pop_vec.push_back(*i);
		}
		for (set<string>::iterator i = param.pop2.begin(); i != param.pop2.end(); ++i){
			output<<'\t'<<*i;
			pop_vec.push_back(*i);
		}
		output<<endl;
	}
	 if(!param.dist_mode){
		if(param.output_per_sample_gt){
			temp_str+="#CHROMOSOME\tPOSITION\tGENE\tMUTATION TYPE\tMUTATION POSITION\tAMINO ACID CHANGE";
			for (int i = 0; i < pop_vec.size(); ++i){
				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j){
					temp_str+='\t'+pop.col_sample_dict[*j];
				}
			}
			var.var_gt_mat.push_back(temp_str);
		}		
	}
	while(getline(input,line)){
		if(line[0]!='#'){
			if(row_cter%10000==0){
				cout<<row_cter<<" variants scanned."<<endl;
			}
			row_cter++;
			if(row_cur==var.var_mat.size()){
				break;
			}
			line_vec=read_char_delim_str(line,'\t');
			if(var.var_mat[row_cur][0]==line_vec[0]){
				if(var.var_mat[row_cur][1]==line_vec[1]){
					if(!param.no_splicing){
						calculate_var_pop_variant_frequency(pop,pop_vec,param,var,row_cur);
					}else{
						if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1){
							calculate_var_pop_variant_frequency(pop,pop_vec,param,var,row_cur);
						}
					}
					row_cur++;
				}else{
					if(stoi(line_vec[1])>stoi(var.var_mat[row_cur][1])){
						row_cur++;
					}					
				}
			}else{
				if(stoi(line_vec[0].substr(7,2))>stoi(var.var_mat[row_cur][0].substr(7,2))){
					row_cur++;
				}
			}
		}
	}
	input.close();
	output.close();
}

void main_analysis_module::variants_discovery_validation(string disc_var,pop_data pop1,string disc_op,string valid_var,pop_data pop2,string valid_op,ann_data ann,input_param param1,input_param param2,var_list& var){
	data_printer dp;
	input_param param_disc=param1,param_valid=param2;

	cout<<"Loading data for the discovery process..."<<endl<<endl;
	var.load_variant_data(disc_var);

	param1.print_input_parameters();

	cout<<"Discovery parameters loaded. Starting exhaustive variant discovery..."<<endl<<endl;
	pop_frequency_analysis(disc_var,disc_op,pop1,param1,var,ann);
	cout<<"Discovery process completed."<<endl<<endl;

	if(param1.output_per_sample_gt){
		cout<<"Output per-sample gt data on variants of interest..."<<endl;
		var.output_gt_data(disc_op+".gt");		
	}

	cout<<"Loading data for the validation process..."<<endl<<endl;

	var.update_var_list();

	if(param2.exhaust_valid_mode){
		param_valid.load_population_label(pop1.pop_dict);
	}

	for (set<string>::iterator i = param2.pop1.begin(); i !=param2.pop1.end(); ++i){
		if(!pop2.pop_col_dict.count(*i)){
			param_valid.pop1.erase(*i);
		}
	}
	for (set<string>::iterator i = param2.pop2.begin(); i !=param2.pop2.end(); ++i){
		if(!pop2.pop_col_dict.count(*i)){
			param_valid.pop2.erase(*i);
		}
	}

	param_valid.print_input_parameters();

	cout<<"Validation parmeters loaded. Start the validtion for target variants..."<<endl<<endl;
	pop_frequency_analysis(valid_var,valid_op,pop2,param_valid,var,ann);

	cout<<"Validation process completed..."<<endl<<endl;

	if(param_valid.output_per_sample_gt){
		cout<<"Output per-sample gt data on variants of interest..."<<endl;
		output_gt_vec.push_back(valid_op+".gt");
		var.output_gt_data(valid_op+".gt");		
	}
}

void main_analysis_module::single_pop_freq_analysis(string var_input_disc,string var_input_val,string& var_output,pop_data pop1,pop_data pop2,ann_data ann,input_param param1,input_param param2,var_list& var){
	input_param param_disc=param1,param_valid=param2;
	data_printer dp;
	string disc_op,valid_op,splicing_status;
	if(param1.no_splicing){
		splicing_status="_no_splicing";
	}else{
		splicing_status="";
	}
	for (unordered_map<string,set<string> >::iterator i = pop1.pop_dict.begin();i != pop1.pop_dict.end(); ++i){
		cout<<"============================================================================================="<<endl;
		if(pop2.pop_col_dict.count(i->first)){
			if(pop2.pop_col_dict[i->first].size()<param2.min_sample){
				cout<<"Only "<<pop2.pop_col_dict[i->first].size()<<" sample from "<<i->first<<" were detected in the validation dataset, failed to meet minimum sample requirement. Analysis will be skipped."<<endl;
			}else{
				cout<<"\nStarting frequency analysis for "<<i->first<<endl;
				param_disc.load_population_label(i->first,pop1.pop_dict);
				param_valid.load_population_label(i->first,pop1.pop_dict);
				disc_op=var_input_disc+"."+i->first+"_"+to_string(param1.pop1_lower)+"_"+to_string(param1.pop2_upper)+"_"+to_string(param1.min_sample)+"_disc_op"+splicing_status;
				valid_op=var_input_val+"."+i->first+"_"+to_string(param2.pop1_lower)+"_"+to_string(param2.pop2_upper)+"_"+to_string(param2.min_sample)+"_valid_op"+splicing_status;
				output_vec.push_back(valid_op);
				variants_discovery_validation(var_input_disc,pop1,disc_op,var_input_val,pop2,valid_op,ann,param_disc,param_valid,var);
			}			
		}else{
			cout<<i->first<<" not found in the validation dataset, analysis will be skipped."<<endl;
		}
	}
	cout<<"============================================================================================="<<endl;
	cout<<"Merging output variant_list for "<<output_vec.size()<<" population."<<endl;
	var_output+=splicing_status;
	var.merge_sort_variant_list(output_vec,var_output);
	cout<<"Merging completed."<<endl<<endl;
	if(param2.output_per_sample_gt){
		cout<<"Merging output genotype list."<<endl;
		var_output+=".gt";
		var.merge_sort_variant_list(output_gt_vec,var_output);
	}
	cout<<"Genotype file merging completed."<<endl<<endl;
}

void main_analysis_module::bi_pop_freq_analysis(string var_input_disc,string var_input_val,string var_output,pop_data pop1,pop_data pop2,ann_data ann,input_param param1,input_param param2,var_list& var){
	input_param param_disc=param1,param_valid=param2;
	data_printer dp;
	string disc_op,valid_op,splicing_status;
	vector<string> pop_vec,pair_vec;
	set<string> pop_set;
	if(param1.no_splicing){
		splicing_status="_no_splicing";
	}else{
		splicing_status="";
	}

	for (unordered_map<string,set<string> >::iterator i = pop1.pop_dict.begin();i != pop1.pop_dict.end(); ++i){
		if(pop2.pop_col_dict.count(i->first)){
			if(pop2.pop_col_dict[i->first].size()>=param2.min_sample){
				pop_vec.push_back(i->first);
				pop_set.insert(i->first);
			}else{
				cout<<"Only "<<pop2.pop_col_dict[i->first].size()<<" sample from "<<i->first<<" were detected in the validation dataset, failed to meet minimum sample requirement. Analysis will be skipped."<<endl;				
			}
		}else{
			cout<<i->first<<" not found in the validation dataset, analysis will be skipped."<<endl;
		}
	}
	cout<<"A total of "<<pop_vec.size()<<" populations available for the bi-population analysis."<<endl;

	for (int i = 0; i < pop_vec.size()-1; ++i){
		for (int j = i+1; j < pop_vec.size(); ++j){
		cout<<"============================================================================================="<<endl;
		cout<<"Start analysis for "<<pop_vec[i]<<"-"<<pop_vec[j]<<"  population pair."<<endl;
		pair_vec.clear();
		pair_vec.push_back(pop_vec[i]);
		pair_vec.push_back(pop_vec[j]);
		param_disc.load_population_label(pair_vec,pop_set);
		param_valid.load_population_label(pair_vec,pop_set);
		disc_op=var_input_disc+"."+pop_vec[i]+"-"+pop_vec[j]+"_"+to_string(param1.pop1_lower)+"_"+to_string(param1.pop2_upper)+"_"+to_string(param1.min_sample)+"_disc_op"+splicing_status;
		valid_op=var_input_val+"."+pop_vec[i]+"-"+pop_vec[j]+"_"+to_string(param2.pop1_lower)+"_"+to_string(param2.pop2_upper)+"_"+to_string(param2.min_sample)+"_valid_op"+splicing_status;
		output_vec.push_back(valid_op);
		variants_discovery_validation(var_input_disc,pop1,disc_op,var_input_val,pop2,valid_op,ann,param_disc,param_valid,var);		
		}
	}

	cout<<"============================================================================================="<<endl;
	cout<<"Merge output variant_list for "<<output_vec.size()<<" population."<<endl;
	var_output=var_output+splicing_status;
	var.merge_sort_variant_list(output_vec,var_output);
	cout<<"Merge completed."<<endl<<endl;
	if(param2.output_per_sample_gt){
		cout<<"Merging output genotype list."<<endl;
		var_output+=".gt";
		var.merge_sort_variant_list(output_gt_vec,var_output);
	}
	cout<<"Genotype file merging completed."<<endl<<endl;
}

#endif