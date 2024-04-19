#include<chrono>
#include<mutex>
#include<algorithm>
#include<thread>
#include<queue>
#include<condition_variable>
#include<filesystem>
#include "utilities.h"
#include "pop_data.h"
#include "var_list.h"
#include "ann_data.h"
#include "input_param.h"
#include "main_analysis_module.h"

using namespace std;

bool g_ready = false;
condition_variable g_cv;
mutex ip_mutex;	
mutex op_mutex;

class thread_analysis_module{
	public:
		var_list var;
		pop_data pop;
		ann_data ann;
		input_param param;
		data_printer dp;

		void load_thread_vec(){
			cout<<th_num<<" cores detected."<<endl;
			for (int i = 0; i < th_num; ++i){
				cter_vec.push_back(0);
			}
		}

		void load_io_address(string str1, string str2){
			input_ad=str1;
			output_ad=str2;
			cout<<"Input file is: "<<input_ad<<endl;
			cout<<"Output files is: "<<output_ad<<endl;
			temp_ad=output_ad+".temp_file";
			cout<<"Temporary data will be stored at: "<<temp_ad<<endl<<endl;	
		}

		pair<int,string> extract_line_idx(string line){
			string temp_str="";
			int idx,s_idx;
			for (int i = 0; i < line.length(); ++i){
				if(line[i]=='\t'){
					idx=stoi(temp_str);
					s_idx=i;
					break;
				}
				temp_str+=line[i];
			}
			temp_str=line.substr(s_idx+1,line.length()-s_idx);
			return make_pair(idx,temp_str);
		}

		void insert_sort_variant_files(string ip_str,string op_str){
			string line,header_str;
			pair<int,string> temp_pair;
			vector<pair<int,string> > sorted_vec(var_num);
			ifstream input_sort;
			ofstream output_sort;
			input_sort.open(ip_str);
			check_file_open_status(input_sort,ip_str);
			int cur=0,cur_idx,temp_int;
			while(getline(input_sort,line)){
				if(line[0]=='#'){
					header_str=line;
				}else{
					temp_pair=extract_line_idx(line);
					cur_idx=temp_pair.first;
					if(cur==0){
						sorted_vec[cur]=temp_pair;
					}else{
						temp_int=cur-1;
						while(cur_idx<sorted_vec[temp_int].first&&temp_int>=0){
							sorted_vec[temp_int+1]=sorted_vec[temp_int];
							--temp_int;
						}
						sorted_vec[temp_int+1]=temp_pair;
					}
					cur++;
				}
			}
			input_sort.close();
			output_sort.open(op_str);
			check_file_open_status(output_sort,op_str);
			output_sort<<header_str<<endl;
			for (int i = 0; i < sorted_vec.size(); ++i){
				output_sort<<sorted_vec[i].second<<endl;
			}
			output_sort.close();
		}	

		unordered_map<string,set<int>> sample_from_valid_population(pop_data& pop,input_param& param){
			vector<int> subsample_idx_vec;
			for(unordered_map<string,set<int> >::iterator i=pop.pop_col_dict.begin();i!=pop.pop_col_dict.end();++i){
				if(i->second.size()>=param.min_sample){
					for (set<int>::iterator j = i->second.begin(); j != i->second.end(); ++j){
						subsample_idx_vec.push_back(*j);
					}
				}
			}
			pair<vector<int>,vector<int> > splitted_idx_vecs=split_vec_into_two(param.sample_frac,param.seed,subsample_idx_vec);
			subsample_idx_vec_disc=splitted_idx_vecs.first;
			subsample_idx_vec_valid=splitted_idx_vecs.second;
			//Load valid pop data
			unordered_set<string> temp_pop_set;
			unordered_map<string,set<int> > subsample_temp_pop_col_dict;
			for (int i = 0; i !=subsample_idx_vec_disc.size(); ++i){
				int sample_idx=subsample_idx_vec_disc[i];
				string sample_name=pop.col_sample_dict[sample_idx];
				string pop_name=pop.sample_dict[sample_name];
				if(!subsample_temp_pop_col_dict.count(pop_name)){
					set<int> temp_set;
					subsample_temp_pop_col_dict[pop_name]=temp_set;
				}
				subsample_temp_pop_col_dict[pop_name].insert(sample_idx);
			}
			unordered_map<string,set<int>> subsample_pop_col_dict_op;
			for(unordered_map<string,set<int> >::iterator i=subsample_temp_pop_col_dict.begin();i!=subsample_temp_pop_col_dict.end();++i){
				if(i->second.size()>=param.min_sample){
					subsample_pop_col_dict_op[i->first]=i->second;
				}
			}
			return subsample_pop_col_dict_op;
		}


		//Frequency analysis module

		void variant_data_reader(int j,pop_data& pop,input_param& param,var_list& var,ann_data& ann,int thread_num){
			input.open(input_ad);
			check_file_open_status(input,input_ad);
			output.open(temp_ad);
			check_file_open_status(output,output_ad);
			int line_cter=1,var_cter;
			string line,chr_str,pos_str,ann_str;
			main_analysis_module main;
			auto start=chrono::high_resolution_clock::now();
			main.output_file_header(output,var,param,pop,pop_vec);			
			cout<<endl<<"Output header generated, start analyzing variant data."<<endl;
			row_cur=0,var_num=0;
			int input_num=0;
			var_num=0;
			
			while(getline(input,line)){
				if(line[0]!='#'){
					input_num++;
					if(line_cter%100000==0){
						cout<<line_cter<<" variants loaded."<<endl;
					}
					line_cter++;
					var_cter=si_q.size();

					while(si_q.size()>=(thread_num-1)*300000){
						this_thread::sleep_for(chrono::milliseconds(1000));
						cout<<si_q.size()<<" in-memory variants currently waiting for processing..."<<endl;
						var_cter=si_q.size();
					}

					if(var.var_mat.size()==0){
						if(param.ann_flag.length()>0){
							ann_str=find_str_after_nth_char(line,7,'\t');
							if(find_str_after_nth_char(ann_str,1,'|').find(param.ann_flag)!=string::npos){
								pair<string,pair<int,int> > temp_pair=make_pair(line,make_pair(line_cter,row_cur));
								// pair<string,int> temp_pair=make_pair(line,line_cter);
								unique_lock<mutex> ul(ip_mutex);
								si_q.push(temp_pair);
								g_ready=true;
								ul.unlock();
								g_cv.notify_all();								
							}
						}else{
							pair<string,pair<int,int> > temp_pair=make_pair(line,make_pair(line_cter,row_cur));
							// pair<string,int > temp_pair=make_pair(line,line_cter);
							unique_lock<mutex> ul(ip_mutex);
							si_q.push(temp_pair);
							g_ready=true;
							ul.unlock();
							g_cv.notify_all();	
						}
					}else{
						if(row_cur==var.var_mat.size()){
							break;
						}

						chr_str=find_str_after_nth_char(line,0,'\t');
						pos_str=find_str_after_nth_char(line,1,'\t');

						int cur_chr_idx=stoi(chr_str.substr(7,2))-4,tar_chr_idx=stoi(var.var_mat[row_cur][0].substr(7,2))-4;

						while(cur_chr_idx>tar_chr_idx){
							row_cur++;
							if(row_cur==var.var_mat.size()){
								break;
							}
							tar_chr_idx=stoi(var.var_mat[row_cur][0].substr(7,2))-4;

						}
						if(row_cur==var.var_mat.size()){
							break;
						}
						
						if(var.var_mat[row_cur][0]==chr_str){
							int cur_pos=stoi(pos_str),tar_pos=stoi(var.var_mat[row_cur][1]);
							while(cur_pos>tar_pos&&var.var_mat[row_cur][0]==chr_str){
								row_cur++;
								if(row_cur==var.var_mat.size()){
									break;
								}
								tar_pos=stoi(var.var_mat[row_cur][1]);
							}
							if(row_cur==var.var_mat.size()){
								break;
							}

							if(var.var_mat[row_cur][1]==pos_str){
								if(param.ann_flag.length()>0){
									ann_str=find_str_after_nth_char(line,7,'\t');
									if(find_str_after_nth_char(ann_str,1,'|').find(param.ann_flag)!=string::npos){
										pair<string,pair<int,int> > temp_pair=make_pair(line,make_pair(line_cter,row_cur));
										unique_lock<mutex> ul(ip_mutex);
										si_q.push(temp_pair);
										g_ready=true;
										ul.unlock();
										g_cv.notify_all();										
									}
								}else{
									pair<string,pair<int,int> > temp_pair=make_pair(line,make_pair(line_cter,row_cur));
									unique_lock<mutex> ul(ip_mutex);
									si_q.push(temp_pair);
									g_ready=true;
									ul.unlock();
									g_cv.notify_all();	
								}
							}
						}						
					}
				}
			}
			cout<<"All "<<input_num<<" variants scanned. Waiting for frequency calculation to be completed..."<<endl;
			var_cter=si_q.size();
			cout<<var_cter<<" variants remained for processing..."<<endl;

			while(si_q.size()!=0){
				this_thread::sleep_for(chrono::milliseconds(10000));
				cout<<si_q.size()<<" variants remained for processing...Processing speed "<<(var_cter-si_q.size())*6<<" variants/min."<<endl;
				var_cter=si_q.size();
			}

			input.close();
			cout<<"Frequency calculation completed. ";
			auto end=chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(end - start);
			cout<<"Time lapsed: "<<duration.count()<<" seconds."<<endl<<endl;
			cout<<"Clearing up all cached queries..."<<endl<<endl;
			this_thread::sleep_for(chrono::milliseconds(3000));
			completed=true;
			g_cv.notify_all();
			g_ready=true;
			output.close();
			cout<<"A total of "<<var_num<<"  variant output."<<endl;
			cout<<"Sorting variant data..."<<endl;
			insert_sort_variant_files(temp_ad,output_ad);
			filesystem::remove(temp_ad);
			cout<<"Sorting completed."<<endl;
		}

		void freq_data_loader(int i,pop_data& pop,input_param& param,var_list& var,ann_data& ann){
			string ip_str,freq_str,gt_str;
			int row_mark,line_mark;
			vector<string> line_vec,pop_rst_vec;
			pair<string,pair<int,int> > q_pair;
			pair<vector<string>,vector<string> > str_vec_pair;
			while(true){
				if(completed){
					break;
				}
				unique_lock<mutex> ul(ip_mutex);	
				if(si_q.size()==0){
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=si_q.front();
					si_q.pop();
					ul.unlock();
					//main analysis
					ip_str=q_pair.first;
					line_mark=q_pair.second.first;
					row_mark=q_pair.second.second;
					main_analysis_module main;
					line_vec=read_char_delim_str(ip_str,'\t');
					vector<int> ind_vec;
					if(line_vec[4].length()!=1){
						goto labelA;
					}
					if(!param.no_splicing){
						str_vec_pair=analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
					}else{
						if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1){
							str_vec_pair=analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
						}
					}
					freq_str=str_vec_pair.first[0];
					gt_str=str_vec_pair.first[1];
					pop_rst_vec=str_vec_pair.second;
					//output data
					if(param.dist_mode){
						if(ind_vec[7]==1){
							ul.lock();
							output<<line_mark<<'\t'<<freq_str<<endl;
							var_num++;
							g_ready=false;
							ul.unlock();							
						}
					}else{
						if(param.exhaust_valid_mode){
							if (ind_vec[0]!=0&&ind_vec[4]==0){
								if(ind_vec[6]==0){
									ul.lock();
									output<<line_mark<<'\t'<<freq_str<<endl;
									var_num++;
									g_ready=false;
									ul.unlock();
								}
							}
						}else{
							if(ind_vec[0]*ind_vec[1]!=0&&ind_vec[0]==param.pop1.size()&&ind_vec[4]==0){
								if(ind_vec[2]*ind_vec[3]!=0){
									ul.lock();
									output<<line_mark<<'\t'<<freq_str<<endl;
									var_num++;
									g_ready=false;
									ul.unlock();
								}
							}
						}
					}
					labelA:;
				}
			}
			completed=true;
			g_cv.notify_all();
			g_ready=true;			
		}

		pair<vector<string>,vector<string> > analyze_pop_freq_per_line(pop_data& pop,input_param& param,var_list& var,ann_data& ann,vector<string>& line_vec,vector<string>& pop_vec,int row_mark,vector<int>& ind_vec){
			int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0,ind6=0;
			int total_num,case_num,var_allele_num;
			int site_total_ind_num=0,site_var_ind_num=0,site_var_num=0;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			unordered_map<string,string> freq_dict;
			string temp_str="",op_str="";
			for (int i = 0; i < pop_vec.size(); ++i){
				total_num=0;
				case_num=0;
				var_allele_num=0;

				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j){
					if(param.depth_file.length()!=0){
						if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j]))){
							total_num+=1;
							int temp;
							if(param.flip){
								temp=check_genotype_ref(line_vec[*j]);
							}else{
								temp=check_genotype(line_vec[*j]);
							}
							case_num+=temp;
						}
					}else{
						if(check_read_depth(line_vec[*j],param.min_depth)){
							total_num+=1;
							int temp;
							if(param.flip){
								temp=check_genotype_ref(line_vec[*j]);
							}else{
								temp=check_genotype(line_vec[*j]);
							}
							if(line_vec[*j][0]=='1'){
								var_allele_num+=1;
							}
							if(line_vec[*j][2]=='1'){
								var_allele_num+=1;
							}
							case_num+=temp;
						}
					}
				}

				site_total_ind_num+=total_num;
				site_var_ind_num+=case_num;
				site_var_num+=var_allele_num;

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
									for (int k = 0; k < var.var_mat[row_mark].size(); ++k){
										if(var.var_mat[row_mark][k]==pop_vec[i]){
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
					if(total_num==0){
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

			double var_freq=(double)site_var_num/(2*site_total_ind_num);
			if(var_freq>param.af&&var_freq<(1-param.af)){
				ind6=1;
			}

			string hline_str=generate_line_header(line_vec);
			op_str=hline_str;
			if(param.dist_mode){
				op_str+=to_string(site_var_num)+'\t'+to_string(site_var_ind_num)+'\t'+to_string(site_total_ind_num);
			}
			for (int i = 0; i < pop_vec.size(); ++i){
				op_str+='\t'+freq_dict[pop_vec[i]];
			}
			//Store gt_data
			//Updating pop_vec	
			if(var.var_mat.size()!=0){
				temp_vec.push_back(var.var_mat[row_mark][0]);
				temp_vec.push_back(var.var_mat[row_mark][1]);
				for (int i = 0; i < var_pop_vec.size(); ++i){
					temp_vec.push_back(var_pop_vec[i]);
				}				
			}		

			vector<string> temp_vec2;
			temp_vec2.push_back(op_str);
			temp_vec2.push_back(temp_str);
			//Store ind variable
			ind_vec.push_back(valid_ind1);
			ind_vec.push_back(valid_ind2);
			ind_vec.push_back(ind1);
			ind_vec.push_back(ind2);
			ind_vec.push_back(ind3);
			ind_vec.push_back(ind4);
			ind_vec.push_back(ind5);
			ind_vec.push_back(ind6);
			return make_pair(temp_vec2,temp_vec);
		}

		//Validation module

		void ev_variant_reader(int j,pop_data& pop,input_param& param,var_list& var,ann_data& ann){
			input.open(input_ad);
			check_file_open_status(input,input_ad);
			output.open(temp_ad);
			check_file_open_status(output,output_ad);
			int line_cter=1,var_cter;
			string line,chr_str,pos_str,ann_str;
			main_analysis_module main;
			auto start=chrono::high_resolution_clock::now();
			main.output_file_header(output,var,param,pop,pop_vec);			
			cout<<endl<<"Output header generated, start analyzing variant data."<<endl;
			row_cur=0;
			int input_num=0;
			while(getline(input,line)){
				if(line[0]!='#'){
					input_num++;
					if(line_cter%100000==0){
						cout<<line_cter<<" variants loaded."<<endl;
					}
					line_cter++;
					if(row_cur==var.var_mat.size()){
						break;
					}
					chr_str=find_str_after_nth_char(line,0,'\t');
					pos_str=find_str_after_nth_char(line,1,'\t');
					int cur_chr_idx=stoi(chr_str.substr(7,2))-4;
					int tar_chr_idx=stoi(var.var_mat[row_cur][0]);
					while(cur_chr_idx>tar_chr_idx){
						row_cur++;
						if(row_cur==var.var_mat.size()){
							break;
						}
						tar_chr_idx=stoi(var.var_mat[row_cur][0]);
					}
					if(row_cur==var.var_mat.size()){
						break;
					}	
					if(stoi(var.var_mat[row_cur][0])==cur_chr_idx){
						int cur_pos=stoi(pos_str),tar_pos=stoi(var.var_mat[row_cur][1]);
						while(cur_pos>tar_pos&&stoi(var.var_mat[row_cur][0])==cur_chr_idx){
							row_cur++;
							if(row_cur==var.var_mat.size()){
									break;
							}
							tar_pos=stoi(var.var_mat[row_cur][1]);
						}
						if(row_cur==var.var_mat.size()){
							break;
						}
						if(var.var_mat[row_cur][1]==pos_str){
							if(param.ann_flag.length()>0){
								ann_str=find_str_after_nth_char(line,7,'\t');
								if(find_str_after_nth_char(ann_str,1,'|').find(param.ann_flag)!=string::npos){
									pair<string,pair<int,int> > temp_pair=make_pair(line,make_pair(line_cter,row_cur));
									unique_lock<mutex> ul(ip_mutex);
									si_q.push(temp_pair);
									g_ready=true;
									ul.unlock();
									g_cv.notify_all();										
								}
							}else{
								pair<string,pair<int,int> > temp_pair=make_pair(line,make_pair(line_cter,row_cur));
								unique_lock<mutex> ul(ip_mutex);
								si_q.push(temp_pair);
								g_ready=true;
								ul.unlock();
								g_cv.notify_all();	
							}
						}
					}	
				}
			}
			cout<<"All "<<input_num<<" variants scanned. Waiting for frequency calculation to be completed..."<<endl;
			var_cter=si_q.size();
			cout<<var_cter<<" variants remained for processing..."<<endl;
			while(si_q.size()!=0){
				this_thread::sleep_for(chrono::milliseconds(10000));
				cout<<si_q.size()<<" variants remained for processing...Processing speed "<<(var_cter-si_q.size())*6<<" variants/min."<<endl;
				var_cter=si_q.size();
			}
			input.close();
			cout<<"Frequency calculation completed. ";
			auto end=chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(end - start);
			cout<<"Time lapsed: "<<duration.count()<<" seconds."<<endl<<endl;
			completed=true;
			g_cv.notify_all();
			g_ready=true;
			output.close();
			cout<<"A total of "<<var_num<<" target variant found."<<endl;
			cout<<"Sorting variant data..."<<endl;
			insert_sort_variant_files(temp_ad,output_ad);
			filesystem::remove(temp_ad);
			cout<<"Sorting completed."<<endl;
		}

		void ev_freq_loader(int i,pop_data& pop,input_param& param,var_list& var,ann_data& ann){
			string ip_str,freq_str,gt_str;
			int row_mark,line_mark;
			vector<string> line_vec,pop_rst_vec;
			pair<string,pair<int,int> > q_pair;
			pair<vector<string>,vector<string> > str_vec_pair;
			while(true){
				if(completed){
					break;
				}
				unique_lock<mutex> ul(ip_mutex);	
				if(si_q.size()==0){
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=si_q.front();
					si_q.pop();
					ul.unlock();
					//main analysis
					ip_str=q_pair.first;
					line_mark=q_pair.second.first;
					row_mark=q_pair.second.second;
					main_analysis_module main;
					line_vec=read_char_delim_str(ip_str,'\t');
					vector<int> ind_vec;				
					if(!param.no_splicing){
						str_vec_pair=ev_analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
					}else{
						if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1){
							str_vec_pair=ev_analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
						}
					}
					freq_str=str_vec_pair.first[0];
					gt_str=str_vec_pair.first[1];
					//output data
					if(param.exhaust_valid_mode)
					{
						if (ind_vec[0]!=0&&ind_vec[4]==0)
						{
							if(ind_vec[6]==0)
							{
								ul.lock();
								output<<line_mark<<'\t'<<freq_str<<endl;
								var_num++;
								g_ready=false;
								ul.unlock();
							}
						}						
					}else if(param.unique_valid_mode)
					{
						if (ind_vec[0]!=0&&ind_vec[4]==0)
						{
							if (ind_vec[2]==0&&ind_vec[6]==0&&ind_vec[3]==1)
							{
								ul.lock();
								output<<line_mark<<'\t'<<freq_str<<endl;
								var_num++;
								g_ready=false;
								ul.unlock();
							}
						}				
					}
				}
			}
		}

		pair<vector<string>,vector<string> > ev_analyze_pop_freq_per_line(pop_data& pop,input_param& param,var_list& var,ann_data& ann,vector<string>& line_vec,vector<string>& pop_vec,int row_mark,vector<int>& ind_vec){
			int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0;
			int total_num,case_num,total_num2,case_num2;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			unordered_map<string,string> freq_dict;
			string temp_str="",op_str="";
			if(param.exhaust_valid_mode)
			{
				for (int i = 0; i < pop_vec.size(); ++i)
				{
					total_num=0;
					case_num=0;
					for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j)
					{
						if(param.depth_file.length()!=0)
						{
							if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j])))
							{
								total_num+=1;
								int temp=check_genotype(line_vec[*j]);
								case_num+=temp;
							}
						}else{
							if(check_read_depth(line_vec[*j],param.min_depth))
							{
								total_num+=1;
								int temp=check_genotype(line_vec[*j]);
								case_num+=temp;
							}
						}
					}
					if(total_num>=param.min_sample)
					{
						freq=(double)case_num/(double)total_num;
						freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
						if(param.pop1.count(pop_vec[i]))
						{
							valid_ind1+=1;
							if (!(freq>=param.pop1_lower&&freq<=param.pop1_upper))
							{
								ind1*=0;
							}else{
								for (int k = 2; k < var.var_mat[row_mark].size(); ++k)
								{
									if(var.var_mat[row_mark][k]==pop_vec[i])
									{
										ind5*=0;
									}								
								}
							}
						}else{
							if(param.pop2.count(pop_vec[i]))
							{
								valid_ind2+=1;
								if(!(freq>=param.pop2_lower&&freq<=param.pop2_upper))
								{
									ind2*=0;
								}
							}
						}
					}else{
						if (total_num==0)
						{
							freq_dict[pop_vec[i]]="0/0";
						}else{
							freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
						}
					}
				}
			}else if(param.unique_valid_mode){
				int pop1_num=0,pop_neither_num=0;
				int cumu_total=0,cumu_case=0;
				for(int i = 0; i < pop_vec.size(); ++i)
				{
					total_num=0;
					case_num=0;
					for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j)
					{
						if(param.depth_file.length()!=0)
						{
							if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j])))
							{
								total_num+=1;
								int temp=check_genotype(line_vec[*j]);
								case_num+=temp;
							}
						}else{
							if(check_read_depth(line_vec[*j],param.min_depth))
							{
								total_num+=1;
								int temp=check_genotype(line_vec[*j]);
								case_num+=temp;
							}
						}
					}
					cumu_total+=total_num;
					cumu_case+=case_num;

					if(param.min_sample_ref==-1||param.min_sample_tar==-1)
					{
						if(total_num>=param.min_sample)
						{
							valid_ind1++;
							freq=(double)case_num/(double)total_num;
							freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
							if (freq>=param.pop1_lower&&freq<=param.pop1_upper)
							{
								pop1_num+=1;
								for (int a=2;a!=var.var_mat[row_mark].size();++a)
								{
									if(pop_vec[i]==var.var_mat[row_mark][a])
									{
										ind5*=0;
									}
								}
							}else{
								if (!(freq>=param.pop2_lower&&freq<=param.pop2_upper))
								{
									pop_neither_num+=1;
								}
							}									
						}else{
							if (total_num==0)
							{
								freq_dict[pop_vec[i]]="0/0";
							}else{
								freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
							}
						}
					}else{
						//Check populations to be validated
						set<string> st;
						for (int a=2;a!=var.var_mat[row_mark].size();++a)
						{
							st.insert(var.var_mat[row_mark][a]);
						}
						//Check for target population
						if(st.count(pop_vec[i]))
						{
							if (total_num>=param.min_sample_tar)
							{
								valid_ind1++;
								freq=(double)case_num/(double)total_num;
								freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
								if (freq>=param.pop1_lower&&freq<=param.pop1_upper)
								{
									pop1_num+=1;
									ind5*=0;
								}								
							}else{
								if (total_num==0)
								{
									freq_dict[pop_vec[i]]="0/0";
								}else{
									freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
								}
							}
						}else{
							if(total_num>=param.min_sample_ref)
							{
								valid_ind1++;
								freq=(double)case_num/(double)total_num;
								freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
								if (!(freq>=param.pop2_lower&&freq<=param.pop2_upper))
								{
									pop_neither_num+=1;
								}
							}else{
								if (total_num==0)
								{
									freq_dict[pop_vec[i]]="0/0";
								}else{
									freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
								}								
							}
						}						
					}
				}

				if(ind5==0)
				{
					if(pop1_num==param.pop_num&&pop_neither_num==0)
					{
						ind1*=0;
						string hline_str=generate_line_header(line_vec);
						op_str=hline_str;
						for (int i = 0; i < pop_vec.size(); ++i)
						{
							op_str+='\t'+freq_dict[pop_vec[i]];
						}
					}					
				}
			}

			if (total_num!=case_num)
			{
				ind3*=0;
			}

			if(ind1==0)
			{
				op_str=var.var_mat[row_mark][0]+'\t'+var.var_mat[row_mark][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10];
				for (int i = 0; i < pop_vec.size(); ++i)
				{
					op_str+='\t'+freq_dict[pop_vec[i]];
				}				
			}
			//Store gt_data
			temp_str=var.var_mat[row_mark][0]+'\t'+var.var_mat[row_mark][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10]+temp_str;
			//Updating pop_vec			
			temp_vec.push_back(var.var_mat[row_mark][0]);
			temp_vec.push_back(var.var_mat[row_mark][1]);
			for (int i = 0; i < var_pop_vec.size(); ++i)
			{
				temp_vec.push_back(var_pop_vec[i]);
			}
			vector<string> temp_vec2;
			temp_vec2.push_back(op_str);
			temp_vec2.push_back(temp_str);
			//Store ind variable
			ind_vec.push_back(valid_ind1);
			ind_vec.push_back(valid_ind2);
			ind_vec.push_back(ind1);
			ind_vec.push_back(ind2);
			ind_vec.push_back(ind3);
			ind_vec.push_back(ind4);
			ind_vec.push_back(ind5);
			return make_pair(temp_vec2,temp_vec);
		}

		//Exhaustive discovery module

		void ed_variant_reader(int j,pop_data& pop,input_param& param,var_list& var,ann_data& ann,int thread_num){
			input.open(input_ad);
			check_file_open_status(input,input_ad);
			output.open(temp_ad);
			check_file_open_status(output,output_ad);
			//Check for sampling option
			if(param.sample_frac<1)
			{
				subsample_pop_col_dict=sample_from_valid_population(pop,param);
			}
			int line_cter=1,var_cter;
			string line,chr_str,pos_str,ann_str;
			main_analysis_module main;
			auto start=chrono::high_resolution_clock::now();
			int valid_subpopulation_num=0;
			if(param.sample_frac<1)
			{
				if(param.analysis_mode==11){
					main.output_file_header(output,var,param,pop,pop_vec);
				}else{
					output<<"#Chromosome"<<'\t'<<"Position"<<'\t'<<"Gene"<<'\t'<<"Variant Type"<<'\t'<<"Mutation Position"<<'\t'<<"Amino Acid Change";
					for (unordered_map<string,set<int> >::iterator i =  subsample_pop_col_dict.begin(); i != subsample_pop_col_dict.end(); ++i)
					{
						if(i->second.size()>=param.min_sample)
						{
							output<<'\t'<<i->first;
							pop_vec.push_back(i->first);
							valid_subpopulation_num++;
						}
					}
					output<<endl;
					cout<<"Number of valid subpopulations:"<<valid_subpopulation_num<<endl;
				}
			}else{
				main.output_file_header(output,var,param,pop,pop_vec);
			}
			if(param.likelihood_file!="")
			{
				output2.open(param.likelihood_file);
				//output likelihood file header
				output2<<"#Chromosome"<<'\t'<<"Position";
				for (int i = 0; i < pop_vec.size(); ++i)
				{
					output2<<'\t'<<pop_vec[i];
				}
				output2<<endl;
			}
			cout<<endl<<"Output header generated, start analyzing variant data."<<endl;
			var_num=0;
			while(getline(input,line))
			{
				if(line[0]!='#')
				{
					if(line_cter%100000==0)
					{
						cout<<line_cter<<" total variants scanned."<<endl;
					}

					line_cter++;
					var_cter=sp_q.size();
					while(sp_q.size()>=(thread_num-1)*200000)
					{
						this_thread::sleep_for(chrono::milliseconds(1000));
						cout<<sp_q.size()<<" in-memory variants currently waiting for processing..."<<endl;
						//cout<<"Processing speed "<<(var_cter-sp_q.size())*60<<" variants/min."<<endl;
						var_cter=sp_q.size();
					}
					if(param.ann_flag.length()>0)
					{
						ann_str=find_str_after_nth_char(line,7,'\t');
						if(find_str_after_nth_char(ann_str,1,'|').find(param.ann_flag)!=string::npos)
						{
							pair<string,int> temp_pair=make_pair(line,line_cter);
							unique_lock<mutex> ul(ip_mutex);
							sp_q.push(temp_pair);
							g_ready=true;
							ul.unlock();
							g_cv.notify_all();								
						}
					}else{
						pair<string,int > temp_pair=make_pair(line,line_cter);
						unique_lock<mutex> ul(ip_mutex);
						sp_q.push(temp_pair);
						g_ready=true;
						ul.unlock();
						g_cv.notify_all();	
					}
				}
			}
			input.close();
			cout<<"All "<<line_cter-1<<" variants scanned. Waiting for frequency calculation to be completed..."<<endl;
			var_cter=sp_q.size();
			cout<<var_cter<<" variants remained for processing..."<<endl;
			while(sp_q.size()!=0)
			{
				this_thread::sleep_for(chrono::milliseconds(10000));
				cout<<sp_q.size()<<" variants remained for processing...Processing speed "<<(var_cter-sp_q.size())*6<<" variants/min."<<endl;
				var_cter=sp_q.size();
			}
			cout<<"Frequency calculation completed. ";
			auto end=chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(end - start);
			cout<<"Time lapsed: "<<duration.count()<<" seconds."<<endl<<endl;
			completed=true;
			g_cv.notify_all();
			g_ready=true;
			output.close();
			if(param.analysis_mode==3||param.analysis_mode==1||param.analysis_mode==5)
			{
				var.sort_var_pop_mat();
				var.output_var_pop_file(param.output_file+".pop_list");
				cout<<endl;
			}
			if(param.analysis_mode!=10&&param.analysis_mode!=11)
			{
				cout<<"A total of "<<eff_var_num<<" effective variant scanned."<<endl;
			}
			
			cout<<"A total of "<<var_num<<"  variant output."<<endl;
			if (param.sample_frac<1)
			{
				cout<<"Output analysis specs for the sampled subset."<<endl;
				ofstream temp_output;
				temp_output.open(output_ad+".sample_partition_list");
				temp_output<<"#Breed\tNumber of sampled sample\tNumber of left-over samples"<<endl;
				// cout<<subsample_pop_col_dict.size()<<endl;
				for (unordered_map<string,set<int>>::iterator i = subsample_pop_col_dict.begin(); i !=subsample_pop_col_dict.end(); ++i)
				{
					temp_output<<i->first<<'\t'<<i->second.size()<<'\t'<<pop.pop_col_dict[i->first].size()-i->second.size()<<endl;
				}
				temp_output<<'*'<<eff_var_num<<endl;
				temp_output.close();
			}
			//Likelihood output stop
			if(param.likelihood_file!="")
			{
				output2.close();
			}
			cout<<"Sorting variant data..."<<endl;
			insert_sort_variant_files(temp_ad,output_ad);
			filesystem::remove(temp_ad);
			cout<<"Sorting completed."<<endl;
		}

		void ed_freq_data_loader(int i,pop_data& pop,input_param& param,var_list& var,ann_data& ann){
			string ip_str,freq_str,gt_str;
			int row_mark,line_mark;
			vector<string> line_vec,pop_rst_vec;
			pair<string,int > q_pair;
			pair<vector<string>,vector<string> > str_vec_pair;
			while(true)
			{
				if(completed)
				{
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0)
				{
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=sp_q.front();
					// cter_vec[i]++;
					sp_q.pop();
					ul.unlock();
					//main analysis
					ip_str=q_pair.first;
					line_mark=q_pair.second;
					main_analysis_module main;
					line_vec=read_char_delim_str(ip_str,'\t');
					vector<int> ind_vec;
					if(!param.INDEL_mode&&line_vec[4].length()!=1)
					{
						goto labelB;
					}
					if(!param.no_splicing)
					{
						str_vec_pair=ed_analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
					}else{
						if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1)
						{
							str_vec_pair=ed_analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
						}
					}
					freq_str=str_vec_pair.first[0];
					//output data
						if(ind_vec[5]==0&&ind_vec[4]==0)
						{
							if(param.analysis_mode==3||param.analysis_mode==1)
							{
								if(ind_vec[0]<=param.max_homo_pop)
								{
									ul.lock();
									output<<line_mark<<'\t'<<freq_str<<endl;
									var_num++;
									g_ready=false;
									ul.unlock();
								}
							}else{
								ul.lock();
								output<<line_mark<<'\t'<<freq_str<<endl;
								var_num++;
								g_ready=false;
								ul.unlock();
							}
						}
					labelB:;
				}
			}
			completed=true;
			g_cv.notify_all();
			g_ready=true;
		}

		pair<vector<string>,vector<string> > ed_analyze_pop_freq_per_line(pop_data& pop,input_param& param,var_list& var,ann_data& ann,vector<string>& line_vec,vector<string>& pop_vec,int row_mark,vector<int>& ind_vec){
			int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0;
			int total_num,case_num,total_num2,case_num2;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			unordered_map<string,string> freq_dict;
			string gt_str="",op_str="";
			for (int i = 0; i < pop_vec.size(); ++i)
			{
				total_num=0;
				case_num=0;
				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j)
				{
					if(param.depth_file.length()!=0)
					{
						if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j])))
						{
							total_num+=1;
						int temp;
						if(param.flip)
						{
							temp=check_genotype_ref(line_vec[*j]);
						}else{
							temp=check_genotype(line_vec[*j]);
						}
							case_num+=temp;
						}
					}else{
						if(check_read_depth(line_vec[*j],param.min_depth))
						{
							total_num+=1;
							int temp;
						if(param.flip)
						{
							temp=check_genotype_ref(line_vec[*j]);
						}else{
							temp=check_genotype(line_vec[*j]);
						}
							case_num+=temp;
						}
					}
				}
				if(total_num>=param.min_sample)
				{
					freq=(double)case_num/(double)total_num;
					freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
					if (freq>=param.pop1_lower&&freq<=param.pop1_upper)
					{
						var_pop_vec.push_back(pop_vec[i]);
						valid_ind1+=1;
						ind4*=0;
					}
				}else{
					if (total_num==0)
					{
						freq_dict[pop_vec[i]]="0/0";
					}else{
						freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
					}
				}
				if (total_num==case_num)
				{
					ind3*=1;
				}else{
					ind3*=0;
				}
			}
			string hline_str=generate_line_header(line_vec);
			op_str=hline_str;
			for (int i = 0; i < pop_vec.size(); ++i)
			{
				op_str+='\t'+freq_dict[pop_vec[i]];
			}
			//Updating pop_vec		
			vector<string> temp_vec2;
			temp_vec2.push_back(op_str);
			temp_vec2.push_back(gt_str);
			//Store ind variable
			ind_vec.push_back(valid_ind1);
			ind_vec.push_back(valid_ind2);
			ind_vec.push_back(ind1);
			ind_vec.push_back(ind2);
			ind_vec.push_back(ind3);
			ind_vec.push_back(ind4);
			ind_vec.push_back(ind5);
			return make_pair(temp_vec2,temp_vec);
		}

		//Shared pattern discovery module 

		void unipop_data_loader(pop_data& pop,input_param& param,var_list& var,ann_data& ann){
			string ip_str,freq_str,gt_str,likelihood_str;
			int line_mark;
			vector<string> line_vec,pop_rst_vec;
			pair<string,int > q_pair;
			pair<vector<string>,vector<string> > str_vec_pair;
			while(true)
			{
				if(completed)
				{
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0)
				{
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=sp_q.front();
					sp_q.pop();
					ul.unlock();

					//main analysis
					ip_str=q_pair.first;
					line_mark=q_pair.second;
					line_vec=read_char_delim_str(ip_str,'\t');
					vector<int> ind_vec;
					if(!param.INDEL_mode&&line_vec[4].length()!=1)
					{
						goto labelA;
					}
					if(param.analysis_mode==5)
					{
						if(!param.no_splicing)
						{
							str_vec_pair=discover_unique_variant(pop,param,var,ann,line_vec,pop_vec,ind_vec);
						}else{
							if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1)
							{
								str_vec_pair=discover_unique_variant(pop,param,var,ann,line_vec,pop_vec,ind_vec);
							}
						}						
					}
					freq_str=str_vec_pair.first[0];
					
					//output data
					ul.lock();
					if(ind_vec[3]==0)
					{
						eff_var_num++;
					}
					//output likelihood
					if(param.likelihood_file!=""&&str_vec_pair.second.size()!=0)
					{
						output2<<str_vec_pair.second[0]<<endl;
					}
					//output freq
					if(ind_vec[2]==0&&ind_vec[4]==0)
					{
						output<<line_mark<<'\t'<<freq_str<<endl;
						var_num++;
					}
					g_ready=false;
					ul.unlock();
					labelA:;
				}
			}
			completed=true;
			g_cv.notify_all();
			g_ready=true;
		}

		pair<vector<string>,vector<string>> discover_unique_variant(pop_data& pop,input_param& param,var_list& var,ann_data& ann,vector<string>& line_vec,vector<string>& pop_vec,vector<int>& ind_vec){
			int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0;
			int total_num,case_num,total_num2,case_num2,pop1_num=0,pop_neither_num=0;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			vector<double> freq_ratio_vec;
			unordered_map<string,string> freq_dict;
			string gt_str="",op_str="";
			int cumu_total=0,cumu_case=0;
			int valid_pop_num=0;
			double max_freq=-1;
			string max_pop="";
			if(param.sample_frac==1)
			{
				//Analyze with full dataset
				for (int i = 0; i < pop_vec.size(); ++i)
				{
					total_num=0;
					case_num=0;
					for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j)
					{
						if(param.depth_file.length()!=0)
						{
							if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j])))
							{
								total_num+=1;
								int temp;
								if(param.flip){
									temp=check_genotype_ref(line_vec[*j]);
									if(pop.pop_col_dict[pop_vec[i]].size()>=3)
									{
										if(line_vec[*j][0]=='0')
										{
											ind2*=0;
										}
									}
								}else{
									temp=check_genotype(line_vec[*j]);
									if(pop.pop_col_dict[pop_vec[i]].size()>=3)
									{
										if(line_vec[*j][2]=='1')
										{
											ind2*=0;
										}
									}
								}
								case_num+=temp;
							}
						}else{
							if(check_read_depth(line_vec[*j],param.min_depth))
							{
								total_num+=1;
								int temp;
								if(param.flip)
								{
									temp=check_genotype_ref(line_vec[*j]);
									if(pop.pop_col_dict[pop_vec[i]].size()>=3)
									{
										if(line_vec[*j][0]=='0')
										{
											ind2*=0;
										}
									}
								}else{
									temp=check_genotype(line_vec[*j]);
									if(pop.pop_col_dict[pop_vec[i]].size()>=3)
									{
										if(line_vec[*j][2]=='1')
										{
											ind2*=0;
										}
									}
								}
								case_num+=temp;
							}
						}
					}
					cumu_total+=total_num;
					cumu_case+=case_num;
					if(total_num>=param.min_sample)
					{
						freq=(double)case_num/(double)total_num;
						freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
						if(!param.combine_all)
						{
							if (freq>=param.pop1_lower&&freq<=param.pop1_upper)
							{
								if(freq>max_freq){
									max_freq=freq;
									max_pop=pop_vec[i];
								}
								var_pop_vec.push_back(pop_vec[i]);
								pop1_num+=1;
							}else{
								if (!(freq>=param.pop2_lower&&freq<=param.pop2_upper))
								{
									pop_neither_num+=1;
								}
							}		

						}
					}else{
						if (total_num==0)
						{
							freq_dict[pop_vec[i]]="0/0";
						}else{
							freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
						}
					}
					// cout<<total_num<<'\t'<<case_num<<endl;	
					if (total_num!=case_num)
					{
						ind3*=0;
					}
				}
			}else{
				//Analyze with sampled dataset
				freq_ratio_vec.clear();
				for (int i = 0; i < pop_vec.size(); ++i)
				{
					if(subsample_pop_col_dict.count(pop_vec[i]))
					{
						total_num=0;
						case_num=0;
						for (set<int>::iterator j=subsample_pop_col_dict[pop_vec[i]].begin(); j!=subsample_pop_col_dict[pop_vec[i]].end(); ++j)
						{
							if(param.depth_file.length()!=0)
							{
								if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j])))
								{
									total_num+=1;
									int temp;
									if(param.flip)
									{
										temp=check_genotype_ref(line_vec[*j]);
										if(line_vec[*j][0]=='0')
										{
											ind2*=0;
										}
									}else{
										temp=check_genotype(line_vec[*j]);
										if(line_vec[*j][2]=='1')
										{
											ind2*=0;
										}
									}
									case_num+=temp;
								}
							}else{
								if(check_read_depth(line_vec[*j],param.min_depth))
								{
									total_num+=1;
									int temp;
									if(param.flip)
									{
										temp=check_genotype_ref(line_vec[*j]);
										if(line_vec[*j][0]=='0')
										{
											ind2*=0;
										}
									}else{
										temp=check_genotype(line_vec[*j]);
										if(line_vec[*j][2]=='1')
										{
											ind2*=0;
										}
									}
									case_num+=temp;
								}
							}
						}

						if(total_num>=param.min_sample)
						{
							cumu_total+=total_num;
							cumu_case+=case_num;							
						}
						if(total_num>=param.min_sample)
						{
							valid_pop_num++;
							freq=(double)case_num/(double)total_num;
							freq_ratio_vec.push_back(freq);
							freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
							if(!param.combine_all)
							{
								if (freq>=param.pop1_lower&&freq<=param.pop1_upper)
								{
									if(freq>max_freq){
										max_freq=freq;
										max_pop=pop_vec[i];
									}
									var_pop_vec.push_back(pop_vec[i]);
									pop1_num+=1;
								}else{
									if (!(freq>=param.pop2_lower&&freq<=param.pop2_upper))
									{
										pop_neither_num+=1;
									}
								}						
							}
						}else{
							freq_ratio_vec.push_back(-1);
							if (total_num==0)
							{
								freq_dict[pop_vec[i]]="0/0";
							}else{
								freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
							}
						}
						if (total_num!=case_num)
						{
							ind3*=0;
						}
					}
				}
			}
			string hline_str=generate_line_header(line_vec);
			if(param.flip)
			{
				hline_str+="|ref";
			}
			//Calculate likelihood ratio for all valid sample
			bool subsample_valid=false;
			vector<double> loglikelihood_vec;
			//Take the average of population number 
			//Output line only when subsample_valid=true
			if(param.likelihood_file!="")
			{
				for (int i = 0; i !=freq_ratio_vec.size(); ++i)
				{
					if(freq_ratio_vec[i]!=-1)//Traverse valid populations only
					{
						double loglikelihood=log(freq_ratio_vec[i]);
						for(int j=0;j<freq_ratio_vec.size();++j)
						{	
							if(freq_ratio_vec[j]!=-1&&j!=i)
							{							
								loglikelihood+=log(1-freq_ratio_vec[j]);								
							}	
						}
						if(loglikelihood>-DBL_MAX){
							subsample_valid=true;
						}	
						double mean_loglikelihood=loglikelihood/valid_pop_num;
						loglikelihood_vec.push_back(mean_loglikelihood);
					}else{
						loglikelihood_vec.push_back(1);
					}
				}
				//generate likelihood output
				string likelihood_str=find_str_after_nth_char(hline_str,0,'\t')+'\t'+find_str_after_nth_char(hline_str,1,'\t');
				for (int i = 0; i != loglikelihood_vec.size(); ++i)
				{
					stringstream temp_stream;
					temp_stream.precision(4);
					temp_stream<<fixed<<loglikelihood_vec[i];
					string op_lklhd=temp_stream.str();
				 	likelihood_str+='\t'+op_lklhd;
				 } 
				 if(cumu_total>=param.eff_sample*param.sample_frac)
				 {
				 	temp_vec.push_back(likelihood_str);
				 }
			}
			//output freq
			hline_str+="|"+line_vec[3]+">"+line_vec[4];
			//output tar_pop
			if(max_freq>0){
				hline_str+="|"+max_pop;
			}
			

			if(param.combine_all)
			{
				op_str=hline_str;
				if(param.pop_num==1)
				{
					for (int i = 0; i < pop_vec.size(); ++i)
					{
						op_str+='\t'+freq_dict[pop_vec[i]];
						vector<string> freq_vec=read_char_delim_str(freq_dict[pop_vec[i]],'/');
						int nom=stoi(freq_vec[0]),denom=stoi(freq_vec[1]);
						if(denom==0)
						{
							continue;
						}
						double freq_1=(double)nom/(double)denom,freq_2=(double)(cumu_case-nom)/(double)(cumu_total-denom);
						if(freq_2>=param.pop2_lower&&freq_2<=param.pop2_upper&&freq_1>=param.pop1_lower&&freq_1<=param.pop1_upper)
						{
							ind1*=0;
							var_pop_vec.push_back(pop_vec[i]);
						}
					}
				}
			}else{
				if(pop_neither_num==0)
				{
					if(param.max_homo_pop>1)
					{
						if(pop1_num>=1&&pop1_num<=param.max_homo_pop)
						{
							op_str=hline_str;
							for (int i = 0; i < pop_vec.size(); ++i)
							{
								op_str+='\t'+freq_dict[pop_vec[i]];
							}
							//Store gt_data
							ind1=0;
						}
					}else if(param.pop_num==pop1_num)
					{
						op_str=hline_str;
						for (int i = 0; i < pop_vec.size(); ++i)
						{
							if(freq_dict.count(pop_vec[i]))
							{
								op_str+='\t'+freq_dict[pop_vec[i]];
							}else{
								op_str+="\tNA";
							}
						}
						//Store gt_data
						ind1=0;
					}				
				}
			}

			vector<string> temp_vec2;
			temp_vec2.push_back(op_str);
			temp_vec2.push_back(gt_str);
			//Store ind variable
			ind_vec.push_back(valid_ind1);
			ind_vec.push_back(valid_ind2);
			ind_vec.push_back(ind1);
			ind_vec.push_back(ind2);
			ind_vec.push_back(ind3);
			ind_vec.push_back(ind4);
			ind_vec.push_back(ind5);
			return make_pair(temp_vec2,temp_vec);
		}

		//Machine learning module
		void ml_variant_reader(int j,pop_data& pop,input_param& param,var_list& var,ann_data& ann,int thread_num)
		{
			input.open(input_ad);
			check_file_open_status(input,input_ad);
			//Load true positive labels(if provided)
			//Create index dict array for different subsampling experiments
			cout<<"Generate partition indices..."<<endl<<endl;
			//Load pop_vec here
			for (unordered_map<string,set<string> >::iterator i =  pop.pop_dict.begin(); i != pop.pop_dict.end(); ++i)
			{
				if(i->second.size()>=param.min_sample)
				{
					pop_vec.push_back(i->first);
				}
			}
			//Load subsample idx dict
			for (int i = 0; i != param.experiment_times; ++i)
			{
				param.seed=i;
				subsample_idx_dict_vec.push_back(sample_from_valid_population(pop,param));
			}
			cout<<"Indices generting complete..."<<endl<<endl;

			//Examine population size
			if(param.verbose)
			{
				for (int i = 0; i < subsample_idx_dict_vec.size(); ++i)
				{
					cout<<"Subsample: "<<i<<endl;
					for (unordered_map<string,set<int>>::iterator j = subsample_idx_dict_vec[i].begin();j != subsample_idx_dict_vec[i].end(); ++j)
					{
						if(j->second.size()>=param.min_sample)
						{
							cout<<'\t'<<j->first<<":"<<j->second.size();
						}
					}
					cout<<endl;
				}
			}

			//Construct llthreshold dict
			//Dim1 thres
			//Dim2 subsamples
			//Dim3 diags
			for (int i = 0; i < 500; ++i)
			{
				double thres=-(double)i/1000;
				for (int j = 0; j != param.experiment_times; ++j)
				{
					llthres_dict[thres].push_back({0,0,0,0});//TP,FP,TN,FN
				}
			}

			cout<<"Diagnostic dictionary constructed."<<endl<<endl;

			int line_cter=1,var_cter;
			string line,chr_str,pos_str,ann_str;
			main_analysis_module main;
			auto start=chrono::high_resolution_clock::now();
			int valid_subpopulation_num=0;
			var_num=0;

			while(getline(input,line))
			{
				if(line[0]!='#')
				{
					if(line_cter%100000==0)
					{
						cout<<line_cter<<" total variants scanned."<<endl;
					}

					line_cter++;
					var_cter=sp_q.size();

					while(sp_q.size()>=(thread_num-1)*200000)
					{
						this_thread::sleep_for(chrono::milliseconds(1000));
						cout<<sp_q.size()<<" in-memory variants currently waiting for processing..."<<endl;
						//cout<<"Processing speed "<<(var_cter-sp_q.size())*60<<" variants/min."<<endl;
						var_cter=sp_q.size();
					}

					pair<string,int> temp_pair=make_pair(line,line_cter);
					unique_lock<mutex> ul(ip_mutex);
					sp_q.push(temp_pair);
					g_ready=true;
					ul.unlock();
					g_cv.notify_all();	
					
				}
			}
			input.close();

			cout<<"All "<<line_cter-1<<" variants scanned. Learning likelihood threshold for classification..."<<endl;
			var_cter=sp_q.size();
			cout<<var_cter<<" variants remained for processing..."<<endl;
			while(sp_q.size()!=0)
			{
				this_thread::sleep_for(chrono::milliseconds(10000));
				cout<<sp_q.size()<<" variants remained for processing...Processing speed "<<(var_cter-sp_q.size())*6<<" variants/min."<<endl;
				var_cter=sp_q.size();
			}

			cout<<"Calculation completed. ";
			auto end=chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(end - start);
			cout<<"Time lapsed: "<<duration.count()<<" seconds."<<endl<<endl;

			completed=true;
			g_cv.notify_all();
			g_ready=true;

			//Output classification boundaries
			map<double,vector<int>> empty_diag_dict;
			vector<map<double,vector<int>>> sample_run_diag_dict(param.experiment_times,empty_diag_dict);
			//Sort diag result by sampling runs
			for (map<double,vector<vector<int>>>::iterator i = llthres_dict.begin(); i !=llthres_dict.end(); ++i)
			{
				for (int j = 0; j !=i->second.size() ; ++j)
				{
					sample_run_diag_dict[j][i->first].push_back(i->second[j][0]);
					sample_run_diag_dict[j][i->first].push_back(i->second[j][1]);
					sample_run_diag_dict[j][i->first].push_back(i->second[j][2]);
					sample_run_diag_dict[j][i->first].push_back(i->second[j][3]);
				}
			}
			// Output diagnostics by sampling runs
			vector<bool> valid_vec(param.experiment_times,true);

			for (int i = 0; i != sample_run_diag_dict.size(); ++i)
			{
				// cout<<"Threshold\tTP\tFP\tTN\tFN\tTPR\tFPR\tPrecision\tRecall\tAccuracy"<<endl;
				for (map<double,vector<int>>::reverse_iterator j = sample_run_diag_dict[i].rbegin(); j != sample_run_diag_dict[i].rend(); ++j)
				{
					int total=j->second[0]+j->second[1]+j->second[2]+j->second[3];
					double TPR=(double) j->second[0]/(j->second[0]+j->second[3]);
					double FPR=(double) j->second[1]/(j->second[2]+j->second[1]);
					double acc=(double) (j->second[0]+j->second[2])/total;
					double precision=(double) j->second[0]/(j->second[0]+j->second[1]);
					double recall=(double) j->second[0]/(j->second[0]+j->second[3]);
					// cout<<j->first;
					// cout<<'\t'<<j->second[0];
					// cout<<'\t'<<j->second[1];
					// cout<<'\t'<<j->second[2];
					// cout<<'\t'<<j->second[3];
					// cout<<'\t'<<TPR;
					// cout<<'\t'<<FPR;
					// cout<<'\t'<<precision;
					// cout<<'\t'<<recall;
					// cout<<'\t'<<acc;
					// cout<<endl;
					if(isnan(TPR)||isnan(FPR))
					{
						valid_vec[i]=false;
					}
				}
			}

			//Exclude ineffective sampling runs
			int valid_run_cnt=0;
			double cumu_mean_llh=0;

			for (int i = 0; i != valid_vec.size(); ++i)
			{
				if(valid_vec[i])
				{
					valid_run_cnt++;
					double optimal_diff=-1;
					double optimal_TPR,optimal_FPR,optimal_thres;
					for (map<double,vector<int>>::reverse_iterator j = sample_run_diag_dict[i].rbegin(); j != sample_run_diag_dict[i].rend(); ++j)
					{
						double TPR=(double) j->second[0]/(j->second[0]+j->second[3]);
						double FPR=(double) j->second[1]/(j->second[2]+j->second[1]);
						double diff=TPR-FPR;
						if(diff>optimal_diff)
						{
							optimal_diff=diff;
							optimal_TPR=TPR;
							optimal_FPR=FPR;
							optimal_thres=j->first;
						}
					}
					cout<<"Subsample run "<<i+1<<':'<<endl;
					cout<<optimal_TPR<<'\t'<<optimal_FPR<<endl;
					cout<<optimal_thres<<'\t'<<optimal_diff<<endl;
					cumu_mean_llh+=optimal_thres;
				}
			}
			cout<<"Valid run count:"<<valid_run_cnt<<endl;
			double mean_optimal_thres=cumu_mean_llh/valid_run_cnt;
			cout<<"Mean llh threshold for classification is: "<<mean_optimal_thres<<endl<<endl<<endl;
			param.mean_llh=mean_optimal_thres;
			param.sample_frac=1;
		}

		void ml_loader(pop_data& pop,input_param& param,var_list& var,ann_data& ann)
		{
			vector<pair<double,string>> max_ll_vec; 
			string signature_pop;
			//Dim1=>subsamples
			//Dim2=>maxlr
			//Dim3=>pop name
			pair<string,int> q_pair;
			// cout<<"ML loader launched."<<endl;
			while(true){
				if(completed){
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0){
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=sp_q.front();
					sp_q.pop();
					ul.unlock();
					//Main analysis
					string ip_str=q_pair.first;
					int line_mark=q_pair.second;
					bool llisvalid=false;
					vector<string> line_vec=read_char_delim_str(ip_str,'\t');
					pair<vector<pair<double,string>>,string> temp_pair=classification_boundary_learner(pop,param,var,ann,line_vec,pop_vec,llisvalid);
					max_ll_vec=temp_pair.first;
					signature_pop=temp_pair.second;
					//Update diagnostics dict
					ul.lock();	
					//Traversing and update logllh_dict
					//Retrieve labels here
					if(llisvalid)
					{
						//Traverse llthres vec
						//{TP,FP,TN,FN}
						for (map<double,vector<vector<int>>>::iterator i = llthres_dict.begin(); i !=llthres_dict.end(); ++i)
						{
							double thres=i->first;
						
							//Traverse different subsample runs
							for (int j = 0; j != i->second.size(); ++j)
							{
								// if(j==0){
								// 	cout<<thres<<'\t'<<max_ll_vec[j].first<<'\t'<<max_ll_vec[j].second<<'\t'<<signature_pop<<endl;
								// }
								
								if(max_ll_vec[j].first==9){
									//Skip
								}else{
									if(max_ll_vec[j].first>thres)
									{
										//signal fired
										if(max_ll_vec[j].second==signature_pop)
										{
											//TP
											i->second[j][0]++;
										}else{
											//FP
											i->second[j][1]++;
										}
									}else{
										//signal not fired
										if(max_ll_vec[j].second=="NA"){
											//Skip
										}else{
											if(signature_pop!="NA")
											{
												//FN
												i->second[j][3]++;
											}else{
												//TN
												// if(j==1)
												// {
													
												// 	if(max_ll_vec[j].first<-10){
												// 		cout<<thres<<'\t'<<max_ll_vec[j].first<<'\t'<<max_ll_vec[j].second<<'\t'<<signature_pop<<endl;
												// 	}
												// }
												i->second[j][2]++;
											}
										}
									}
								}

								// if(j==0){
								// 	cout<<'\t'<<i->second[j][0];
								// 	cout<<'\t'<<i->second[j][1]<<'\t';
								// 	cout<<'\t'<<i->second[j][2]<<'\t';
								// 	cout<<'\t'<<i->second[j][3]<<endl;									
								// }

							}
						}
					}
					g_ready=false;
					ul.unlock();											
				}				
			}
			completed=true;
			g_cv.notify_all();
			g_ready=true;
		}

		pair<vector<pair<double,string>>,string> classification_boundary_learner(pop_data& pop,input_param& param,var_list& var,ann_data& ann,vector<string>& line_vec,vector<string>& pop_vec,bool & llisvalid)
		{
			vector<pair<double,string>> rst_mat;
			vector<double> freq_ratio_vec;
			vector<pair<double,string>> max_ll_vec;
			vector<string> pop_subsample_vec;
			string signature_pop="NA";
			int total_num=0,case_num=0;
			int cumu_total=0,cumu_case=0;
			int pop_tar_num=0,pop_neither_num=0;
			double freq;
			bool not_all_var=false;
			bool not_all_ref=false;
			//Label variants based on provided parameters
			// cout<<"Calculating var freq."<<endl;
			for (int i = 0; i < pop_vec.size(); ++i)
			{
				//Discover signatures if needed
				//Label true or false positive signatures based on hard clipper(if labels are not provided)
				total_num=0;
				case_num=0;
				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j)
				{
					if(param.depth_file.length()!=0)
					{

						if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j])))
						{
							total_num+=1;
							int temp;
							if(param.flip){
								temp=check_genotype_ref(line_vec[*j]);
							}else{
								temp=check_genotype(line_vec[*j]);
							}
							case_num+=temp;
						}
					}else{
						if(check_read_depth(line_vec[*j],param.min_depth))
						{
							total_num+=1;
							int temp;
							if(param.flip)
							{
								temp=check_genotype_ref(line_vec[*j]);
							}else{
								temp=check_genotype(line_vec[*j]);
							}
							case_num+=temp;
						}
					}
				}
				cumu_total+=total_num;
				cumu_case+=case_num;
				if(total_num>=param.min_sample)
				{
					freq=(double)case_num/(double)total_num;
					if(!param.combine_all)
					{
						if (freq>=param.pop1_lower&&freq<=param.pop1_upper)
						{
							signature_pop=pop_vec[i];
							pop_tar_num+=1;
						}else{
							if (!(freq>=param.pop2_lower&&freq<=param.pop2_upper))
							{
								pop_neither_num+=1;
							}
						}		
					}
				}
				//not all are variants
				if (total_num!=case_num)
				{
					not_all_var=true;
				}	
				if(total_num!=0)
				{
					not_all_ref=true;
				}
			}

			//Determine if signature pop qualifies
			if(pop_tar_num==param.pop_num&&pop_neither_num==0&&cumu_total>=param.eff_sample&&not_all_var&&not_all_ref)
			{
				//Valid signature
			}else{
				//Not a signature
				signature_pop="NA";
			}
			// cout<<signature_pop<<endl;

			//*******************************************************************//
			//Calculate log likelihood function for different subsamples
			for (int j = 0; j !=subsample_idx_dict_vec.size(); ++j)
			{
				int valid_pop_num=0;
				pop_subsample_vec.clear();
				freq_ratio_vec.clear();
				
				//Traverse all subsamples by population label
				for (int i = 0; i !=pop_vec.size(); ++i)
				{
					//Examine for valid subsampled population
					if(subsample_idx_dict_vec[j].count(pop_vec[i])&&subsample_idx_dict_vec[j][pop_vec[i]].size()>=param.min_sample)
					{
						
						total_num=0;
						case_num=0;
						for (set<int>::iterator k=subsample_idx_dict_vec[j][pop_vec[i]].begin(); k!=subsample_idx_dict_vec[j][pop_vec[i]].end(); ++k)
						{
							if(param.depth_file.length()!=0)
							{
								if(check_read_depth(line_vec[*k],min(10.0,pop.col_depth_dict[*k])))
								{
									total_num+=1;
									int temp;
									if(param.flip)
									{
										temp=check_genotype_ref(line_vec[*k]);
									}else{
										temp=check_genotype(line_vec[*k]);
									}
									case_num+=temp;
								}
							}else{
								if(check_read_depth(line_vec[*k],param.min_depth))
								{
									total_num+=1;
									int temp;
									if(param.flip)
									{
										temp=check_genotype_ref(line_vec[*k]);
									}else{
										temp=check_genotype(line_vec[*k]);
									}
									case_num+=temp;
								}
							}
						}
						if(total_num>=param.min_sample)
						{
							//Calcualte total effective sample number at the locus
							cumu_total+=total_num;
							cumu_case+=case_num;	
							//Calculate frequency for valid populations
							valid_pop_num++;
							freq=(double)case_num/(double)total_num;
							pop_subsample_vec.push_back(pop_vec[i]);
							freq_ratio_vec.push_back(freq);	
							// if(j==0){
							// 	cout<<'\t'<<freq;
							// }		
						}else{
							pop_subsample_vec.push_back("NA");
							freq_ratio_vec.push_back(-1);
						}		
					}
				}

				// if(j==0){
				// 	cout<<endl<<valid_pop_num<<endl;
				// }

				//Calculate likelihood ratio for all valid sample
				vector<double> loglikelihood_vec;
				double max_llh=-DBL_MAX;
				string pop_name="NA";
				for (int i = 0; i !=freq_ratio_vec.size(); ++i)
				{
					if(freq_ratio_vec[i]!=-1)//Traverse valid populations only
					{
						double loglikelihood=log(freq_ratio_vec[i]);
						for(int j=0;j<freq_ratio_vec.size();++j)
						{	
							if(freq_ratio_vec[j]!=-1&&j!=i)
							{											
								loglikelihood+=log(1-freq_ratio_vec[j]);
							}	
						}

						if(loglikelihood>-DBL_MAX)
						{
							llisvalid=true;
						}	
						double mean_loglikelihood=loglikelihood/valid_pop_num;

						if(mean_loglikelihood>max_llh)
						{
							max_llh=mean_loglikelihood;
							pop_name=pop_subsample_vec[i];
						}
						loglikelihood_vec.push_back(mean_loglikelihood);
					}else{
						loglikelihood_vec.push_back(1);
					}
				}



				//Output maxllh
				if(llisvalid)
				{
					if(cumu_total>=param.eff_sample*param.sample_frac)
					{
						// if(j==0){
						// cout<<endl<<max_llh<<":"<<pop_name<<":"<<llisvalid;
						// cout<<endl;							
						// }

						max_ll_vec.push_back({max_llh,pop_name});
					}else{
						max_ll_vec.push_back({9,"NA"});
					}					
				}else{
					max_ll_vec.push_back({9,"NA"});
				}
			}

			return {max_ll_vec,signature_pop};
		}

		//Learning subsamples with fixed parameters

		void llh_ml_loader(pop_data& pop,input_param& param,var_list& var,ann_data& ann)
		{
			pair<string,string> lr_freq_pair;
			pair<string,int> q_pair;
			// cout<<"LLH ml loader launched."<<endl;
			while(true)
			{
				if(completed)
				{
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0)
				{
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=sp_q.front();
					sp_q.pop();
					ul.unlock();
					//main analysis
					string ip_str=q_pair.first;
					int line_mark=q_pair.second;
					vector<string> line_vec=read_char_delim_str(ip_str,'\t');
					bool is_signature=false;
					lr_freq_pair=signature_discovery_ml(pop,param,var,ann,line_vec,pop_vec,is_signature);
					string lr_str=lr_freq_pair.first;
					string freq_str=lr_freq_pair.second;
					//output data
					ul.lock();	
					if(is_signature)
					{
						//Output frequency
						if(freq_str!="")
						{
							output<<line_mark<<'\t'<<freq_str<<endl;
							var_num++;
						}
						//Output likelihood
						output2<<lr_str<<endl;
					}
					g_ready=false;
					ul.unlock();									
				}				
			}
			completed=true;
			g_cv.notify_all();
			g_ready=true;
		}

		pair<string,string> signature_discovery_ml(pop_data& pop,input_param& param,var_list& var,ann_data& ann,vector<string>& line_vec,vector<string>& pop_vec,bool& is_signature)
		{
			int total_num,case_num,total_num2,case_num2,pop1_num=0,pop_neither_num=0;
			int valid_pop_num=0;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			vector<double> freq_ratio_vec;
			unordered_map<string,string> freq_dict;
			string op_str="";
			int cumu_total=0,cumu_case=0;
			double mean_ll_thres=param.mean_llh;
			//Analyze with full dataset
			for (int i = 0; i < pop_vec.size(); ++i)
			{
				total_num=0;
				case_num=0;
				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j)
				{	
					if(param.depth_file.length()!=0)
					{
						if(check_read_depth(line_vec[*j],min(10.0,pop.col_depth_dict[*j])))
						{
							total_num+=1;
							int temp;
							if(param.flip){
								temp=check_genotype_ref(line_vec[*j]);
							}else{
								temp=check_genotype(line_vec[*j]);
							}
							case_num+=temp;
						}
					}else{

						if(check_read_depth(line_vec[*j],param.min_depth))
						{
							total_num+=1;
							int temp;
							if(param.flip)
							{
								temp=check_genotype_ref(line_vec[*j]);
							}else{
								temp=check_genotype(line_vec[*j]);
							}
							case_num+=temp;
						}
					}
				}
				if(total_num>=param.min_sample)
				{
					cumu_total+=total_num;
					cumu_case+=case_num;							
				}
				if(total_num>=param.min_sample)
				{
					valid_pop_num++;
					freq=(double)case_num/(double)total_num;
					freq_ratio_vec.push_back(freq);
					freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
					if(!param.combine_all)
					{
						if (freq>=param.pop1_lower&&freq<=param.pop1_upper)
						{
							var_pop_vec.push_back(pop_vec[i]);
							pop1_num+=1;
						}else{
							if (!(freq>=param.pop2_lower&&freq<=param.pop2_upper))
							{
								pop_neither_num+=1;
							}
						}						
					}
				}else{
					freq_ratio_vec.push_back(-1);
					if (total_num==0)
					{
						freq_dict[pop_vec[i]]="0/0";
					}else{
						freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
					}
				}
			}
			//Calculate likelihood for all valid sample
			bool subsample_valid=false;
			vector<double> loglikelihood_vec;
			for (int i = 0; i !=freq_ratio_vec.size(); ++i)
			{
				if(freq_ratio_vec[i]!=-1)//Traverse valid populations only
				{
					double loglikelihood=log(freq_ratio_vec[i]);

					for(int j=0;j<freq_ratio_vec.size();++j)
					{
						if(freq_ratio_vec[j]!=-1&&j!=i)
						{							
							loglikelihood+=log(1-freq_ratio_vec[j]);								
						}	
					}
					if(loglikelihood>-DBL_MAX)
					{
						subsample_valid=true;
					}	
					double mean_loglikelihood=loglikelihood/valid_pop_num;
					loglikelihood_vec.push_back(mean_loglikelihood);
					// cout<<'\t'<<freq_ratio_vec[i]<<":"<<loglikelihood_vec[i];
				}else{
					loglikelihood_vec.push_back(1);
				}
			}
			// cout<<endl;
			//Find the max ll pop
			double max_ll=-DBL_MAX;
			int max_ll_idx=0;
			for (int i = 0; i!=loglikelihood_vec.size(); ++i)
			{
				if(loglikelihood_vec[i]!=1&&loglikelihood_vec[i]>max_ll){
					max_ll=loglikelihood_vec[i];
					max_ll_idx=i;
				}
			}
			//Compare max ll and the ll threshold
			double ll_thres=mean_ll_thres;
			string likelihood_str;
			// cout<<max_ll<<'\t'<<ll_thres<<'\t'<<valid_pop_num<<endl;
			if(max_ll>=ll_thres)
			{
				// cout<<max_ll<<'\t'<<ll_thres<<'\t'<<valid_pop_num<<endl;
				//Signature found
				is_signature=true;
				string hline_str=generate_line_header(line_vec);
				if(param.flip)
				{
					hline_str+="|ref";
				}
				//Generate likelihood output
				likelihood_str=find_str_after_nth_char(hline_str,0,'\t')+'\t'+find_str_after_nth_char(hline_str,1,'\t');
				for (int i = 0; i != loglikelihood_vec.size(); ++i)
				{
					stringstream temp_stream;
					temp_stream.precision(4);
					temp_stream<<fixed<<loglikelihood_vec[i];
					string op_lklhd=temp_stream.str();
				 	likelihood_str+='\t'+op_lklhd;
				 } 
				//Generate freq output
				hline_str+="|"+line_vec[3]+">"+line_vec[4];
				hline_str+="|"+pop_vec[max_ll_idx];
				op_str=hline_str;
				for (int i = 0; i < pop_vec.size(); ++i)
				{
					if(freq_dict.count(pop_vec[i]))
					{
						op_str+='\t'+freq_dict[pop_vec[i]];
					}else{
						op_str+="\tNA";
					}
				}

			}
			return make_pair(likelihood_str,op_str);
		}

		//Quantitative STR discovery module 

		void STR_data_loader(pop_data& pop,input_param& param){
			int line_mark;
			vector<string> line_vec;
			string ip_str,op_str;
			while(true)
			{
				if(completed)
				{
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0)
				{
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					pair<string,int> q_pair=sp_q.front();
					sp_q.pop();
					ul.unlock();
					//main analysis
					ip_str=q_pair.first;
					line_mark=q_pair.second;
					line_vec=read_char_delim_str(ip_str,'\t');
					//Read STR length
					int ind=0;
					bool is_str=true;
					if(param.analysis_mode==9)
					{
						op_str=analyze_STR_data(pop,param,line_vec,pop_vec,ind,is_str);	
					}
					//output data
					ul.lock();
					if(is_str)
					{
						eff_var_num++;
					}
					if(ind==0)
					{
						output<<line_mark<<'\t'<<op_str<<endl;
						var_num++;
					}
					g_ready=false;
					ul.unlock();
				}
			}
			completed=true;
			g_cv.notify_all();
			g_ready=true;						
		}

		string analyze_STR_data(pop_data& pop,input_param& param,vector<string>& line_vec,vector<string>& pop_vec,int& ind,bool& is_str){
			string op_str="",breed_STR;
			map<int,int> gt_repsz_dict;
			//Detecting if the locus is STR
			string ref_allele=line_vec[3];
			string alt_allele_field=line_vec[4];
			vector<string> alt_alleles=read_char_delim_str(alt_allele_field,',');
			alt_alleles.push_back(ref_allele);
			if(alt_alleles.size()<=2)
			{
				is_str=false;
				ind=1;
				return op_str;					
			}
			if(alt_alleles.size()>2)
			{
				if(!isSTR(alt_alleles))
				{
					is_str=false;
					ind=1;
					return op_str;					
				}
			}
			//Alter ref and alt allele representation
			line_vec[3]=alt_alleles.back();
			stringstream ss;
			ss<<alt_alleles[0];
			for (int i = 1; i !=alt_alleles.size()-1; ++i)
			{
				ss<<','<<alt_alleles[i];
			}
			line_vec[4]=ss.str();
			//Filter by repetitive element length
			int unit_length=find_str_after_nth_char(line_vec[3],0,')').length()-find_str_after_nth_char(line_vec[3],0,'(').length()-1;
			if(unit_length<param.min_rep_size)
			{
				ind=1;
				is_str=false;
				return op_str;
			}
			//Special escape
			if(is_str&&param.StrPrint)
			{
				op_str=line_vec[0]+'\t'+line_vec[1]+'\t'+line_vec[2]+'\t'+line_vec[3]+'\t'+line_vec[4]+'\t'+line_vec[6];
				return op_str;
			}
			// Read rep size
			string repsz_str;
			int repsz;
			repsz_str=find_str_after_nth_char(line_vec[3],1,')');
			repsz=stoi(repsz_str);
			gt_repsz_dict[0]=repsz;
			alt_alleles=read_char_delim_str(line_vec[4],',');
			for (int i = 0; i < alt_alleles.size(); ++i)
			{
				repsz_str=find_str_after_nth_char(alt_alleles[i],1,')');
				repsz=stoi(repsz_str);
				gt_repsz_dict[i+1]=repsz;
			}
			///Determin classification boundaries
			int rep_bound=param.rep_num;
			if(param.kmeans)
			{
				vector<double> rep_vec;
				for (map<int,int>::iterator i = gt_repsz_dict.begin(); i != gt_repsz_dict.end(); ++i)
				{
					rep_vec.push_back(i->second);
				}
				sort(rep_vec.begin(),rep_vec.end());
				rep_bound=kmeans_find_boundary(rep_vec,20);
			}
			//start analysis
			int pop_target_num=0;
			int pop_neither_num=0;
			for (int i = 0; i < pop_vec.size(); ++i)
			{
				breed_STR="";
				int total_num=0;
				int case_num=0;
				vector<double> repeat_vec;
				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j)
				{
					if(line_vec[*j][0]!='.')
					{
						//check_depth
						int allele_1=line_vec[*j][0]-'0',allele_2=line_vec[*j][2]-'0';
						string depth_str=read_char_delim_str(line_vec[*j],':')[1];
						vector<string> depth_vec=read_char_delim_str(depth_str,',');
						int total_depth;
						if(allele_1==allele_2)
						{
							total_depth=stoi(depth_vec[allele_1]);
						}else{
							total_depth=stoi(depth_vec[allele_1])+stoi(depth_vec[allele_2]);
						}
						if(total_depth>=param.min_depth)
						{
							// Calculate repeat number
							// repeat_vec.push_back((double)gt_repsz_dict[allele_1]);
							// repeat_vec.push_back((double)gt_repsz_dict[allele_2]);
							total_num++;
							int rep_num1=gt_repsz_dict[allele_1];
							int rep_num2=gt_repsz_dict[allele_2];
							if(param.isEXP)
							{
								if(rep_num1>=rep_bound&&rep_num2>=rep_bound)
								{
									case_num++;
								}								
							}else{
								if(rep_num1<rep_bound&&rep_num2<rep_bound)
								{
									case_num++;
								}
							}
						}
					}
				}
				double freq=(double) case_num/total_num;
				breed_STR=to_string(case_num)+'/'+to_string(total_num);
				if(total_num>=param.min_sample)
				{
					// Output mean repeat
					// stringstream temp_stream;
					// double mean=calc_mean(repeat_vec);
					// double se=calc_deviation(repeat_vec);
					// temp_stream<<fixed<<setprecision(4)<<mean<<"+-"<<se<<'('<<total_num<<')';
					// breed_STR=temp_stream.str();

					if(freq>=param.pop1_lower)
					{
						pop_target_num++;
					}else if(freq<param.pop1_lower&&freq>param.pop2_upper){
						pop_neither_num++;
					}
				}
				op_str+='\t'+breed_STR;
			}

			if(!param.isGS)
			{
				if(pop_neither_num!=0)
				{
					ind=1;
					return op_str;
				}
			}

			if(pop_target_num==0)
			{
				ind=1;
				return op_str;
			}

			if(param.max_homo_pop>1)
			{
				if(pop_target_num>param.max_homo_pop)
				{
					ind=1;
					return op_str;
				}
			}else{
				if(pop_target_num!=param.pop_num)
				{
					ind=1;
					return op_str;
				}
			}
			//Generate header line
			string hline_str=generate_line_header(line_vec);
			//Output allele type by rep number
			alt_alleles=read_char_delim_str(line_vec[4],',');
			alt_alleles.push_back(line_vec[3]);
			map<int,string> low_rep_dict,high_rep_dict;
			string low_rep_alleles="",high_rep_alleles="";
			for (int i = 0; i < alt_alleles.size(); ++i)
			{
				int cp_num=stoi(find_str_after_nth_char(alt_alleles[i],1,')'));
				if (cp_num>=rep_bound)
				{
					high_rep_dict[cp_num]=alt_alleles[i];
				}else{
					low_rep_dict[cp_num]=alt_alleles[i];
				}
			}
			for (map<int,string>::iterator i = low_rep_dict.begin(); i != low_rep_dict.end(); ++i)
			{
				if(low_rep_alleles=="")
				{
					low_rep_alleles+=i->second;
				}else{
					low_rep_alleles+=','+i->second;
				}				
			}
			for (map<int,string>::iterator i = high_rep_dict.begin(); i != high_rep_dict.end(); ++i)
			{
				if(high_rep_alleles=="")
				{
					high_rep_alleles+=i->second;
				}else{
					high_rep_alleles+=','+i->second;
				}				
			}
			hline_str=hline_str+"|"+low_rep_alleles+"|"+high_rep_alleles;
			op_str=hline_str+op_str;
			return op_str;
		}
		//Variants aggregation analysis module

		void gene_data_reader(pop_data pop,input_param param,var_list& var,ann_data ann){
			input.open(input_ad);
			check_file_open_status(input,input_ad);
			output.open(temp_ad);
			check_file_open_status(output,output_ad);
			int line_cter=1,var_cter;
			string line,chr_str,pos_str,ann_str;
			main_analysis_module main;
			auto start=chrono::high_resolution_clock::now();
			main.output_file_header(output,var,param,pop,pop_vec);			
			cout<<endl<<"Output header generated, start analyzing variant data."<<endl;
			var_num=0;
			while(getline(input,line))
			{
				if(line[0]!='#')
				{
					if(line_cter%100000==0)
					{
						cout<<line_cter<<" variants processed."<<endl;
					}
					line_cter++;
					if(param.ann_flag.length()>0)
					{
						ann_str=find_str_after_nth_char(line,7,'\t');
						if(find_str_after_nth_char(ann_str,1,'|').find(param.ann_flag)!=string::npos)
						{
							pair<string,int> temp_pair=make_pair(line,line_cter);
							unique_lock<mutex> ul(ip_mutex);
							sp_q.push(temp_pair);
							g_ready=true;
							ul.unlock();
							g_cv.notify_all();										
						}
					}else{
						pair<string,int > temp_pair=make_pair(line,line_cter);
						unique_lock<mutex> ul(ip_mutex);
						sp_q.push(temp_pair);
						g_ready=true;
						ul.unlock();
						g_cv.notify_all();	
					}
				}
			}
			input.close();
			cout<<"All "<<line_cter-1<<" variants scanned. Waiting for frequency calculation to be completed..."<<endl;
			var_cter=sp_q.size();
			cout<<var_cter<<" variants remained for processing..."<<endl;
			while(sp_q.size()!=0)
			{
				this_thread::sleep_for(chrono::milliseconds(10000));
				cout<<sp_q.size()<<" variants remained for processing...Processing speed "<<(var_cter-sp_q.size())*6<<" variants/min."<<endl;
				var_cter=sp_q.size();
			}
			cout<<"Frequency calculation completed. ";
			auto end=chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(end - start);
			cout<<"Time lapsed: "<<duration.count()<<" seconds."<<endl<<endl;
			completed=true;
			g_cv.notify_all();
			g_ready=true;
			//output per_gene_data
			for(unordered_map<string,unordered_map<int,int> >::iterator i=gene_mutation_dict.begin();i!=gene_mutation_dict.end();++i){
				output<<i->first<<endl;
				for (int j = 0; j < pop_vec.size(); ++j)
				{
					
				}
			}
			output.close();
			cout<<"A total of "<<var_num<<"  variant output."<<endl;
		}

		void gene_data_loader(pop_data pop,input_param param,var_list& var,ann_data ann){
			string ip_str,gene;
			int line_mark;
			vector<string> line_vec,pop_rst_vec;
			pair<string,int > q_pair;
			pair<string,unordered_map<string,int> > mutation_pair;
			unordered_map<int,int> mutation_dict;
			while(true)
			{
				if(completed)
				{
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0)
				{
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=sp_q.front();
					sp_q.pop();
					ul.unlock();
					//main analysis
					ip_str=q_pair.first;
					line_mark=q_pair.second;
					line_vec=read_char_delim_str(ip_str,'\t');
					vector<int> ind_vec;
					if(param.analysis_mode==7)
					{
						if(!param.no_splicing)
						{
							mutation_dict=count_variant_num_by_gene(pop,param,var,ann,line_vec,pop_vec);
						}else{
							if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1)
							{
								mutation_dict=count_variant_num_by_gene(pop,param,var,ann,line_vec,pop_vec);
							}
						}						
					}
					gene=find_str_after_nth_char(line_vec[7],3,'|');
					//output data
					//Need check
					ul.lock();
					for (int i = 8; i < line_vec.size(); ++i)
					{
						gene_mutation_dict[gene][i]+=mutation_dict[i];
					}
					var_num++;
					g_ready=false;
					ul.unlock();
				}						
			}
		}

		unordered_map<int,int> count_variant_num_by_gene(pop_data pop,input_param param,var_list& var,ann_data ann,vector<string> line_vec,vector<string> pop_vec){
			unordered_map<int,int> mutation_dict;
			for (int i = 8; i < line_vec.size(); ++i)
			{
				mutation_dict[i]=check_genotype(line_vec[i]);
			}
			return mutation_dict;
		}

		//Analysis launcher
		
		void multi_thread_freq_analysis(int t){
			if(param.analysis_mode==3||param.analysis_mode==1){
				thread t1(&thread_analysis_module::ed_variant_reader,this,0,ref(pop),ref(param),ref(var),ref(ann),t);
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::ed_freq_data_loader,this,i,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}	
			}else if(param.analysis_mode==2||param.analysis_mode==8){
				thread t1(&thread_analysis_module::ev_variant_reader,this,0,ref(pop),ref(param),ref(var),ref(ann));
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::ev_freq_loader,this,i,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}					
			}else if(param.analysis_mode==6||param.analysis_mode==4){
				thread t1(&thread_analysis_module::variant_data_reader,this,0,ref(pop),ref(param),ref(var),ref(ann),t);
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::freq_data_loader,this,i,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}				
			}else if(param.analysis_mode==5){
				thread t1(&thread_analysis_module::ed_variant_reader,this,0,ref(pop),ref(param),ref(var),ref(ann),t);
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::unipop_data_loader,this,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}	
			}else if(param.analysis_mode==7){
				thread t1(&thread_analysis_module::gene_data_reader,this,ref(pop),ref(param),ref(var),ref(ann));
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::gene_data_loader,this,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}	
			}else if(param.analysis_mode==9){
				thread t1(&thread_analysis_module::ed_variant_reader,this,0,ref(pop),ref(param),ref(var),ref(ann),t);				
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::STR_data_loader,this,ref(pop),ref(param)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}	
			}else if(param.analysis_mode==10){
				thread t1(&thread_analysis_module::ed_variant_reader,this,0,ref(pop),ref(param),ref(var),ref(ann),t);	
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::llh_ml_loader,this,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}	
			}else if(param.analysis_mode==11){
				//Learn classification boundary
				thread t1(&thread_analysis_module::ml_variant_reader,this,0,ref(pop),ref(param),ref(var),ref(ann),t);	
				for (int i = 1; i < t; ++i)
				{
					th_vec.push_back(thread(&thread_analysis_module::ml_loader,this,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i)
				{
					th_vec[i].join();
				}	
				completed=false;
				g_ready=false;

				//Discovery with learned boundary
				cout<<"Classification boundary learned, starting variant discovery..."<<endl;
				completed=false;
				g_ready=false;
				vector<thread> th_vec2;
				thread t2(&thread_analysis_module::ed_variant_reader,this,0,ref(pop),ref(param),ref(var),ref(ann),t);	
				for (int i = 1; i < t; ++i)
				{
					th_vec2.push_back(thread(&thread_analysis_module::llh_ml_loader,this,ref(pop),ref(param),ref(var),ref(ann)));
				}
				t2.join();
				for (int i = 0; i < th_vec2.size(); ++i)
				{
					th_vec2[i].join();
				}	
			}
		}


	private:
		ifstream input,temp_input;
		ofstream output,temp_output,output2;
		string input_ad,output_ad,temp_ad;
		bool completed=false;
		int var_num=0,row_cur=0,eff_var_num=0;
		unordered_map<string,unordered_map<int,int>> gene_mutation_dict;
		vector<int> subsample_idx_vec_disc,subsample_idx_vec_valid;
		unordered_map<string,set<int>> subsample_pop_col_dict;
		map<double,vector<vector<int>>> llthres_dict;//TP,FP,TN,FN
		vector<unordered_map<string,set<int>>> subsample_idx_dict_vec;
		unsigned int th_num=thread::hardware_concurrency();
		vector<int> cter_vec;
		vector<string> pop_vec;
		vector<string> subsample_pop_vec;
		map<double,vector<vector<int>>> logllh_dict;
		//dim1 logllh threshold
		//dim2 subsamples
		//dim3 diags
		queue<string> s_q;
		queue<pair<string,int>> sp_q;		
		queue<pair<string,pair<int,int>>> si_q;
		vector<thread> th_vec;
};


