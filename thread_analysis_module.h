#include<chrono>
#include<mutex>
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
			input.open(ip_str);
			int cur=0,cur_idx,temp_int;
			while(getline(input,line)){
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
			input.close();
			output.open(op_str);
			check_file_open_status(output,op_str);
			output<<header_str<<endl;
			for (int i = 0; i < sorted_vec.size(); ++i){
				output<<sorted_vec[i].second<<endl;
			}
			output.close();
		}	

		//Frequency analysis module

		void variant_data_reader(int j,pop_data pop,input_param param,var_list& var,ann_data ann){
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
			while(getline(input,line)){
				if(line[0]!='#'){
					if(line_cter%100000==0){
						cout<<line_cter<<" variants loaded."<<endl;
					}
					if(row_cur==var.var_mat.size()){
						break;
					}
					line_cter++;
					chr_str=find_str_after_nth_char(line,0,'\t');
					pos_str=find_str_after_nth_char(line,1,'\t');
					// cout<<chr_str<<'\t'<<pos_str<<endl;
					if(var.var_mat[row_cur][0]==chr_str){
						// cout<<"marker 1"<<endl;
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
							row_cur++;			
						}else{
							// cout<<"marker 2"<<endl;
							if(stoi(pos_str)>stoi(var.var_mat[row_cur][1])){
								row_cur++;
							}		
						}
						// cout<<"marker 3"<<endl;
					}else{
						if(stoi(chr_str.substr(7,2))>stoi(var.var_mat[row_cur][0].substr(7,2))){
							row_cur++;
						}
					}
				}
			}
			input.close();
			cout<<"All "<<row_cur<<" variants scanned. Waiting for frequency calculation to be completed..."<<endl;
			var_cter=si_q.size();
			cout<<var_cter<<" variants remained for processing..."<<endl;
			while(si_q.size()!=0){
				this_thread::sleep_for(chrono::milliseconds(10000));
				cout<<si_q.size()<<" variants remained for processing...Processing speed "<<(var_cter-si_q.size())*6<<" variants/min."<<endl;
				var_cter=si_q.size();
			}
			cout<<"Frequency calculation completed. ";
			auto end=chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(end - start);
			cout<<"Time lapsed: "<<duration.count()<<" seconds."<<endl<<endl;
			completed=true;
			g_cv.notify_all();
			g_ready=true;
			//Print out per-core processing lines
			// for (int i = 0; i < cter_vec.size(); ++i){
			// 	cout<<"Core "<<i<<" processed "<<cter_vec[i]<<" lines."<<endl;
			// }
			output.close();
			cout<<"A total of "<<var_num<<"  variant output."<<endl;
			cout<<"Sorting variant data..."<<endl;
			insert_sort_variant_files(temp_ad,output_ad);
			cout<<"Sorting completed."<<endl;
		}

		void freq_data_loader(int i,pop_data pop,input_param param,var_list& var,ann_data ann){
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
					cter_vec[i]++;
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
						ul.lock();
						output<<line_mark<<'\t'<<freq_str<<endl;
						var_num++;
						g_ready=false;
						ul.unlock();
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
				}
			}
		}

		pair<vector<string>,vector<string> > analyze_pop_freq_per_line(pop_data pop,input_param param,var_list& var,ann_data ann,vector<string> line_vec,vector<string> pop_vec,int row_mark,vector<int>& ind_vec){
			int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0;
			int total_num,case_num,total_num2,case_num2;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			unordered_map<string,string> freq_dict;
			string temp_str="",op_str="";
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
			op_str=var.var_mat[row_mark][0]+'\t'+var.var_mat[row_mark][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10];
			for (int i = 0; i < pop_vec.size(); ++i){
				op_str+='\t'+freq_dict[pop_vec[i]];
			}
			//Store gt_data
			temp_str=var.var_mat[row_mark][0]+'\t'+var.var_mat[row_mark][1]+'\t'+read_char_delim_str(line_vec[7],'|')[3]+'\t'+read_char_delim_str(line_vec[7],'|')[1]+'\t'+read_char_delim_str(line_vec[7],'|')[13]+'\t'+read_char_delim_str(line_vec[7],'|')[10]+temp_str;
			//Updating pop_vec			
			temp_vec.push_back(var.var_mat[row_mark][0]);
			temp_vec.push_back(var.var_mat[row_mark][1]);
			for (int i = 0; i < var_pop_vec.size(); ++i){
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

		void ed_variant_reader(int j,pop_data pop,input_param param,var_list& var,ann_data ann){
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
			while(getline(input,line)){
				if(line[0]!='#'){
					if(line_cter%100000==0){
						cout<<line_cter<<" variants processed."<<endl;
					}
					line_cter++;
					if(param.ann_flag.length()>0){
						ann_str=find_str_after_nth_char(line,7,'\t');
						if(find_str_after_nth_char(ann_str,1,'|').find(param.ann_flag)!=string::npos){
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
			while(sp_q.size()!=0){
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
			if(param.analysis_mode==3){
				var.sort_var_pop_mat();
			}
			cout<<"A total of "<<var_num<<"  variant output."<<endl;
			cout<<"Sorting variant data..."<<endl;
			insert_sort_variant_files(temp_ad,output_ad);
			filesystem::remove(temp_ad);
			cout<<"Sorting completed."<<endl;
		}

		void ed_freq_data_loader(int i,pop_data pop,input_param param,var_list& var,ann_data ann){
			string ip_str,freq_str,gt_str;
			int row_mark,line_mark;
			vector<string> line_vec,pop_rst_vec;
			pair<string,int > q_pair;
			pair<vector<string>,vector<string> > str_vec_pair;
			while(true){
				if(completed){
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0){
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
					if(!param.no_splicing){
						str_vec_pair=ed_analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
					}else{
						if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1){
							str_vec_pair=ed_analyze_pop_freq_per_line(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
						}
					}
					freq_str=str_vec_pair.first[0];
					gt_str=str_vec_pair.first[1];
					pop_rst_vec=str_vec_pair.second;
					pop_rst_vec.push_back(to_string(line_mark));
					//output data
					if(ind_vec[5]==0&&ind_vec[4]==0){
						if(param.analysis_mode==3){
							ul.lock();
							var.var_pop_mat.push_back(pop_rst_vec);
							output<<line_mark<<'\t'<<freq_str<<endl;
							var_num++;
							g_ready=false;
							ul.unlock();
						}else{
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

		pair<vector<string>,vector<string> > ed_analyze_pop_freq_per_line(pop_data pop,input_param param,var_list& var,ann_data ann,vector<string> line_vec,vector<string> pop_vec,int row_mark,vector<int>& ind_vec){
			int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0;
			int total_num,case_num,total_num2,case_num2;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			unordered_map<string,string> freq_dict;
			string gt_str="",op_str="";
			for (int i = 0; i < pop_vec.size(); ++i){
				total_num=0;
				case_num=0;
				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j){
					if(line_vec[*j][0]=='0'||line_vec[*j][0]=='1'){
						total_num+=1;
						case_num+=check_genotype(line_vec[*j]);
						gt_str+='\t'+to_string(check_genotype(line_vec[*j]));
					}else{
						if(line_vec[*j][0]=='.'){
							gt_str+="\t2";
						}else{
							gt_str+="\t3";
						}
					}
				}
				if(total_num>=param.min_sample){
					freq=(double)case_num/(double)total_num;
					freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
					if (freq>=param.pop1_lower&&freq<=param.pop1_upper){
						var_pop_vec.push_back(pop_vec[i]);
						ind4*=0;
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
			string hline_str=generate_line_header(line_vec);
			op_str=hline_str;
			for (int i = 0; i < pop_vec.size(); ++i){
				op_str+='\t'+freq_dict[pop_vec[i]];
			}
			//Store gt_data
			gt_str=hline_str+gt_str;
			//Updating pop_vec			
			temp_vec.push_back(line_vec[0]);
			temp_vec.push_back(line_vec[1]);
			for (int i = 0; i < var_pop_vec.size(); ++i){
				temp_vec.push_back(var_pop_vec[i]);
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

		//Shared pattern discovery module 

		void unipop_data_loader(pop_data pop,input_param param,var_list& var,ann_data ann){
			string ip_str,freq_str,gt_str;
			int row_mark,line_mark;
			vector<string> line_vec,pop_rst_vec;
			pair<string,int > q_pair;
			pair<vector<string>,vector<string> > str_vec_pair;
			while(true){
				if(completed){
					break;
				}
				unique_lock<mutex> ul(ip_mutex);
				if(sp_q.size()==0){
					g_cv.wait(ul,[](){return g_ready;});
				}else{
					q_pair=sp_q.front();
					// cter_vec[i]++;
					sp_q.pop();
					ul.unlock();
					//main analysis
					ip_str=q_pair.first;
					line_mark=q_pair.second;
					line_vec=read_char_delim_str(ip_str,'\t');
					vector<int> ind_vec;
					if(param.analysis_mode==5){
						if(!param.no_splicing){
							str_vec_pair=discover_unique_variant(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
						}else{
							if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1){
								str_vec_pair=discover_unique_variant(pop,param,var,ann,line_vec,pop_vec,row_mark,ind_vec);
							}
						}						
					}else{
						if(!param.no_splicing){

						}else{
							if (ann.gene_transcript_dict[read_char_delim_str(line_vec[7],'|')[3]].size()==1){

							}
						}		
					}

					freq_str=str_vec_pair.first[0];
					gt_str=str_vec_pair.first[1];
					pop_rst_vec=str_vec_pair.second;
					pop_rst_vec.push_back(to_string(line_mark));
					//output data
					if(ind_vec[2]==0&&ind_vec[4]==0){
						ul.lock();
						var.var_pop_mat.push_back(pop_rst_vec);
						output<<line_mark<<'\t'<<freq_str<<endl;
						var_num++;
						g_ready=false;
						ul.unlock();
					}
				}
			}
		}

		pair<vector<string>,vector<string> > discover_unique_variant(pop_data pop,input_param param,var_list& var,ann_data ann,vector<string> line_vec,vector<string> pop_vec,int row_mark,vector<int>& ind_vec){
			int ind1=1,ind2=1,ind3=1,ind4=1,ind5=1,valid_ind1=0,valid_ind2=0;
			int total_num,case_num,total_num2,case_num2,pop1_num=0,pop_neither_num=0;
			double freq;
			vector<string> var_pop_vec,temp_vec;
			unordered_map<string,string> freq_dict;
			string gt_str="",op_str="";
			for (int i = 0; i < pop_vec.size(); ++i){
				total_num=0;
				case_num=0;

				for (set<int>::iterator j=pop.pop_col_dict[pop_vec[i]].begin(); j!=pop.pop_col_dict[pop_vec[i]].end(); ++j){
					if(line_vec[*j][0]=='0'||line_vec[*j][0]=='1'){
						total_num+=1;
						case_num+=check_genotype(line_vec[*j]);
						gt_str+='\t'+to_string(check_genotype(line_vec[*j]));
					}else{
						if(line_vec[*j][0]=='.'){
							gt_str+="\t2";
						}else{
							gt_str+="\t3";
						}
					}
				}
				if(total_num>=param.min_sample){
					freq=(double)case_num/(double)total_num;
					freq_dict[pop_vec[i]]=to_string(case_num)+"/"+to_string(total_num);
					if (freq>=param.pop1_lower&&freq<=param.pop1_upper){
						var_pop_vec.push_back(pop_vec[i]);
						pop1_num+=1;
					}else{
						if (!(freq>=param.pop2_lower&&freq<=param.pop2_upper)){
							pop_neither_num+=1;
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

			if(pop1_num==param.pop_num&&pop_neither_num==0){
				string hline_str=generate_line_header(line_vec);
				op_str=hline_str;
				for (int i = 0; i < pop_vec.size(); ++i){
					op_str+='\t'+freq_dict[pop_vec[i]];
				}
				//Store gt_data
				gt_str=hline_str+gt_str;
				ind1=0;
			}
			//Updating pop_vec			
			temp_vec.push_back(line_vec[0]);
			temp_vec.push_back(line_vec[1]);
			for (int i = 0; i < var_pop_vec.size(); ++i){
				temp_vec.push_back(var_pop_vec[i]);
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

		//Analysis launcher
		
		void multi_thread_freq_analysis(int t){
			if(param.analysis_mode==3||param.analysis_mode==1){
				thread t1(&thread_analysis_module::ed_variant_reader,this,0,pop,param,ref(var),ann);
				for (int i = 1; i < t; ++i){
					th_vec.push_back(thread(&thread_analysis_module::ed_freq_data_loader,this,i,pop,param,ref(var),ann));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i){
					th_vec[i].join();
				}	
			}else if(param.analysis_mode==2){
				thread t1(&thread_analysis_module::variant_data_reader,this,0,pop,param,ref(var),ann);
				for (int i = 1; i < t; ++i){
					th_vec.push_back(thread(&thread_analysis_module::freq_data_loader,this,i,pop,param,ref(var),ann));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i){
					th_vec[i].join();
				}					
			}else if(param.analysis_mode==5){
				thread t1(&thread_analysis_module::ed_variant_reader,this,0,pop,param,ref(var),ann);
				for (int i = 1; i < t; ++i){
					th_vec.push_back(thread(&thread_analysis_module::unipop_data_loader,this,pop,param,ref(var),ann));
				}
				t1.join();
				for (int i = 0; i < th_vec.size(); ++i){
					th_vec[i].join();
				}	
			}
		}


	private:
		ifstream input,temp_input;
		ofstream output,temp_output;
		string input_ad,output_ad,temp_ad;
		bool completed=false;
		int var_num=0,row_cur=0;
		unsigned int th_num=thread::hardware_concurrency();
		vector<int> cter_vec;
		vector<string> pop_vec;
		queue<string> s_q;
		queue<pair<string,int> > sp_q;		
		queue<pair<string,pair<int,int> > > si_q;
		vector<thread> th_vec;
};


