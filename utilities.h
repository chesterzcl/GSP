#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<set>
#include<cmath>
#include<float.h>
#include<cstdlib>
#include<unordered_map>
#include<map>
#include<sstream>
#include<iomanip>
#include<unordered_set>
using namespace std;

void check_file_open_status(ifstream& input,string input_address){
	if(!input){
		cout<<"Failed to open the input file: "<< input_address<<endl<<endl;
		abort();
	}
	cout<<"File: "<<input_address<<" successfully opened."<<endl;
}

void check_file_open_status(ofstream& output,string output_address){
	if(!output){
		cout<<"Failed to open the output file: "<< output_address<<endl<<endl;
		abort();
	}
	cout<<"File: "<<output_address<<" successfully opened."<<endl;
}

vector<string> read_char_delim_str(string& line,char delim){
	vector<string> output;
	string temp="";
	for (int i = 0; i < line.size(); ++i){
		if (line[i]==delim){
			output.push_back(temp);
			temp="";
		}else{
			temp+=line[i];
		}
	}
	output.push_back(temp);
	return output;
}

string generate_line_header(vector<string> line_vec){
	vector<string> ann_vec=read_char_delim_str(line_vec[7],'|');
	string lheader_str;
	if(ann_vec.size()>1){
		lheader_str=line_vec[0]+'\t'+line_vec[1]+'\t'+ann_vec[3]+'\t'+ann_vec[1]+'\t'+ann_vec[13]+'\t'+ann_vec[10];
	}else{
		lheader_str=line_vec[0]+'\t'+line_vec[1];
	}
	return lheader_str;
}

string find_str_after_nth_char(string& str,int n,char delim){
	int cter=0,idx;
	string op_str;
	if(n>0){
		for (int i = 0; i < str.length(); ++i){
			if(str[i]==delim){
				cter++;
			}
			if(cter==n){
				idx=i+1;
				while(str[idx]!=delim&&idx<str.length()){
					op_str+=str[idx];
					idx++;
				}
				break;
			}
		}		
	}else{
		for (int i = 0; i < str.length(); ++i){
			if(str[i]==delim){
				break;
			}
			op_str+=str[i];
		}
	}
	return op_str;
}

int check_genotype(string& geno_field){
	if(geno_field[0]=='1'){
		return 1;
	}else{
		return 0;
	}
}

int check_genotype_ref(string& geno_field){
	if(geno_field[2]=='0'){
		return 1;
	}else{
		return 0;
	}
}

bool check_read_depth(string& geno_field,double min_depth){
	string depth_str=find_str_after_nth_char(geno_field,1,':');
	if(depth_str=="."){
		return false;
	}
	int rdp=stoi(find_str_after_nth_char(depth_str,0,','));
	int adp=stoi(find_str_after_nth_char(depth_str,1,','));
	if(rdp+adp<min_depth){
		return false;
	}
	return true;
}

template <class T> pair<vector<T>,vector<T>> split_vec_into_two(double ratio,int seed,vector<T>& ip_vec){
	int N=ip_vec.size();
	int N1=N*ratio;
	int N2=N-N1;
	vector<T> v1,v2;
	unordered_set<T> st;
	srand(seed);
	while(st.size()<N1){
		int idx=rand()%N;
		if(!st.count(idx)){
			v1.push_back(ip_vec[idx]);
			st.insert(idx);
		}
	}	
	for(int i=0;i!=N;++i){
		if(!st.count(i)){
			v2.push_back(ip_vec[i]);
		}
	}
	return make_pair(v1,v2);
}

// map<string,double> calculate_likelihood(vector<string>& line_vec, vector<string>& sample_idx_vec, unordered_map<string,set<int>>& sample_pop_idx_dict){
// 	for (unordered_map<string,set<int>>::iterator i = sample_pop_idx_dict.begin(); i != sample_pop_idx_dict.end(); ++i){
// 		int ;
// 		int ;
// 		for (set<int>::iterator j= i->second.begin(); j != i->second.end(); ++j){
// 			int idx=*j;
// 			string gene_str=line_vec[idx];
			
// 		}
		
// 	}

// }

string Find_diff_between_two_str_head(string str1,string str2){
	string long_str,short_str,diff="";
	if(str1.length()>str2.length()){
		long_str=str1;
		short_str=str2;
	}else if(str1.length()<str2.length()){
		long_str=str2;
		short_str=str1;
	}else{
		return diff;
	}
	//Try align from head
	for (int i = 0; i < short_str.length(); ++i){
		if(long_str[i]!=short_str[i]){
			return diff;
		}
	}
	diff=long_str.substr(short_str.length());
	return diff;
}

bool isRep(string& target,string& unit){
	int replen=unit.length();
	if(target.length()%replen==0){
		int rep_num=target.length()/replen;
		for (int i = 0; i <rep_num;++i){
			for (int j = 0; j <replen; ++j){
				if(target[i*replen+j]!=unit[j]){
					return false;
				}				
			}
		}
	}else{
		return false;
	}
	return true;
}

bool isSTR(vector<string>& input_vec){
	//Sort alleles by length
	map<int,string> dict;
	for (int i = 0; i < input_vec.size(); ++i){
		if(input_vec[i]=="*"){
			return false;
		}
		if(dict.count(input_vec[i].size())){
			return false;
		}
		dict[input_vec[i].length()]=input_vec[i];
	}
	vector<string> ordered_vec;
	for (map<int,string>::iterator i = dict.begin(); i != dict.end(); ++i){
		ordered_vec.push_back(i->second);
	}
	//Find differences between two strings
	vector<string> diff_vec;
	for (int i = 0; i < ordered_vec.size()-1; ++i){
		string diff=Find_diff_between_two_str_head(ordered_vec[i],ordered_vec[i+1]);
		if(diff.length()==0){
			return false;
		}
		diff_vec.push_back(diff);
	}
	string rep_unit="";
	bool found=false;
	vector<int> rep_num_vec;
	//Check if any is the rep of something(starting from the head,up to 6 nt)
	int first_element_len=diff_vec[0].length();
	for (int replen = 1; replen <= min(first_element_len,6); ++replen){
		string cand_rep=diff_vec[0].substr(0,replen);
		if(isRep(diff_vec[0],cand_rep)){
			rep_num_vec.clear();
			rep_num_vec.push_back(0);
			rep_num_vec.push_back(diff_vec[0].length()/cand_rep.length());
			//Check if all other diff is the rep of that
			for (int i = 1; i != diff_vec.size(); ++i){
				if(isRep(diff_vec[i],cand_rep)){
					rep_num_vec.push_back(diff_vec[i].length()/cand_rep.length());
				}else{
					goto labelC;
				}
			}
			found=true;
			rep_unit=cand_rep;
			break;
		}
		labelC:;
	}

	if(!found){
		return false;
	}
	//Modify input vec into str representation(also detecting the rep unit in the common sequence)
	int cumu_rep_num=0;
	for (int i = 0; i < ordered_vec.size(); ++i){
		int idx;
		for (int j = 0; j != input_vec.size(); ++j){
			if(ordered_vec[i]==input_vec[j]){
				idx=j;
			}
		}
		cumu_rep_num+=rep_num_vec[i];
		string STR=ordered_vec[0]+"("+rep_unit+")"+to_string(cumu_rep_num);
		input_vec[idx]=STR;
	}
	return true;
}

//t-test module
double calc_mean(vector<double>& arr){
   double sum = 0;
   for (int i = 0; i < arr.size(); i++)
      sum = sum + arr[i];
   return (double) sum / arr.size();
}
//calculating standard deviation
double calc_deviation(vector<double>& arr){
   double sum = 0;
   for (int i = 0; i < arr.size(); i++)
      sum = sum + (arr[i] - calc_mean(arr)) * (arr[i] - calc_mean(arr));
   return sqrt((double) sum / (arr.size() - 1));
}
//finding t-test statistics of two data
double calc_ttest(vector<double>& arr1,vector<double>& arr2){
   double mean1 = calc_mean(arr1);
   double mean2 = calc_mean(arr2);
   double sd1 = calc_deviation(arr1);
   double sd2 = calc_deviation(arr2);
   double t_test = (mean1 - mean2) / sqrt((sd1 * sd1) / arr1.size() + (sd2 * sd2) / arr2.size());
   return t_test;
}

int kmeans_find_boundary(vector<double>& input_vec,int iter){
	double center1=input_vec.front(),center2=input_vec.back();
	set<int> cluster1,cluster2;
	int iter_cter=0;
	while(iter_cter<=iter){
		//Asign cluster
		cluster1.clear();
		cluster2.clear();
		for (int i = 0; i < input_vec.size(); ++i){
			if(abs(input_vec[i]-center1)>=abs(input_vec[i]-center2)){
				cluster2.insert(input_vec[i]);
			}else{
				cluster1.insert(input_vec[i]);
			}
		}
		//Calculate cluster mean
		center1=0;
		for (set<int>::iterator i =cluster1.begin(); i != cluster1.end(); ++i){
			center1+=*i;
		}
		center1=(double) center1/cluster1.size();
		center2=0;
		for (set<int>::iterator i =cluster2.begin(); i != cluster2.end(); ++i){
			center2+=*i;
		}
		center2=(double) center2/cluster2.size();	
		iter_cter++;
	}
	return *(cluster2.begin());
}

class data_printer{
	public:
		template <class T> void print(vector<T> a){
		cout<<"Elements in the container:";
			for (int i = 0; i < a.size(); ++i){
				cout<<'\t'<<a[i];
			}
			cout<<endl;
		}

		void print(vector<vector<string> > a){
			for (int i = 0; i < a.size(); ++i){
				print(a[i]);
			}
		}

		void print(set<int> a){
			for (set<int>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<*i;
			}
			cout<<endl;
		}
		void print(set<string> a){
			for (set<string>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<*i;
			}
			cout<<endl;
		}
		void print(unordered_map<string, string> a){
			for (unordered_map<string, string>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<i->first<<'\t'<<i->second;
			}
			cout<<endl;
		}
		void print(unordered_map<string, int> a){
			for (unordered_map<string, int>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<i->first<<'\t'<<i->second;
			}
			cout<<endl;
		}
		void print(unordered_map<string, set<string> > a){
			for (unordered_map<string, set<string> >::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<":";
				print(i->second);
			}
		}
		void print(unordered_map<string, set<int> > a){
			for (unordered_map<string, set<int> >::iterator i = a.begin(); i != a.end(); ++i){
				cout<<i->first<<":";
				print(i->second);
			}
		}
		void print(unordered_map<string, double> a){
			for (unordered_map<string, double>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<i->first<<'\t'<<i->second;
			}
			cout<<endl;
		}
		void print(unordered_map<int, int> a){
			for (unordered_map<int, int>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<i->first<<'\t'<<i->second;
			}
			cout<<endl;
		}
		void print(unordered_map<int, string> a){
			for (unordered_map<int, string>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<i->first<<'\t'<<i->second;
			}
			cout<<endl;
		}
		void print(unordered_map<int, double> a){
			for (unordered_map<int, double>::iterator i = a.begin(); i != a.end(); ++i){
				cout<<'\t'<<i->first<<'\t'<<i->second;
			}
			cout<<endl;
		}
};



