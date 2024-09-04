#include "utilities.h"
#include "pop_data.h"
#include "var_list.h"
#include "ann_data.h"
#include "input_param.h"
#include "main_analysis_module.h"
#include "thread_analysis_module.h"

using namespace std;

int main(int argc, char const *argv[]){
	cout<<"Engine initiated."<<endl;
	if(argc<2){
		cout<<"Insufficient number of arguments provided."<<endl;
		cout<<"Engine shut down."<<endl;
		return 0;
	}
	input_param param;
	param.reset_freq_param();
	param.read_parameters(argc,argv);
	string ann_base_input,var_pop,var_ip,var_op;
	thread_analysis_module tp;
	tp.param=param;
	if(tp.param.ann_file!=""){
		tp.ann.load_annotation_db(tp.param.ann_file);
	}
	cout<<"Launching analysis module."<<endl;
	tp.pop.load_population_data(tp.param.pop_file);
	tp.pop.index_data(tp.param.pop_file);
	if(tp.param.depth_file.length()!=0){
		tp.pop.load_sample_depth_sumstats(tp.param.depth_file);
	}
	// tp.param.print_input_parameters();
	tp.load_io_address(tp.param.vcf_file,tp.param.output_file);
	tp.param.launch_analysis_module();
	if(tp.param.analysis_mode==2||tp.param.analysis_mode==8){
		tp.param.load_population_label(tp.pop.pop_dict);
		tp.var.load_ev_variant_data(tp.param.var_list_file);
	}else if(tp.param.analysis_mode==4){
		if(tp.param.var_list_file!=""){
			tp.var.load_variant_data(tp.param.var_list_file);
		}
	}
	tp.multi_thread_freq_analysis(tp.param.thread_num);
	// if(tp.param.analysis_mode==11){
	// 	//Discovery with learned boundary
	// 	cout<<"Classification boundary learned, starting variant discovery..."<<endl;
	// 	tp.param.analysis_mode=10;
	// 	tp.multi_thread_freq_analysis(tp.param.thread_num);
	// }
	cout<<"Variant discovery completed."<<endl;
	return 0;
}




