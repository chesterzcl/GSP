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
		return 0;
	}
	input_param param;
	param.read_parameters(argc,argv);
	string ann_base_input,var_pop,var_ip,var_op;
	thread_analysis_module tp;
	tp.param=param;
	tp.ann.load_annotation_db(tp.param.ann_file);
	cout<<"Start variant discovery."<<endl;
	tp.var.load_variant_data(tp.param.vcf_file);
	tp.pop.load_population_data(tp.param.pop_file);
	tp.pop.index_data(tp.param.pop_file);
	tp.param.print_input_parameters();
	tp.load_io_address(tp.param.vcf_file,tp.param.output_file);
	tp.multi_thread_freq_analysis(tp.param.thread_num);
	cout<<"Variant discovery completed."<<endl;


	return 0;
}







