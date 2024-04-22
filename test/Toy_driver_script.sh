#!/bin/bash

tool_dir=$(dirname $(pwd))
dir=$(pwd)


eff_sample=300
upper=0.9
lower=0.1

${tool_dir}/GSEXP \
unique \
-i ${dir}/Toy.vcf \
-o ${dir}/Toy_fixed_ref.txt \
-p ${dir}/Toy_pop_label.txt \
-t 2 \
--flip \
--effsample ${eff_sample} \
--mind 10 \
--mins 3 \
--popnum 1 \
--p1l ${upper} --p1u 1.0 --p2l 0.0 --p2u ${lower} 

${tool_dir}/GSEXP \
lhdisc \
-i ${dir}/Toy.vcf \
-o ${dir}/Toy_ml_ref.txt \
-p ${dir}/Toy_pop_label.txt \
-t 2 \
--flip \
--lh-thres -0.0031 \
--effsample ${eff_sample} \
--mind 10 \
--mins 3 \
--popnum 1 \
--p1l ${upper} --p1u 1.0 --p2l 0.0 --p2u ${lower} 

${tool_dir}/Postprocessing ${dir} Toy_fixed_ref
${tool_dir}/Postprocessing ${dir} Toy_ml_ref

echo "====================================================="
echo "Number of signature discovered with the adaptive discovery module:"
echo $(($(wc -l<${dir}/Toy_ml_ref_filtered_300.txt)-1))
echo "-----------------------------------------------------"
echo "Number of signature discovered with the fixed-frequency based discovery module:"
echo $(($(wc -l<${dir}/Toy_fixed_ref_filtered_300.txt)-1))
echo "====================================================="

