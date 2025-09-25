#ifndef SEGMENT_INTEGRATED_ANALYSIS_H
#define SEGMENT_INTEGRATED_ANALYSIS_H

#include "main_analysis_module.h"
#include "segment_analysis_module.h"
#include "segment_output_module.h"
#include "utilities.h"
#include "pop_data.h"
#include "var_list.h"
#include "ann_data.h"
#include "input_param.h"

using namespace std;

// Extends main_analysis_module to add segment-first discovery capability
class segment_integrated_analysis : public main_analysis_module {
private:
    segment_analysis_module seg_analyzer;
    segment_output_module seg_output;
    vector<variant_locus> variant_loci;
    vector<genomic_segment> identified_segments;

public:
    // Constructor that initializes segment analyzer with parameters
    segment_integrated_analysis(const input_param& param) 
        : seg_analyzer(param.seg_bandwidth, param.seg_density_threshold, 
                      param.seg_merge_distance, param.seg_min_variants,
                      param.seg_likelihood_threshold, param.seg_adaptive_bandwidth) {}

    // Enhanced frequency analysis that includes likelihood scoring and regional density
    void segment_population_frequency_analysis(string input_address, string output_address, 
                                             pop_data pop, input_param param, 
                                             var_list& var, ann_data ann) {
        
        // Phase 1: Load VCF data and apply frequency filters first (sequential hierarchy)
        input.open(input_address);
        check_file_open_status(input, input_address);
        
        cout << "Phase 1: Loading variants and applying frequency filters..." << endl;
        
        variant_loci.clear();
        int row_cur = 0, row_cter = 1;
        vector<string> pop_vec;
        
        // Get population vector (reuse existing logic)
        for (unordered_map<string,set<string> >::iterator i = pop.pop_dict.begin(); 
             i != pop.pop_dict.end(); ++i) {
            if (i->second.size() >= param.min_sample) {
                pop_vec.push_back(i->first);
            }
        }
        
        // Process VCF file and collect variant data with sequential filtering
        while(getline(input, line)) {
            if(line[0] != '#') {
                if(row_cter % 10000 == 0) {
                    cout << row_cter << " variants processed." << endl;
                }
                row_cter++;
                
                if(row_cur == var.var_mat.size()) break;
                
                line_vec = read_char_delim_str(line, '\t');
                
                // Match VCF line with target variant
                if(var.var_mat[row_cur][0] == line_vec[0] && 
                   var.var_mat[row_cur][1] == line_vec[1]) {
                    
                    // SEQUENTIAL THRESHOLD HIERARCHY:
                    // Step 1: Apply frequency gates first (cheap filter)
                    double frequency = calculate_variant_frequency(pop, pop_vec, param, row_cur);
                    bool freq_pass = (frequency >= param.pop1_lower && frequency <= param.pop1_upper);
                    
                    if (freq_pass) {
                        // Step 2: Calculate likelihood only for frequency survivors
                        double likelihood_score = calculate_variant_likelihood(pop, pop_vec, param, row_cur);
                        
                        // Create variant_locus object for survivors
                        string variant_info = read_char_delim_str(line_vec[7], '|')[3] + "|" +
                                            read_char_delim_str(line_vec[7], '|')[1] + "|" +
                                            read_char_delim_str(line_vec[7], '|')[10];
                        
                        variant_locus variant(line_vec[0], stoi(line_vec[1]), 
                                            likelihood_score, frequency, variant_info);
                        variant_loci.push_back(variant);
                    }
                    
                    row_cur++;
                }
            }
        }
        input.close();
        
        cout << "Phase 1 completed. " << variant_loci.size() << " frequency-filtered variants loaded." << endl;
        
        // Phase 2: Calculate likelihood scores and apply φ(Lj) during KDE
        cout << "Phase 2: Calculating regional density scores with φ(Lj) gates..." << endl;
        seg_analyzer.calculate_regional_densities(variant_loci);
        cout << "Phase 2 completed. Regional density scores calculated." << endl;
        
        // Phase 3: Apply density threshold (final filter in hierarchy)
        cout << "Phase 3: Applying density threshold filter..." << endl;
        // Note: frequency filtering already applied, now just check density
        for (auto& variant : variant_loci) {
            variant.passes_criteria = (variant.regional_density >= seg_analyzer.regional_density_threshold);
        }
        
        int passing_variants = 0;
        for (const auto& variant : variant_loci) {
            if (variant.passes_criteria) passing_variants++;
        }
        cout << "Phase 3 completed. " << passing_variants << " variants pass sequential hierarchy." << endl;
        
        // Phase 4: Identify genomic segments
        cout << "Phase 4: Identifying genomic segments..." << endl;
        identified_segments = seg_analyzer.identify_segments(variant_loci);
        cout << "Phase 4 completed. " << identified_segments.size() << " segments identified." << endl;
        
        // Phase 5: Generate output files
        cout << "Phase 5: Writing output files..." << endl;
        seg_output.write_all_outputs(output_address, variant_loci, identified_segments, pop_vec);
        cout << "Phase 5 completed. All output files generated." << endl;
    }

    // Calculate likelihood score using integration priority system:
    // Priority: SigMl > SigLh > SigFreq pseudo-likelihood
    // Li = max_over_alleles_and_populations(log_likelihood)
    double calculate_variant_likelihood(pop_data& pop, vector<string>& pop_vec, 
                                      input_param& param, int row_cur) {
        
        double max_likelihood = 0.0;
        
        // LIKELIHOOD INTEGRATION PRIORITY SYSTEM:
        if (param.ml_mode && !param.likelihood_file.empty()) {
            // Priority 1: Use SigMl scores if available
            max_likelihood = get_ml_likelihood_score(row_cur, param);
            
        } else if (param.lh_mode) {
            // Priority 2: Use SigLh likelihood calculations if available
            max_likelihood = get_lh_likelihood_score(pop, pop_vec, param, row_cur);
            
        } else {
            // Priority 3: SigFreq pseudo-likelihood fallback
            max_likelihood = get_frequency_pseudo_likelihood(pop, pop_vec, param, row_cur);
        }
        
        return max_likelihood;
    }

private:
    // Priority 1: Get ML-predicted likelihood scores
    double get_ml_likelihood_score(int row_cur, input_param& param) {
        // TODO: Integrate with existing ML likelihood file loading
        // For now, placeholder that should read from param.likelihood_file
        
        // NOTE: Keep all scores on same scale (prefer log)
        // If linear ML outputs, apply logit transform or adjust threshold accordingly
        
        return 1.0; // Placeholder - implement file reading
    }
    
    // Priority 2: Get SigLh likelihood calculations  
    double get_lh_likelihood_score(pop_data& pop, vector<string>& pop_vec, 
                                  input_param& param, int row_cur) {
        // TODO: Integrate with existing SigLh algorithms
        // Should reuse existing likelihood calculation infrastructure
        
        double max_likelihood = 0.0;
        
        // Calculate likelihood across populations using proper log-likelihood
        for (const string& pop_name : pop_vec) {
            int total_samples = 0;
            int variant_samples = 0;
            
            for (set<int>::iterator j = pop.pop_col_dict[pop_name].begin(); 
                 j != pop.pop_col_dict[pop_name].end(); ++j) {
                
                if(line_vec[*j][0] == '0' || line_vec[*j][0] == '1') {
                    total_samples++;
                    variant_samples += check_genotype(line_vec[*j]);
                }
            }
            
            if (total_samples > 0) {
                double freq = (double)variant_samples / (double)total_samples;
                
                // TODO: Replace with proper log-likelihood calculation
                // This should use existing SigLh mathematical framework
                double pop_likelihood = calculate_log_likelihood(freq, total_samples, variant_samples);
                
                max_likelihood = max(max_likelihood, pop_likelihood);
            }
        }
        
        return max_likelihood;
    }
    
    // Priority 3: SigFreq pseudo-likelihood (fallback)
    double get_frequency_pseudo_likelihood(pop_data& pop, vector<string>& pop_vec, 
                                         input_param& param, int row_cur) {
        // IMPORTANT: Document scale! This should be consistent with Li scale
        
        double max_likelihood = 0.0;
        
        for (const string& pop_name : pop_vec) {
            int total_samples = 0;
            int variant_samples = 0;
            
            for (set<int>::iterator j = pop.pop_col_dict[pop_name].begin(); 
                 j != pop.pop_col_dict[pop_name].end(); ++j) {
                
                if(line_vec[*j][0] == '0' || line_vec[*j][0] == '1') {
                    total_samples++;
                    variant_samples += check_genotype(line_vec[*j]);
                }
            }
            
            if (total_samples > 0) {
                double freq = (double)variant_samples / (double)total_samples;
                
                // Pseudo-likelihood: frequency deviation from neutral
                // Scale documented: range [0, 1] for compatibility with other Li sources
                double pop_likelihood = abs(freq - 0.5) * 2.0;
                
                max_likelihood = max(max_likelihood, pop_likelihood);
            }
        }
        
        return max_likelihood;
    }
    
    // Helper: Calculate proper log-likelihood (placeholder for SigLh integration)
    double calculate_log_likelihood(double freq, int total_samples, int variant_samples) {
        // TODO: Implement proper binomial log-likelihood or integrate with existing SigLh
        // For now, return scaled frequency deviation as placeholder
        return abs(freq - 0.5) * 2.0;
    }
    
    // Calculate variant frequency (reuse existing logic)
    double calculate_variant_frequency(pop_data& pop, vector<string>& pop_vec, 
                                     input_param& param, int row_cur) {
        int total_samples = 0;
        int variant_samples = 0;
        
        for (const string& pop_name : pop_vec) {
            for (set<int>::iterator j = pop.pop_col_dict[pop_name].begin(); 
                 j != pop.pop_col_dict[pop_name].end(); ++j) {
                
                if(line_vec[*j][0] == '0' || line_vec[*j][0] == '1') {
                    total_samples++;
                    variant_samples += check_genotype(line_vec[*j]);
                }
            }
        }
        
        if (total_samples > 0) {
            return (double)variant_samples / (double)total_samples;
        }
        return 0.0;
    }
    
    // Get identified segments for external access
    const vector<genomic_segment>& get_segments() const {
        return identified_segments;
    }
    
    // Get variant loci data for external access
    const vector<variant_locus>& get_variant_loci() const {
        return variant_loci;
    }
    
public:
    // Print segment analysis parameters
    void print_segment_parameters(const input_param& param) {
        cout << "=== Segment-First Discovery Parameters ===" << endl;
        cout << "Gaussian kernel bandwidth: " << param.seg_bandwidth << " bp" << endl;
        cout << "Regional density threshold: " << param.seg_density_threshold << endl;
        cout << "Segment merge distance: " << param.seg_merge_distance << " bp" << endl;
        cout << "Minimum variants per segment: " << param.seg_min_variants << endl;
        cout << "Target frequency range: [" << param.pop1_lower << ", " << param.pop1_upper << "]" << endl;
        cout << "Minimum samples: " << param.min_sample << endl;
        cout << "=========================================" << endl << endl;
    }
    
    // Validate segment analysis parameters
    bool validate_segment_parameters(const input_param& param) {
        if (param.seg_bandwidth <= 0) {
            cout << "Error: Segment bandwidth must be positive." << endl;
            return false;
        }
        if (param.seg_density_threshold < 0) {
            cout << "Error: Regional density threshold cannot be negative." << endl;
            return false;
        }
        if (param.seg_merge_distance < 0) {
            cout << "Error: Segment merge distance cannot be negative." << endl;
            return false;
        }
        if (param.seg_min_variants < 1) {
            cout << "Error: Minimum variants per segment must be at least 1." << endl;
            return false;
        }
        return true;
    }
};

#endif // SEGMENT_INTEGRATED_ANALYSIS_H