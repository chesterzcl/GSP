#ifndef SEGMENT_OUTPUT_MODULE_H
#define SEGMENT_OUTPUT_MODULE_H

#include <fstream>
#include <iomanip>
#include "segment_analysis_module.h"
#include "utilities.h"

using namespace std;

class segment_output_module {
public:
    // Write per-locus results with likelihood and density scores
    void output_locus_results(const string& output_file, 
                             const vector<variant_locus>& variants,
                             const vector<string>& population_names) {
        
        ofstream output(output_file);
        check_file_open_status(output, output_file);
        
        // Write header
        output << "#CHROMOSOME\tPOSITION\tLIKELIHOOD_SCORE\tFREQUENCY\t"
               << "REGIONAL_DENSITY\tPASSES_CRITERIA\tVARIANT_INFO";
        
        // Add population frequency columns
        for (const string& pop : population_names) {
            output << "\t" << pop;
        }
        output << endl;
        
        // Write variant data
        for (const auto& variant : variants) {
            output << variant.chromosome << "\t"
                   << variant.position << "\t"
                   << fixed << setprecision(6) << variant.likelihood_score << "\t"
                   << fixed << setprecision(4) << variant.frequency << "\t"
                   << fixed << setprecision(6) << variant.regional_density << "\t"
                   << (variant.passes_criteria ? "PASS" : "FAIL") << "\t"
                   << variant.variant_info;
            
            // Add population-specific frequency data if available
            // This would be populated from the main analysis pipeline
            output << endl;
        }
        
        output.close();
    }
    
    // Write segment results
    void output_segment_results(const string& output_file, 
                               const vector<genomic_segment>& segments) {
        
        ofstream output(output_file);
        check_file_open_status(output, output_file);
        
        // Write header
        output << "#CHROMOSOME\tSTART_POS\tEND_POS\tLENGTH_BP\tVARIANT_COUNT\t"
               << "MEAN_LIKELIHOOD\tMEAN_DENSITY\tVARIANT_POSITIONS" << endl;
        
        // Write segment data
        for (const auto& segment : segments) {
            int segment_length = segment.end_position - segment.start_position + 1;
            
            output << segment.chromosome << "\t"
                   << segment.start_position << "\t"
                   << segment.end_position << "\t"
                   << segment_length << "\t"
                   << segment.variant_count << "\t"
                   << fixed << setprecision(6) << segment.mean_likelihood << "\t"
                   << fixed << setprecision(6) << segment.mean_density << "\t";
            
            // List all variant positions in this segment
            for (size_t i = 0; i < segment.variants.size(); ++i) {
                if (i > 0) output << ",";
                output << segment.variants[i].position;
            }
            output << endl;
        }
        
        output.close();
    }
    
    // Write detailed segment analysis with per-variant info
    void output_detailed_segments(const string& output_file, 
                                 const vector<genomic_segment>& segments) {
        
        ofstream output(output_file);
        check_file_open_status(output, output_file);
        
        // Write header
        output << "#SEGMENT_ID\tCHROMOSOME\tSEGMENT_START\tSEGMENT_END\t"
               << "VARIANT_POSITION\tLIKELIHOOD_SCORE\tFREQUENCY\t"
               << "REGIONAL_DENSITY\tVARIANT_INFO" << endl;
        
        // Write detailed segment data
        for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
            const auto& segment = segments[seg_idx];
            string segment_id = "SEG_" + to_string(seg_idx + 1);
            
            for (const auto& variant : segment.variants) {
                output << segment_id << "\t"
                       << segment.chromosome << "\t"
                       << segment.start_position << "\t"
                       << segment.end_position << "\t"
                       << variant.position << "\t"
                       << fixed << setprecision(6) << variant.likelihood_score << "\t"
                       << fixed << setprecision(4) << variant.frequency << "\t"
                       << fixed << setprecision(6) << variant.regional_density << "\t"
                       << variant.variant_info << endl;
            }
        }
        
        output.close();
    }
    
    // Write summary statistics
    void output_segment_summary(const string& output_file, 
                               const vector<genomic_segment>& segments,
                               const vector<variant_locus>& all_variants) {
        
        ofstream output(output_file);
        check_file_open_status(output, output_file);
        
        // Calculate summary statistics
        int total_variants = all_variants.size();
        int passing_variants = 0;
        int variants_in_segments = 0;
        
        for (const auto& variant : all_variants) {
            if (variant.passes_criteria) passing_variants++;
        }
        
        for (const auto& segment : segments) {
            variants_in_segments += segment.variant_count;
        }
        
        // Write summary
        output << "# Segment-First Discovery Analysis Summary" << endl;
        output << "Total variants analyzed: " << total_variants << endl;
        output << "Variants passing dual thresholds: " << passing_variants << endl;
        output << "Variants organized into segments: " << variants_in_segments << endl;
        output << "Total segments identified: " << segments.size() << endl;
        output << "Percentage of passing variants in segments: " 
               << fixed << setprecision(2) 
               << (passing_variants > 0 ? (100.0 * variants_in_segments / passing_variants) : 0.0) 
               << "%" << endl;
        
        if (!segments.empty()) {
            // Calculate segment size statistics
            vector<int> segment_sizes;
            vector<int> segment_lengths;
            
            for (const auto& segment : segments) {
                segment_sizes.push_back(segment.variant_count);
                segment_lengths.push_back(segment.end_position - segment.start_position + 1);
            }
            
            sort(segment_sizes.begin(), segment_sizes.end());
            sort(segment_lengths.begin(), segment_lengths.end());
            
            output << endl << "# Segment Statistics" << endl;
            output << "Mean variants per segment: " 
                   << fixed << setprecision(1) 
                   << (accumulate(segment_sizes.begin(), segment_sizes.end(), 0.0) / segments.size()) << endl;
            output << "Median variants per segment: " 
                   << segment_sizes[segments.size() / 2] << endl;
            output << "Mean segment length (bp): " 
                   << fixed << setprecision(0)
                   << (accumulate(segment_lengths.begin(), segment_lengths.end(), 0.0) / segments.size()) << endl;
            output << "Median segment length (bp): " 
                   << segment_lengths[segments.size() / 2] << endl;
        }
        
        output.close();
    }
    
    // Main output function that creates all output files
    void write_all_outputs(const string& base_output_name,
                          const vector<variant_locus>& variants,
                          const vector<genomic_segment>& segments,
                          const vector<string>& population_names) {
        
        // Write per-locus results
        string locus_file = base_output_name + ".locus";
        output_locus_results(locus_file, variants, population_names);
        cout << "Per-locus results written to: " << locus_file << endl;
        
        // Write segment summary
        string segment_file = base_output_name + ".segments";
        output_segment_results(segment_file, segments);
        cout << "Segment results written to: " << segment_file << endl;
        
        // Write detailed segment analysis
        string detailed_file = base_output_name + ".segments.detailed";
        output_detailed_segments(detailed_file, segments);
        cout << "Detailed segment analysis written to: " << detailed_file << endl;
        
        // Write summary statistics
        string summary_file = base_output_name + ".summary";
        output_segment_summary(summary_file, segments, variants);
        cout << "Analysis summary written to: " << summary_file << endl;
    }
};

#endif // SEGMENT_OUTPUT_MODULE_H