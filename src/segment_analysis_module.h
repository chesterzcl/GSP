#ifndef SEGMENT_ANALYSIS_MODULE_H
#define SEGMENT_ANALYSIS_MODULE_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>
#include "utilities.h"
#include "pop_data.h"
#include "var_list.h"
#include "input_param.h"

using namespace std;

// Structure to hold variant data with likelihood and density scores
struct variant_locus {
    string chromosome;
    int position;
    double likelihood_score;    // Li - likelihood score for this locus
    double frequency;          // Variant frequency 
    double regional_density;   // Si - regional density score
    bool passes_criteria;      // Whether it passes dual thresholds
    string variant_info;       // Additional variant annotation
    
    variant_locus(string chr, int pos, double lh_score, double freq, string info = "") 
        : chromosome(chr), position(pos), likelihood_score(lh_score), 
          frequency(freq), regional_density(0.0), passes_criteria(false), 
          variant_info(info) {}
};

// Structure for genomic segments
struct genomic_segment {
    string chromosome;
    int start_position;
    int end_position;
    int variant_count;
    double mean_likelihood;
    double mean_density;
    vector<variant_locus> variants;
    
    genomic_segment(string chr, int start, int end) 
        : chromosome(chr), start_position(start), end_position(end), 
          variant_count(0), mean_likelihood(0.0), mean_density(0.0) {}
};

class segment_analysis_module {
public:
    // Configuration parameters
    double bandwidth_sigma;           // Gaussian kernel bandwidth (default: 10000 bp)
    double regional_density_threshold; // Minimum Si score for segment formation
    double segment_merge_distance;    // Distance to merge nearby segments
    double likelihood_threshold;      // Threshold T for φ(Lj) filter function
    bool adaptive_bandwidth;          // Whether to use adaptive bandwidth calculation
    int min_variants_per_segment;     // Minimum variants required for a segment
    
    // Constructor with default parameters
    segment_analysis_module(double sigma = 10000.0, double density_thresh = 0.1, 
                          double merge_dist = 50000.0, int min_vars = 3, 
                          double lh_thresh = 0.1, bool adaptive = false)
        : bandwidth_sigma(sigma), regional_density_threshold(density_thresh),
          segment_merge_distance(merge_dist), min_variants_per_segment(min_vars),
          likelihood_threshold(lh_thresh), adaptive_bandwidth(adaptive) {}

    // Calculate normalized Gaussian kernel density using the corrected formula:
    // Si = Σj exp(-dij²/2σ²) φ(Lj) Lj / Σj exp(-dij²/2σ²) φ(Lj)
    // φ(Lj) gate applied to neighbors only - central locus included regardless
    double calculate_regional_density(const variant_locus& target_variant, 
                                    const vector<variant_locus>& all_variants) {
        
        if (all_variants.empty()) return 0.0;
        
        double weighted_sum = 0.0;      // Numerator: Σj exp(-dij²/2σ²) φ(Lj) Lj
        double weight_sum = 0.0;        // Denominator: Σj exp(-dij²/2σ²) φ(Lj)
        
        // Calculate density contribution from each variant
        for (const auto& variant : all_variants) {
            // Skip if different chromosome - no cross-chromosome influence
            if (variant.chromosome != target_variant.chromosome) continue;
            
            // Calculate distance between variants (self-contribution: dii = 0)
            double distance = abs(variant.position - target_variant.position);
            
            // Calculate Gaussian kernel weight (note: 1/(σ√(2π)) cancels out)
            double exponent = -(distance * distance) / (2.0 * bandwidth_sigma * bandwidth_sigma);
            double gaussian_weight = exp(exponent);
            
            // Apply φ(Lj) threshold gate to neighbors only, not to central locus
            double phi_Lj = 1.0; // Default: include
            if (distance > 0) { // If not self (neighbor)
                phi_Lj = (variant.likelihood_score > likelihood_threshold) ? 1.0 : 0.0;
            } else { // Self contribution (dii = 0)
                // Include central locus regardless of threshold to stabilize Si
                phi_Lj = 1.0;
            }
            
            // Add to numerator and denominator
            weighted_sum += gaussian_weight * phi_Lj * variant.likelihood_score;
            weight_sum += gaussian_weight * phi_Lj;
        }
        
        // Return normalized density: Si = numerator / denominator
        return (weight_sum > 0.0) ? (weighted_sum / weight_sum) : 0.0;
    }

    // Calculate regional density for all variants in a chromosome region
    void calculate_regional_densities(vector<variant_locus>& variants) {
        
        // Sort variants by chromosome and position for efficient processing
        sort(variants.begin(), variants.end(), 
             [](const variant_locus& a, const variant_locus& b) {
                 if (a.chromosome != b.chromosome) return a.chromosome < b.chromosome;
                 return a.position < b.position;
             });
        
        // EDGE CASE CHECK: Verify likelihood thresholds before processing
        if (!check_likelihood_threshold_edge_case(variants)) {
            // Set all densities to NaN if all below threshold
            for (auto& variant : variants) {
                variant.regional_density = std::numeric_limits<double>::quiet_NaN();
            }
            return;
        }

        // Group variants by chromosome for adaptive bandwidth calculation
        map<string, double> chromosome_bandwidths;
        if (adaptive_bandwidth) {
            set<string> chromosomes;
            for (const auto& variant : variants) {
                chromosomes.insert(variant.chromosome);
            }
            
            for (const string& chrom : chromosomes) {
                // EDGE CASE CHECK: Skip empty chromosomes
                if (!handle_density_calculation_edge_cases(variants, chrom)) {
                    continue;
                }
                
                chromosome_bandwidths[chrom] = calculate_adaptive_bandwidth(variants, chrom);
                cout << "Adaptive bandwidth for " << chrom << ": " << chromosome_bandwidths[chrom] << " bp" << endl;
            }
        }
        
        // For each variant, calculate its regional density
        for (size_t i = 0; i < variants.size(); ++i) {
            
            // Use adaptive or fixed bandwidth
            double effective_sigma = adaptive_bandwidth ? 
                chromosome_bandwidths[variants[i].chromosome] : bandwidth_sigma;
            
            // EDGE CASE HANDLING: Clamp sigma for chromosomes
            // If σ > chrom_length/6, clamp to chrom_length/6 so 3σ fits
            effective_sigma = handle_bandwidth_edge_cases(variants, variants[i].chromosome, effective_sigma);
            
            // Define search window around current variant (±3σ for 99.7% coverage)
            int window_size = static_cast<int>(3.0 * effective_sigma);
            int min_pos = variants[i].position - window_size;
            int max_pos = variants[i].position + window_size;
            
            // Collect variants within the window on same chromosome
            vector<variant_locus> window_variants;
            for (const auto& variant : variants) {
                if (variant.chromosome == variants[i].chromosome && 
                    variant.position >= min_pos && variant.position <= max_pos) {
                    window_variants.push_back(variant);
                }
            }
            
            // Temporarily update bandwidth_sigma for this calculation
            double original_sigma = bandwidth_sigma;
            bandwidth_sigma = effective_sigma;
            
            // Calculate regional density for this variant
            variants[i].regional_density = calculate_regional_density(variants[i], window_variants);
            
            // Restore original bandwidth
            bandwidth_sigma = original_sigma;
        }
    }

    // Apply dual threshold filtering (frequency + regional density)
    void apply_dual_threshold_filter(vector<variant_locus>& variants,
                                   double min_frequency, double max_frequency) {
        
        // If regional_density_threshold is set to auto-calculate (e.g., -1), 
        // use robust default: median(Si) + 1 × MAD(Si)
        double density_thresh = regional_density_threshold;
        if (density_thresh < 0) {
            density_thresh = calculate_robust_density_threshold(variants);
        }
        
        for (auto& variant : variants) {
            bool freq_pass = (variant.frequency >= min_frequency && 
                            variant.frequency <= max_frequency);
            // Note: density threshold must match the scale of Li
            bool density_pass = (variant.regional_density >= density_thresh);
            
            variant.passes_criteria = freq_pass && density_pass;
        }
    }

    // Calculate adaptive bandwidth: σ = clamp(D95/2, 5kb, 50kb)
    double calculate_adaptive_bandwidth(const vector<variant_locus>& variants, const string& chromosome) {
        vector<int> distances;
        
        // Collect inter-hit distances for variants passing frequency filters on this chromosome
        vector<variant_locus> chrom_variants;
        for (const auto& variant : variants) {
            if (variant.chromosome == chromosome) {
                chrom_variants.push_back(variant);
            }
        }
        
        if (chrom_variants.size() < 2) {
            return bandwidth_sigma; // Fall back to default
        }
        
        // Sort by position
        sort(chrom_variants.begin(), chrom_variants.end(), 
             [](const variant_locus& a, const variant_locus& b) {
                 return a.position < b.position;
             });
        
        // Calculate consecutive distances
        for (size_t i = 1; i < chrom_variants.size(); ++i) {
            int distance = chrom_variants[i].position - chrom_variants[i-1].position;
            distances.push_back(distance);
        }
        
        if (distances.empty()) return bandwidth_sigma;
        
        // Calculate 95th percentile (D95)
        sort(distances.begin(), distances.end());
        size_t p95_index = (size_t)(0.95 * distances.size());
        if (p95_index >= distances.size()) p95_index = distances.size() - 1;
        
        double d95 = distances[p95_index];
        double adaptive_sigma = d95 / 2.0;
        
        // Clamp to [5kb, 50kb] range
        adaptive_sigma = max(5000.0, min(50000.0, adaptive_sigma));
        
        return adaptive_sigma;
    }

    // Identify contiguous segments from filtered variants
    vector<genomic_segment> identify_segments(const vector<variant_locus>& variants) {
        
        vector<genomic_segment> segments;
        vector<variant_locus> passing_variants;
        
        // Collect only variants that pass dual criteria
        for (const auto& variant : variants) {
            if (variant.passes_criteria) {
                passing_variants.push_back(variant);
            }
        }
        
        if (passing_variants.empty()) return segments;
        
        // Sort by chromosome and position
        sort(passing_variants.begin(), passing_variants.end(), 
             [](const variant_locus& a, const variant_locus& b) {
                 if (a.chromosome != b.chromosome) return a.chromosome < b.chromosome;
                 return a.position < b.position;
             });
        
        // Build segments
        genomic_segment current_segment(passing_variants[0].chromosome, 
                                       passing_variants[0].position,
                                       passing_variants[0].position);
        current_segment.variants.push_back(passing_variants[0]);
        
        for (size_t i = 1; i < passing_variants.size(); i++) {
            const auto& variant = passing_variants[i];
            
            // Same chromosome and within merge distance
            if (variant.chromosome == current_segment.chromosome &&
                (variant.position - current_segment.end_position) <= segment_merge_distance) {
                
                current_segment.variants.push_back(variant);
                current_segment.end_position = variant.position;
                
            } else {
                // Finalize current segment if it meets minimum variant requirement
                if (current_segment.variants.size() >= min_variants_per_segment) {
                    finalize_segment_statistics(current_segment);
                    segments.push_back(current_segment);
                }
                
                // Start new segment
                current_segment = genomic_segment(variant.chromosome, variant.position, variant.position);
                current_segment.variants.push_back(variant);
            }
        }
        
        // Add final segment
        if (current_segment.variants.size() >= min_variants_per_segment) {
            finalize_segment_statistics(current_segment);
            segments.push_back(current_segment);
        }
        
        return segments;
    }

private:
    // Calculate robust default threshold: median(Si) + 1 × MAD(Si)
    // Uses only frequency-filtered loci (but before density pass) as per design decision 6
    double calculate_robust_density_threshold(const vector<variant_locus>& variants) {
        vector<double> densities;
        
        // Collect densities from all variants (these are already frequency-filtered in sequential hierarchy)
        for (const auto& variant : variants) {
            densities.push_back(variant.regional_density);
        }
        
        if (densities.empty()) return 0.1; // Fallback default
        
        // Handle sparse chromosomes (Y/MT) - fallback to fixed threshold if sample count < 200
        if (densities.size() < 200) {
            cout << "Warning: Sparse chromosome detected (" << densities.size() << " variants). Using fixed threshold." << endl;
            return 0.1; // or consider widening σ
        }
        
        sort(densities.begin(), densities.end());
        
        // Calculate median
        double median = (densities.size() % 2 == 0) ? 
            (densities[densities.size()/2 - 1] + densities[densities.size()/2]) / 2.0 :
            densities[densities.size()/2];
        
        // Calculate MAD (Median Absolute Deviation)
        vector<double> abs_deviations;
        for (double density : densities) {
            abs_deviations.push_back(abs(density - median));
        }
        sort(abs_deviations.begin(), abs_deviations.end());
        
        double mad = (abs_deviations.size() % 2 == 0) ? 
            (abs_deviations[abs_deviations.size()/2 - 1] + abs_deviations[abs_deviations.size()/2]) / 2.0 :
            abs_deviations[abs_deviations.size()/2];
        
        // Handle MAD=0 edge case
        if (mad == 0.0) {
            cout << "Warning: MAD=0 detected. Using fallback threshold: median + 0.1" << endl;
            return median + 0.1; // Fixed epsilon in Li's units
        }
        
        return median + 1.0 * mad;
    }

    // Handle bandwidth edge cases for chromosomes
    double handle_bandwidth_edge_cases(const vector<variant_locus>& variants, 
                                      const string& chromosome, double sigma) {
        // Find chromosome length by getting max position
        int max_pos = 0;
        for (const auto& variant : variants) {
            if (variant.chromosome == chromosome) {
                max_pos = max(max_pos, variant.position);
            }
        }
        
        if (max_pos == 0) return sigma; // No variants found
        
        // If σ > chrom_length/6, clamp to chrom_length/6 so 3σ fits
        double max_sigma = max_pos / 6.0;
        if (sigma > max_sigma) {
            cout << "Warning: Bandwidth " << sigma << " too large for chromosome " 
                 << chromosome << " (length ~" << max_pos << "). Clamping to " << max_sigma << endl;
            return max_sigma;
        }
        
        return sigma;
    }

    // Handle edge cases during density calculation
    bool handle_density_calculation_edge_cases(const vector<variant_locus>& variants, 
                                              const string& chromosome) {
        // Edge case 1: Empty chromosomes
        bool found_variants = false;
        for (const auto& variant : variants) {
            if (variant.chromosome == chromosome) {
                found_variants = true;
                break;
            }
        }
        
        if (!found_variants) {
            cout << "Warning: Empty chromosome " << chromosome << " detected. Skipping density calculation." << endl;
            return false;
        }
        
        // Edge case 2: Very sparse data (< 3 variants)
        int variant_count = 0;
        for (const auto& variant : variants) {
            if (variant.chromosome == chromosome) {
                variant_count++;
            }
        }
        
        if (variant_count < 3) {
            cout << "Warning: Very sparse chromosome " << chromosome << " (" << variant_count 
                 << " variants). Density calculation may be unreliable." << endl;
        }
        
        return true;
    }

    // Check for all variants below likelihood threshold
    bool check_likelihood_threshold_edge_case(const vector<variant_locus>& variants) {
        int above_threshold = 0;
        for (const auto& variant : variants) {
            if (variant.likelihood_score > likelihood_threshold) {
                above_threshold++;
            }
        }
        
        if (above_threshold == 0) {
            cout << "Error: All variants below likelihood threshold " << likelihood_threshold 
                 << ". Regional density will be undefined (NaN). Consider lowering threshold." << endl;
            return false;
        }
        
        if (above_threshold < 10) {
            cout << "Warning: Only " << above_threshold << " variants above likelihood threshold. "
                 << "Density calculation may be unstable." << endl;
        }
        
        return true;
    }

private:
    // Calculate summary statistics for a segment
    void finalize_segment_statistics(genomic_segment& segment) {
        segment.variant_count = segment.variants.size();
        
        if (segment.variant_count == 0) return;
        
        // Calculate mean likelihood and density
        double likelihood_sum = 0.0;
        double density_sum = 0.0;
        
        for (const auto& variant : segment.variants) {
            likelihood_sum += variant.likelihood_score;
            density_sum += variant.regional_density;
        }
        
        segment.mean_likelihood = likelihood_sum / segment.variant_count;
        segment.mean_density = density_sum / segment.variant_count;
    }
};

#endif // SEGMENT_ANALYSIS_MODULE_H