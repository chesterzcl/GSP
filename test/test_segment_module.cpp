#include "../src/segment_analysis_module.h"
#include "../src/segment_output_module.h"
#include <iostream>
#include <cassert>

using namespace std;

// Test function for normalized Gaussian kernel density calculation
void test_gaussian_kernel_density() {
    cout << "Testing normalized Gaussian kernel density calculation..." << endl;
    
    segment_analysis_module analyzer(1000.0, 0.1, 5000.0, 2); // 1kb bandwidth
    
    // Create test variants
    vector<variant_locus> test_variants = {
        variant_locus("chr1", 1000, 0.8, 0.3, "test1"),  // High likelihood
        variant_locus("chr1", 1500, 0.6, 0.4, "test2"),  // Medium likelihood  
        variant_locus("chr1", 2000, 0.9, 0.5, "test3"),  // High likelihood
        variant_locus("chr1", 5000, 0.2, 0.1, "test4")   // Low likelihood, distant
    };
    
    // Test density calculation for middle variant
    variant_locus target = test_variants[1]; // Position 1500
    double density = analyzer.calculate_regional_density(target, test_variants);
    
    cout << "Target variant at position " << target.position << endl;
    cout << "Calculated regional density: " << density << endl;
    
    // Density should be > 0 since there are nearby variants with likelihood scores
    assert(density > 0.0);
    cout << "✓ Gaussian kernel density calculation test passed." << endl << endl;
}

// Test function for segment identification
void test_segment_identification() {
    cout << "Testing segment identification..." << endl;
    
    segment_analysis_module analyzer(1000.0, 0.3, 2000.0, 2); // Low density threshold
    
    // Create test variants with calculated densities
    vector<variant_locus> test_variants = {
        variant_locus("chr1", 1000, 0.8, 0.3, "test1"),
        variant_locus("chr1", 1200, 0.7, 0.4, "test2"),
        variant_locus("chr1", 1400, 0.9, 0.5, "test3"),
        variant_locus("chr1", 5000, 0.6, 0.2, "test4"), // Distant variant
        variant_locus("chr1", 5200, 0.8, 0.3, "test5")
    };
    
    // Calculate regional densities
    analyzer.calculate_regional_densities(test_variants);
    
    // Apply filtering (assume all pass frequency criteria)
    for (auto& variant : test_variants) {
        variant.passes_criteria = (variant.regional_density >= 0.3);
    }
    
    // Identify segments
    vector<genomic_segment> segments = analyzer.identify_segments(test_variants);
    
    cout << "Number of segments identified: " << segments.size() << endl;
    
    for (size_t i = 0; i < segments.size(); ++i) {
        cout << "Segment " << (i+1) << ": chr" << segments[i].chromosome 
             << ":" << segments[i].start_position 
             << "-" << segments[i].end_position
             << " (" << segments[i].variant_count << " variants)" << endl;
    }
    
    // Should identify at least one segment
    assert(segments.size() > 0);
    cout << "✓ Segment identification test passed." << endl << endl;
}

// Test function for output functionality
void test_output_functionality() {
    cout << "Testing output functionality..." << endl;
    
    segment_output_module output_module;
    
    // Create test data
    vector<variant_locus> test_variants = {
        variant_locus("chr1", 1000, 0.8, 0.3, "gene1|SNP|A>G"),
        variant_locus("chr1", 1200, 0.7, 0.4, "gene1|SNP|C>T"),
        variant_locus("chr1", 1400, 0.9, 0.5, "gene1|INDEL|DEL")
    };
    
    // Set regional densities and criteria
    test_variants[0].regional_density = 0.45;
    test_variants[0].passes_criteria = true;
    test_variants[1].regional_density = 0.52;
    test_variants[1].passes_criteria = true;
    test_variants[2].regional_density = 0.48;
    test_variants[2].passes_criteria = true;
    
    // Create test segment
    genomic_segment test_segment("chr1", 1000, 1400);
    test_segment.variants = test_variants;
    test_segment.variant_count = 3;
    test_segment.mean_likelihood = 0.8;
    test_segment.mean_density = 0.48;
    
    vector<genomic_segment> test_segments = {test_segment};
    vector<string> pop_names = {"EUR", "ASN", "AFR"};
    
    // Test output files creation
    string test_output = "test_output";
    output_module.write_all_outputs(test_output, test_variants, test_segments, pop_names);
    
    cout << "✓ Output functionality test completed." << endl;
    cout << "Check generated files: test_output.locus, test_output.segments, etc." << endl << endl;
}

// Main test function
int main() {
    cout << "=== GSP Segment-First Discovery Module Tests ===" << endl << endl;
    
    try {
        test_gaussian_kernel_density();
        test_segment_identification();
        test_output_functionality();
        
        cout << "=== All Tests Passed Successfully! ===" << endl;
        cout << "The segment-first discovery module is ready for integration." << endl;
        
    } catch (const exception& e) {
        cout << "Test failed with error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}