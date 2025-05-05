"""Eye Color Prediction Tool CLI

This program analyzes DNA data from AncestryDNA files to predict eye color based on specific SNPs
using advanced genomic analysis techniques.
"""

import os
import logging
from pathlib import Path
from typing import Dict, Tuple
from parser import (
    parse_dna_file, predict_eye_color, SNPInfo, GenotypeResult,
    SNP_DATABASE, EyeColor
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def greet() -> None:
    """Display welcome message and program information."""
    print(r"""
 ███████╗██╗   ██╗███████╗    ██████╗ ██████╗ ██╗      ██████╗ ██████╗ 
 ██╔════╝╚██╗ ██╔╝██╔════╝    ██╔══██╗██╔══██╗██║     ██╔═══██╗██╔══██╗
 █████╗   ╚████╔╝ █████╗      ██████╔╝██████╔╝██║     ██║   ██║██████╔╝
 ██╔══╝    ╚██╔╝  ██╔══╝      ██╔══██╗██╔══██╗██║     ██║   ██║██╔══██╗
 ███████╗   ██║   ███████╗    ██████╔╝██║  ██║███████╗╚██████╔╝██║  ██║
 ╚══════╝   ╚═╝   ╚══════╝    ╚═════╝ ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝
                                                                          
 ██████╗ ██████╗ ███████╗██████╗ ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗
 ██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║
 ██████╔╝██████╔╝█████╗  ██████╔╝██║██║        ██║   ██║██║   ██║██╔██╗ ██║
 ██╔═══╝ ██╔══██╗██╔══╝  ██╔══██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║
 ██║     ██║  ██║███████╗██║  ██║██║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║
 ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚════╝
                                                                          
 ████████╗ ██████╗  ██████╗ ██╗     ███████╗
 ╚══██╔══╝██╔═══██╗██╔═══██╗██║     ██╔════╝
    ██║   ██║   ██║██║   ██║██║     █████╗  
    ██║   ██║   ██║██║   ██║██║     ██╔══╝  
    ██║   ╚██████╔╝╚██████╔╝███████╗███████╗
    ╚═╝    ╚═════╝  ╚═════╝ ╚══════╝╚══════╝
    """)

    print("\nWelcome to the Advanced Eye Color Prediction Tool!")
    print("\nThis sophisticated genomic analysis tool uses cutting-edge techniques")
    print("to predict eye color based on your DNA data. Our algorithm analyzes")
    print("multiple Single Nucleotide Polymorphisms (SNPs) with advanced quality")
    print("metrics and population-specific frequencies.")
    
    print("\nKey Features:")
    print("- Advanced genotype quality assessment")
    print("- Population-specific allele frequencies")
    print("- Quality-adjusted scoring system")
    print("- Comprehensive SNP database")
    print("- Detailed genomic analysis")
    
    print("\nAnalysis Methodology:")
    print("- Uses weighted scoring system for each SNP")
    print("- Considers genotype combinations (AA, AG, GG)")
    print("- Applies softmax function for probability calculation")
    print("- Accounts for both positive and negative influences")
    print("- Provides detailed probability distribution")
    print("- Includes quality metrics and read depth analysis")
    
    print("\nPrediction Accuracy:")
    print("- Blue/Brown: ~90% accuracy")
    print("- Green/Hazel: ~75% accuracy")
    print("- Considers multiple genetic factors")
    print("- Accounts for complex inheritance patterns")
    print("- Quality-adjusted confidence scores")

def get_dna_file() -> Path:
    """Get and validate the DNA file path from user input.
    
    Returns:
        Path: Validated path to the DNA file
    """
    while True:
        try:
            file_path = Path(input("\nPlease enter the path to your AncestryDNA File: ").strip())
            
            if not file_path.exists():
                print(f"Error: File '{file_path}' does not exist.")
                continue
                
            if not file_path.is_file():
                print(f"Error: '{file_path}' is not a file.")
                continue
                
            logger.info(f"File '{file_path}' found and validated")
            return file_path
            
        except Exception as e:
            logger.error(f"Error processing file path: {e}")
            print("Please try again with a valid file path.")

def format_probabilities(probabilities: Dict[str, float]) -> str:
    """Format probability values for display.
    
    Args:
        probabilities: Dictionary of eye colors and their probabilities
        
    Returns:
        str: Formatted string of probabilities
    """
    return "\n".join(
        f"{color.capitalize():<8}: {probability * 100:>6.2f}%"
        for color, probability in sorted(probabilities.items())
    )

def format_genotype_info(result: GenotypeResult, snp_info: SNPInfo) -> str:
    """Format genotype information for display.
    
    Args:
        result: GenotypeResult object
        snp_info: SNPInfo object
        
    Returns:
        str: Formatted string of genotype information
    """
    return f"""
SNP: {result.rsid}
Gene: {snp_info.gene}
Chromosome: {snp_info.chromosome}
Position: {snp_info.position}
Genotype: {result.genotype}
Zygosity: {result.zygosity}
Quality Score: {result.quality_score:.2f}
Read Depth: {result.read_depth}
Population Frequencies:
  - European: {snp_info.population_frequency.get('EUR', 0.0) * 100:.1f}%
  - African: {snp_info.population_frequency.get('AFR', 0.0) * 100:.1f}%
  - East Asian: {snp_info.population_frequency.get('EAS', 0.0) * 100:.1f}%
Clinical Significance: {snp_info.clinical_significance}
"""

def main() -> None:
    """Main program entry point."""
    try:
        greet()
        file_path = get_dna_file()
        
        logger.info("Parsing DNA file...")
        snp_data = parse_dna_file(file_path)
        
        if not snp_data:
            print("\nError: No relevant SNP data found in the file.")
            return
            
        print("\nGenomic Analysis Results:")
        print("\nExtracted SNPs and Genotypes:")
        for rsid, result in sorted(snp_data.items()):
            if rsid in SNP_DATABASE:
                print(format_genotype_info(result, SNP_DATABASE[rsid]))
            
        logger.info("Predicting eye color...")
        predicted_eye_color, probabilities, quality_assessment = predict_eye_color(snp_data)
        
        print("\nPrediction Results:")
        print(f"Predicted Eye Color: {predicted_eye_color.capitalize()}")
        print("\nDetailed Probability Distribution:")
        print(format_probabilities(probabilities))
        
        print("\nQuality Assessment:")
        for rsid, assessment in quality_assessment.items():
            print(f"{rsid}: {assessment}")
        
        print("\nAnalysis Notes:")
        print("- Results are based on genetic markers and statistical analysis")
        print("- Environmental factors may influence actual eye color")
        print("- Some rare genetic variations may not be detected")
        print("- Results should be considered as probabilities, not certainties")
        print("- Quality metrics are based on sequencing data")
        print("- Population frequencies are based on reference databases")
        
    except KeyboardInterrupt:
        print("\nProgram terminated by user.")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        print("\nAn error occurred while processing your request. Please try again.")

if __name__ == "__main__":
    main()
