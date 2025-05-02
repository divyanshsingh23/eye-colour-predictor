"""cli entry point for the code
"""
from parser import parse_dna_file , predict_eye_color 
import os

def greet():

    print(r""" 
 _______ ____    ____  _______    .______   .______       _______  _______   __    ______ .___________. __    ______   .__   __. 
|   ____|\   \  /   / |   ____|   |   _  \  |   _  \     |   ____||       \ |  |  /      ||           ||  |  /  __  \  |  \ |  | 
|  |__    \   \/   /  |  |__      |  |_)  | |  |_)  |    |  |__   |  .--.  ||  | |  ,----'`---|  |----`|  | |  |  |  | |   \|  | 
|   __|    \_    _/   |   __|     |   ___/  |      /     |   __|  |  |  |  ||  | |  |         |  |     |  | |  |  |  | |  . `  | 
|  |____     |  |     |  |____    |  |      |  |\  \----.|  |____ |  '--'  ||  | |  `----.    |  |     |  | |  `--'  | |  |\   | 
|_______|    |__|     |_______|   | _|      | _| `._____||_______||_______/ |__|  \______|    |__|     |__|  \______/  |__| \__| 
                                                                                                                       """)

    
    print("Welcome to the Eye Color Prediction Tool!")
    print("This program takes into consideration the following SNPs to predict eye color:")
    print("- rs12913832 (HERC2)")
    print("- rs1800407 (OCA2)")
    print("- rs1129038 (HERC2)")
    print("- rs12203592 (IRF4)")
    print("- rs16891982 (SLC45A2)")
    print("- rs1426654 (SLC24A5)")
    print("- rs11636232 (TYRP1)")
    print("- rs12896399 (SLC24A4)")
    print("- rs1393350 (TYR)")
    print("- rs4778138 (OCA2)")
    print("\nBy entering your AncestryDNA file, this tool will analyze the SNPs and predict your likely eye color.")
    print("\nLet's get started! ")

def get_dna_file():
    import os 

    while True:
        file_path = input("Please enter the path to your AncestryDNA File")

        if os.path.isfile(file_path):
            print("File found")
            return file_path
        else:
            print("File not found")


def main():
    greet()
    
    file_path = input("Please enter the path to your AncestryDNA File: ").strip()

    if not os.path.exists(file_path):
        print("Error: File does not exist.")
        return

    print("File found")
    print(f"File '{file_path}' has been successfully loaded into the program!")

    snp_data = parse_dna_file(file_path)

    print("\nExtracted Eye-Color Related SNPs:")
    for rsid, genotype in snp_data.items():
        print(f"{rsid}: {genotype}")

    predicted_eye_color, probabilities = predict_eye_color(snp_data)

    print("\nPredicted Eye Color: ", predicted_eye_color)
    print("Probabilities for each eye color:")
    for color, probability in probabilities.items():
        print(f"{color.capitalize()}: {probability * 100:.2f}%")
  
if __name__ == "__main__":
    main()
