# StructurePrompting
=======

Title: Structure Prompting with USPEX program

I. Introduction
The framework begins with a recommendation system with a backbone from Dempster- Shafer theory to estimate configurations constructed by chemical composition and crystal symmetry which are high potentially stable in nature. A genetic algorithm-powered structure generator to prompt hypothetical skeletons that match the recommended configurations and highly expected to optimize multi-geometrical and physical properties contraints. Finally, high-potential skeletons will be validated for physical properties using first-principle calculations.



II. Installation
1. Create a conda virtual environment:
conda create --name <env> --file requirements.txt


2. Install Universal Structure Predictor: Evolutionary Xtallography (USPEX)
Please follow instruction from USPEX program in:
https://uspex-team.org/online_utilities/uspex_manual_release/EnglishVersion/uspex_manual_english/sect0011.html


3. Install Vienna Ab initio Simulation Package (VASP)
Please follow instruction to setting VASP calculation in:
https://www.vasp.at/wiki/index.php/Installing_VASP.5.X.X


III. Usage
1. Dempster-Shafer Theory based Recommendation System  
Explain the use of Dempster-Shafer theory in the recommendation system.



2. Structure Prompting 
Example of running Structure Prompting program is shown in "example" folder.
Bash script of running the example is shown in run.sh. Details description is shown as in follow.

# # # Save all generated structures.
python $code_dir/collect_structures.py example/Sm1Fe12

# # # Prepare Orbital-Field matrix representation for all generated structures.
python $code_dir/descriptor.py example/Sm1Fe12

# # # Estimate formation energy, saturated magnetization for all generated structures.
python $code_dir/process_features.py example/Sm1Fe12
python $code_dir/collect_enthalpy_force.py example/Sm1Fe12


3. Prompting Hypothetical Skeletons 
# # # Prompt new hypothetical structures in the next generation of USPEX.
python $code_dir/seeding.py example/Sm1Fe12


