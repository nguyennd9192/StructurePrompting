import pandas as pd
import numpy as np

if __name__ == '__main__':
    # Load dataset of materials structures
    
    
    for job in jobs:
		ofm_ene_file = get_ofm_ene_csv_dir(job=job)
		idv_file = get_uspex_dir(job=job, uspex_file="results1/Individuals_magmom_ene.csv")
		save_fig_at = get_uspex_dir(job=job, uspex_file="ML/uspex_population")
		symmetry_seeding(job=job, uspex_file=ofm_ene_file, 
				idv_file=idv_file, save_fig_at=save_fig_at)