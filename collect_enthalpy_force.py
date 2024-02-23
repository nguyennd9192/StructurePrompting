from lib import *
from ase.phasediagram import PhaseDiagram
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import re 
if __name__ == '__main__':
	is_uspex_pd = True
	is_oqmd_pd = False

	if is_uspex_pd:
		for job in jobs:
			eles = add_element(job)

			idv = get_uspex_dir(job=job, uspex_file="results1/Individuals")
			idv_csv, comps_array, comps_list, phase_corners = Individuals2csv(filename=idv)
			origin_file = get_uspex_dir(job=job, uspex_file="results1/origin")
			origin2csv(filename=origin_file)
			
			ene_df = pd.read_csv(idv_csv, index_col="ID")
			indexes = ene_df.index
			vol = ene_df["Volume"].values
			natom = ene_df[comps_list].sum(axis=1)

			mmp = result_dir + "{0}/results1/magmomProperties".format(job)
			if os.path.isfile(mmp):
				# # to get magmom_pa
				magmom_pv = get_magmom_v(filename=mmp)
				n_atom_has_mag = len(magmom)
				

				magmom 	= magmom_pv * vol[:n_atom_has_mag]
				for index in indexes[:n_atom_has_mag]:
					for element in comps_array:
						if element in LANTHANIDE_MAG.keys():
							n_lanthanide = ene_df.loc[index, "Composition{}".format(element)]
							magmom[index] += n_lanthanide*J4f[element]*gJ4f[element] # # 0.714 with Sm


				ene_df.loc[indexes[:n_atom_has_mag],"magmom"] = magmom
				ene_df.loc[indexes[:n_atom_has_mag],"magmom_pv"] = magmom_pv
				ene_df.loc[indexes[:n_atom_has_mag],"magmom_pa"] = ene_df["magmom"] / natom

			# # # to get energy_substance_pa
			energy_substance = copy.copy(ene_df["Enthalpy"].values)
			for comp_th, comp in enumerate(comps_array):
				# print (toten_simple[comp], ene_df["Composition{0}".format(comp)])
				energy_substance -= toten_simple[comp] * ene_df["Composition{0}".format(comp)]
			ene_df["energy_substance_pa"] = energy_substance / natom

			# # # save and plot phase diag with properties
			magmom_ene_csv = get_uspex_dir(job=job, uspex_file="results1/Individuals_magmom_ene.csv")
			ene_df.to_csv(magmom_ene_csv) 


			# # get ofm file
			ofm_file = get_ofm_csv_dir(job=job)
			ofm_df = pd.read_csv(ofm_file, index_col="name")
			n_ofm = len(ofm_df)

			ID = [norm_uspex_id(full_path=k, rmvs=["/symstruct/"]) for k in ofm_df.index]
			ofm_df["ID"] = ID

			# # # merge and  save df 
			ofm_ene_file = get_ofm_ene_csv_dir(job=job)

			makedirs(ofm_ene_file)
			merge_df = ene_df.merge(ofm_df, how='inner', on='ID')
			save_symm = copy.copy(merge_df["SYMM"].values)

			merge_df = merge_df.loc[:, (merge_df != merge_df.iloc[0]).any()]
			remove_cols = ["Unnamed: 0", "saveat"] # "Origin", "Gen"
			for c in remove_cols:
				merge_df.drop(c, axis=1, inplace=True)
			merge_df = merge_df.set_index('ID')
			merge_df = merge_df[~merge_df.index.duplicated(keep='first')]
			merge_df["SYMM"] = save_symm
			assert len(ofm_df) == len(ene_df) == len(merge_df)

			merge_df.to_csv(ofm_ene_file)
