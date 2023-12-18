from ase.phasediagram import PhaseDiagram
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import re 
from general_lib_37 import *
import csv
from a4_phase_diag import Individuals2csv, create_phasediag


def create_enthalpies_plot(filename):
	df = pd.read_csv(filename, index_col="ID")
	ele1, ele2, ele3 = "Sm", "Fe", "Ti"
	eles = [ele1, ele2, ele3]
	refs = []

	fig = plt.figure(figsize=(8, 8))

	indexes = df.index
	for idx in indexes:
		try:
			enthalpy_values = list(map(float, df.loc[idx, :].values))
			if len(enthalpy_values) == 4:
				plt.plot(enthalpy_values, linestyle='-.', alpha=0.8)
			plt.annotate(idx, xy=(0, enthalpy_values[0]), size=12)
		except Exception as e:
			pass
	plt.ylabel('enthalpy')
	plt.xlabel('vr step')
	# plt.ylim([-130, -90])

	# plt.title(title, **title_font)
	plt.tight_layout(pad=1.1)

	# # # setting only for Individuals # # #
	saveat = filename.replace("/uspex_calc/","/figure/").replace(".csv",".pdf")

	makedirs(saveat)
	plt.savefig(saveat, transparent=False)

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
				magmom = get_magmom_v(filename=mmp)
				n_atom_has_mag = len(magmom)
				ene_df.loc[indexes[:n_atom_has_mag],"magmom"] = magmom * vol[:n_atom_has_mag]
				ene_df.loc[indexes[:n_atom_has_mag],"magmom_pv"] = magmom
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


	if is_oqmd_pd: 
		for ele in ["Zr", "N", "As", "P", "Ti"]: # "Co", "Cu", "Ga", "B"
			comps_array = copy.copy(org_comp)
			comps_array.append(ele)
			cmp_name = "".join(comps_array)

			# # structures from original oqmd 
			oqmd_file = input_dir + "/oqmd_csv/{}.csv".format(cmp_name)
			savedir = result_dir + "/figure/oqmd_csv/{}".format(cmp_name)
			df = pd.read_csv(oqmd_file, index_col="ID")

			df = pd.read_csv(oqmd_file, index_col=0)
			create_phasediag(filename=oqmd_file, 
				comps_array=comps_array, comps_list=["Composition{}".format(k) for k in comps_array],
				savedir=savedir, phase_corners=[[1, 0, 0], [0, 1, 0] , [0, 0, 1]],
				z_values=df["energy_substance_pa"].values, vcenter=0.0,
				v_range_dict=dict({"delta_e": [-0.1, 0.1],
					"energy_substance_pa": [-0.1, 0.1]}),
				is_show_all=True)











































