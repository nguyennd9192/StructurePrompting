import os, sys, qmpy, pickle, shutil

sc_dir = "/media/nguyen/work/SmFe12_screening/source_code/"
for ld, subdirs, files in os.walk(sc_dir):
	if os.path.isdir(ld) and ld not in sys.path:
		sys.path.append(ld)

import re, glob
import numpy as np
import pandas as pd
from general_lib import get_subdirs, get_basename, makedirs
from features import OFMFeature, XRDFeature
from plot import scatter_plot, scatter_plot_2, scatter_plot_3


from manifold_processing import Preprocessing
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.cif import CifWriter

import pymatgen as pm
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt


alljobs = ["Sm-Fe10-M2_wyckoff", "Sm-Fe11-M_wyckoff_2", 
			"Sm2-Fe21-M3_std6000",

			# "Sm2-Fe21-M3", # # temporal remove this job out of consideration, replacing by _std6000

			"Sm2-Fe22-M2", "Sm2-Fe23-M", 
			"check_Mo_2-22-2", 
			"CuAlZnTi_Sm2-Fe22-M2", "CuAlZnTi_Sm2-Fe21-M3"

			"Sm2-Fe23-Ga1", "Sm2-Fe22-Ga2", "Sm2-Fe21-Ga3",
			"Sm2-Fe23-Si1", "Sm2-Fe22-Si2", "Sm2-Fe21-Si3"
			]

prefix_local_idx = "_lcid_"


def get_local_vector(struct_dir_file, save_csv, n_neighbors):
	data = []
	with open(struct_dir_file) as file:
		struct_reprs = file.read().splitlines()

	for struct_repr in struct_reprs:
		with open(struct_repr, "rb") as tmp_f:
			struct_repr_obj = pickle.load(tmp_f) #  encoding='latin1'  
			data.append(struct_repr_obj)

	structs_locals = np.array([file.__getattribute__('locals')
							 for file in data])
	structs_locals_xyz = np.array([file.__getattribute__('local_xyz')
							 for file in data])
	structs_center_atoms = np.array([file.__getattribute__('atoms')
							 for file in data])
	# structs_wyckoff_all_sites = np.array([file.__getattribute__('wyckoff_all_sites')
	#                          for file in data])
	# structure_files = np.array([file.__getattribute__('name')
	#                          for file in data])
	feature_names = np.array(data[0].__getattribute__('feature_name'))
	X_all_locals = np.concatenate(structs_locals, axis=0)

	# # remove all zeros features
	non_zeros = np.where(X_all_locals.any(axis=0))[0]
	X_all_locals = X_all_locals[:, non_zeros]
	feature_names = feature_names[non_zeros]
	# # get local index
	all_local_idxes = []
	for ith, struct_repr in enumerate(struct_reprs):
		# # loop all structures
		struct_locals = structs_locals[ith]
		struct_locals_xyz = structs_locals_xyz[ith]

		# # read origin structure
		# struture_file = structure_files[ith]
		# print (struture_file)

		# with open(struture_file, 'r') as f:
		#     inp = f.read()
		# struct = pm.Structure.str(input_string=inp, fmt='poscar')
		# origin_lattice =  struct.lattice
		
		# print origin_lattice

		n_locals = struct_locals.shape[0]
		locals_idxes = []
		
		# # loop all locals in this structure
		for k in range(n_locals):
			this_idx = "{0}{1}{2}".format(struct_repr, prefix_local_idx, k)
			this_local_xyz = struct_locals_xyz[k]
			locals_idxes.append(this_idx)

			# this_struct =  pm.Structure(lattice=origin_lattice,
			#     species=this_local_xyz["atoms"][:n_neighbors], 
			#     coords=this_local_xyz["coords"][:n_neighbors])
			# poscar = Poscar(structure=this_struct)
			# poscar.write_file("test.poscar")
		all_local_idxes = np.concatenate([all_local_idxes, locals_idxes], axis=0)
	
	# # save to csv
	df = pd.DataFrame(X_all_locals, 
		index=all_local_idxes, columns=feature_names)

	print np.concatenate(structs_center_atoms)
	print len(np.ravel(structs_center_atoms)), len(all_local_idxes)
	df["center"] = np.concatenate(structs_center_atoms) # concatenate is ok with Nd-Fe-B
	# df["wyckoff_all_sites"] = np.ravel(structs_wyckoff_all_sites)

	print all_local_idxes
	df = df.dropna()
	makedirs(save_csv)
	df.to_csv(save_csv)
	print "(n_structs, n_locals, ofm_dims):", X_all_locals.shape
	print "Number of local indexes: ", all_local_idxes.shape


def get_local_xyz(local_adress, n_neighbors, saveat):
	print "local_adress: ", local_adress

	# # last digit is local index
	local_idx = local_adress[local_adress.find(prefix_local_idx) + len(prefix_local_idx):]

	# # the rest is structure representation object
	struct_repr = local_adress[: local_adress.find(prefix_local_idx)]

	# # read structure representation e.g. ABC.ofm_no_d
	with open(struct_repr, "rb") as tmp_f:
	  struct_repr_obj = pickle.load(tmp_f) #  encoding='latin1'


	# # get structure dir
	structure_file = struct_repr_obj.__getattribute__('name')
	# # get locals
	locals_str = struct_repr_obj.__getattribute__('locals')
	all_locals_xyz = struct_repr_obj.__getattribute__('local_xyz')
	center_atoms = struct_repr_obj.__getattribute__('atoms')


	locals_xyz = all_locals_xyz[int(local_idx)]

	n_locals = locals_str.shape[0]

	print (structure_file)
	# read origin structure
	with open(structure_file, 'r') as f:
	  inp = f.read()
	struct = pm.Structure.from_str(input_string=inp, fmt='poscar')
	origin_lattice =  struct.lattice


	# # sav to local structure file, center of local denoted by Neon
	export_species = locals_xyz["atoms"] #[:n_neighbors]
	print export_species[0], center_atoms[int(local_idx)]
	
	export_species[0] = "C"

	if export_species[0] == "Nd":
	  export_species[0] = "N"
	elif export_species[0] == "Fe":
	  export_species[0] = "F"
	elif export_species[0] == "B":
	  export_species[0] = "O"

	this_struct =  pm.Structure(lattice=origin_lattice,
		species=export_species, 
		coords=locals_xyz["coords"], # [:n_neighbors]
		to_unit_cell=False,
		coords_are_cartesian=True) # locals_xyz stored in cartesian
	# poscar = Poscar(structure=this_struct)
	cif = CifWriter(struct=this_struct)

		
	saveat = saveat[:saveat.find("id_ene_f_")]
	saveat.replace("from_", "")
	saveat +=".cif"
	makedirs(saveat)
	# poscar.write_file(saveat)
	cif.write_file(saveat)


def get_sub_idx(a, b, c, d, idx):
  ax = a[idx]
  bx = b[idx]
  cx = [round(k, 3) for k in c[idx]]

  dx = [round(k, 3) for k in d[idx]]

  return ax, bx, cx, dx


def processing_map(X, config, title, saveat, save_localdir, 
			full_local_idx, annot_name, color_array,
			marker_array, Xp=None):
	# # manifold learning

	processing = Preprocessing()
	# preprocessor = processing.__getattribute__("")
	processing.similarity_matrix = X
	
	# # ['iso_map', 'standard', 'locallyLinearEmbedding',
	# # 'hessianEigenmapping', 'tsne', 'mds']

	# X_trans_non_normalize, _, dim_reduc_method = processing.iso_map(**config_isomap)
	X_trans_non_normalize, _, dim_reduc_method, ranking_results = processing.tsne(**config)


	# # normalize dim1, dim2 of X transform
	scaler = MinMaxScaler()
	X_trans = scaler.fit_transform(X_trans_non_normalize)

	ntrains = X_trans.shape[0]
	# # fit and transform

	if Xp is not None:
		# # dim_reduc_method: dimensional reduction method, fit already with X
		Xp_trans_non_normalize = dim_reduc_method.transform(Xp)
		nprojects = Xp_trans_non_normalize.shape[0]

		# # scaler transform the Xp projecting matrix, then concatenate
		X_trans = np.concatenate([X_trans, scaler.transform(Xp_trans_non_normalize)])
		color_array = np.concatenate([color_array, ["blue"]*nprojects])
		save_localdir += ".project"

	
	x = X_trans[:, 0]
	y = X_trans[:, 1]


	xlb, xub = 0.25, 0.75
	ylb, yub = 0.25, 0.75
	n_neighbors = 7

	x = X_trans[:, 0]
	y = X_trans[:, 1]
	scatter_plot_3(x=x, y=y, 
			# xvlines=[xlb, xub], yhlines=[ylb, yub], 
			xvlines=None, yhlines=None, 
			s=100, alpha=0.2, 
			# title=title,
			sigma=None, mode='scatter', 
			name=None,  # all_local_idxes
			x_label='Dim 1', y_label="Dim 2", 
			save_file=saveat,
			interpolate=False, color_array=color_array, 
			preset_ax=None, linestyle='-.', marker=marker_array)

	scatter_plot_2(x=x, y=y, 
			# xvlines=[xlb, xub], yhlines=[ylb, yub], 
			xvlines=None, yhlines=None, 
			s=80, alpha=0.2, 
			# title=title,
			sigma=None, mode='scatter', 
			name=None,  # all_local_idxes
			x_label='Dim 1', y_label="Dim 2", 
			save_file=saveat.replace(".pdf", "mix.pdf"),
			interpolate=False, color_array=color_array, 
			preset_ax=None, linestyle='-.', marker=marker_array)
	
	# concern_x_idx = np.where( x > 0 and x < 10 )[0]
	if save_localdir is not None:
	  extension = "cif" # poscar
	  color_array = np.array(color_array)

	  this_idx = (x < xlb) & (y < ylb)

	  # local_xyz_saves = ["{0}/region11/{1}.{2}".format(save_localdir, get_basename(k), extension) for k in concern_locals]
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region11/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)

	  # local_xyz_saves = ["{0}/region12/{1}.{2}".format(save_localdir, get_basename(k), extension) for k in concern_locals]
	  this_idx = (x > xlb) & (x < xub) & (y < ylb)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region12/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)

	  # local_xyz_saves = ["{0}/region13/{1}.{2}".format(save_localdir, get_basename(k), extension) for k in concern_locals]

	  this_idx = (x > xub) & (y < ylb)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region13/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)


	  # local_xyz_saves = ["{0}/region21/{1}.{2}".format(save_localdir, get_basename(k), extension) for k in concern_locals]
	  this_idx = (y > ylb) & (y < yub) & (x < xlb)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region21/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)


	  this_idx = (y > ylb) & (y < yub) & (x > xlb) & (x < xub)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region22/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)


	  this_idx = (y > ylb) & (y < yub) & (x > xub)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region23/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)


	  this_idx = (y > yub) & (x < xlb)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region31/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)


	  this_idx = (y > yub) & (x > xlb) & (x < xub)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region32/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)


	  this_idx = (y > yub) & (x > xub)
	  concern_locals, color_code, xs, ys = get_sub_idx(a=full_local_idx, b=color_array, 
									  c=x, d=y, idx=this_idx)
	  local_xyz_saves = ["{0}/region33/x{1}-y{2}.{3}.{4}.{5}".format(
		save_localdir, xs[ith], ys[ith], color_code[ith], get_basename(k),  extension) 
			  for ith, k in enumerate(concern_locals)]
	  for local_adress, saveat in zip(concern_locals, local_xyz_saves):
		  get_local_xyz(local_adress=local_adress,
			  n_neighbors=n_neighbors, saveat=saveat)

	  print np.sort(concern_locals)
	  print len(concern_locals)

	return ranking_results

def process_localname(name):
	return_name = name.replace(".cif", "").replace("local_str_map/", "").replace("concern_region/", "")

	regs = ["11", "12", "13", "21", "22", "23", "31", "32", "33"]
	for reg in regs:
		return_name = return_name.replace("region{0}/".format(reg), "")

	return return_name

def get_job(address):
	for job in alljobs:
		if job in address:
			return job

def get_idx(basename, job, list_indexes):
	for ref in list_indexes:
		if job in ref and basename in ref:
			return ref
	return None


def local_map_stats(struct_dir_file, save_localdir, saveat):
	with open(struct_dir_file) as file:
		struct_reprs = file.read().splitlines()
		print struct_reprs

	full_regions = glob.glob(save_localdir + "/*")

	regions = [get_basename(k) for k in full_regions]

	# data = dict({k: dict() for k in regions})
	# data = dict()
	df = pd.DataFrame(columns=regions)

	target_variable = "energy_substance_pa"
	df[target_variable] = "non_converged"
	for region, full_region in zip(regions, full_regions):
		listdirs = glob.glob(full_region + "/*")
		# listdirs = [get_basename(k) for k in listdirs]

		listdirs = [process_localname(k) for k in listdirs]

		for local_adress in listdirs:
			tmp = [m.start() for m in re.finditer('_', local_adress)]

			# # last digit is local index
			# local_idx = local_adress[max(tmp) + 1:]
			# struct_repr = local_adress[: max(tmp)]

			# # last digit is local index
			local_idx = local_adress[local_adress.find(prefix_local_idx) + len(prefix_local_idx):]
			# # the rest is structure representation object
			struct_repr = local_adress[: local_adress.find(prefix_local_idx)]

			# # get dft_stab
			struct_basename = get_basename(struct_repr)
			struct_basename = struct_basename[:struct_basename.find(".")]

			struct_job = get_job(address=local_adress)
			
			try:
				tmp_df = pd.read_csv("{0}/result/standard/{1}.csv".format(maindir, struct_job), index_col=0)
				idx = "{0}/result/standard/{1}/{2}_cnvg_cnvg".format(maindir, struct_job, struct_basename)
				
				idx = get_idx(basename=struct_basename, 
					job=struct_job, list_indexes=tmp_df.index)
			except Exception as e:
				tmp_df = pd.DataFrame()
				idx = None
				print "Cannot read standard result file:", 
				pass      
			
			if struct_repr not in df.index:
				df.loc[struct_repr, :] = 0

			df.loc[struct_repr, region] += 1
			if idx is not None: 
				dft_stab = tmp_df.loc[idx, target_variable] # 
				df.loc[struct_repr, target_variable] = dft_stab
			else:
				df.loc[struct_repr, target_variable] = "non_cnvg"     

			# if struct_repr not in data.keys():
			#     data[struct_repr] = dict({region: 0})

			# if region not in data[struct_repr]:
			#     data[struct_repr][region] = 0

			# data[struct_repr][region] += 1

			# if struct_repr in data[region].keys():
			# else:

	# df = pd.DataFrame.dict(data, orient="index") # columns
	df.to_csv(saveat)
	print df

def df2X(filename, concern_center, pv):

  df = pd.read_csv(filename, index_col=0)
  if concern_center != "all":
	  df = df[df["center"] == concern_center]

  df.drop("center", 1, inplace=True)

  print df.columns
  if pv is None:
	pv = df.columns
  df = df[pv]
  X = df.values # values

  scaler = MinMaxScaler()
  X = scaler.fit_transform(X)

  return X, df.index

def tsne_search_params(config_tsne, job, concern_center, 
		annot_name, color_array, full_local_idx,
		marker_array):
  full_results = []
  pps = np.arange(1, 1000, 200) # 10, 50, 10
  ees = np.arange(1, 1000, 200)

  for pp in pps: # 20, 50, 80, 100, 200, 400, 800, 1500, 2000
	for ee in ees: # , 20, 40, 50, 100, 200
	  config_tsne["perplexity"] = pp
	  config_tsne["early_exaggeration"] = ee
		  
	  savefig = "{0}/result/local_str_map/{1}/tsne/{2}/{1}_{2}".format(maindir, job, concern_center)
	  savefig += "_pp{0}_ee{1}.pdf".format(pp, ee)


	  rst = dict({})

	  title = "locals from {}".format(get_basename(struct_dir_file))
	  tmp_rst = processing_map(X=X, title=title, saveat=savefig,
		full_local_idx=full_local_idx,
		  config=config_tsne,
		  annot_name=annot_name,
		  color_array=color_array,
		  marker_array=marker_array,
		  save_localdir=None, Xp=None)
	  rst["perplexity"] = pp
	  rst["early_exaggeration"] = ee
	  for k, v in tmp_rst.items():
		rst[k] = np.array(v)

	  print "pp, ee:", pp, ee
	  full_results.append(rst)
	  rst_df = pd.DataFrame(full_results)
	  rst_df.to_csv("{0}/result/local_str_map/{1}/{2}_ranking.csv".format(maindir, 
				job, concern_center))
		
if __name__ == '__main__':
	# # FOR INITIAL SUBSTITUTION STRUCTURES
	maindir = "/media/nguyen/work/SmFe12_screening"

	NdFeB_stab_ids = [
	"B10Ce2Ni1_id_2025086", "B6Ce1Cr2_id_94775", "B14Ho4Ni1_id_2109156", "B2Co2Dy1_id_2014474", 
	  "B4Ce1Ni1_id_2023368", "B4Ce1Cr1_id_2023392", "B5Eu2Os3_id_180411", "B4Eu2Rh5_id_183842", 
	  "B1Ce1Co4_id_2014377", "B6Ce1Ni12_id_2078436", "B4Co21Nd5_id_126928", "B6Co19Nd5_id_125302", 
	  "B1Nd1Ni4_id_2070064", "B2Ce3Ni13_id_1778822", "B4Ce3Co11_id_2015412", "B6Ce2Re3_id_1987748", 
	  "B4Ce1Ru4_id_2074934", "B2Ce1Ir2_id_180315", "B6Eu3Rh8_id_1771853", "B3Ce2Co7_id_2016628"
		  ]

	oqmd_NdFeB_stab = ["Nd2Fe14B", "B6Fe12Nd1_id_1925663id", "B6Fe2Nd5_id_613227id", 
		"B4Fe4Nd1_id_1775947id"]  
	oqmd_NdFeB_unstab = ["B1Fe2Nd1_id_1561590id", "B1Fe1Nd2_id_1611342id", 
	  "B6Fe3Nd4_id_672682id", "B2Fe1Nd1_id_1327316id"]
	jobs = [
			# "all_Mo", 
	
			# "Sm2-Fe23-Si1", 
			# "Sm2-Fe22-Si2", 
			# "Sm2-Fe21-Ga3",
			# "Sm2-Fe23-Ga1", 
			# "Sm2-Fe22-Ga2",
			# "Sm2-Fe21-M3_std6000"
			# "all_Ga"

			# "Sm2-Fe22-Ga2|with_org",
			# "Sm2-Fe21-Ga3|with_org"

			# "Sm-Fe11-M_wyckoff_2", 
			# "Sm-Fe10-M2_wyckoff",
			# "Sm2-Fe21-M3",
			# "Sm2-Fe22-M2|with_org", 

			# "Sm2-Fe23-M", 
			# "oqmd_Sm-Fe-Mo",
			# "check_Mo_2-22-2",
			# "CuAlZnTi_Sm2-Fe22-M2", "CuAlZnTi_Sm2-Fe21-M3"

			# "Nd-Fe-B",
			# "Nd-Fe-B|stable"
			"oqmd_NdFeB+Nd-Fe-B" # # for the NdFeB local structure paper
			]

	# # for the NdFeB local structure paper
	pv = [ #'s2-ofcenter', 
		'p1-ofcenter', 'd6-ofcenter', 'f4-ofcenter', 's2-s2', 
		'd6-s2', 'f4-s2', 's2-d6', 'd6-d6', 'd6-f4',
		'p1-d6', 'd6-p1','f4-d6', 's2-f4', 'f4-f4', 'p1-f4', 'f4-p1',
		's2-p1', 'p1-p1',  'p1-s2',
			] 

	# # for the Sm2-Fe21-M3_std6000 
	# pv = ['s1-ofcenter', 's2-ofcenter', 'd5-ofcenter', 'd6-ofcenter',
	#    'f6-ofcenter', 's1-s1', 's2-s1', 'd5-s1', 'd6-s1', 'f6-s1',
	#    's1-s2', 's2-s2', 'd5-s2', 'd6-s2', 'f6-s2', 's1-d5', 's2-d5',
	#    'd5-d5', 'd6-d5', 'f6-d5', 's1-d6', 's2-d6', 'd5-d6', 'd6-d6',
	#    'f6-d6', 's1-f6', 's2-f6', 'd5-f6', 'd6-f6']

	# # Sm2-Fe21-Ga3
	# pv = ['s2-ofcenter', 'p1-ofcenter', 'd6-ofcenter', 'd10-ofcenter',
	#    'f6-ofcenter', 's2-s2', 'p1-s2', 'd6-s2', 'd10-s2', 'f6-s2',
	#    's2-p1', 'p1-p1', 'd6-p1', 'd10-p1', 'f6-p1', 's2-d6', 'p1-d6',
	#    'd6-d6', 'd10-d6', 'f6-d6', 's2-d10', 'p1-d10', 'd6-d10',
	#    'd10-d10', 'f6-d10', 's2-f6', 'p1-f6', 'd6-f6', 'd10-f6']

	# # take all variable
	pv= None
	
	feature_type = "ofm1_no_d"
	# # FOR isolated job: e.g. 1-11-1, 1-10-2, 2-23-1, 2-22-2, ofm1_no_d
	jobs = ["{0}|{1}".format(job, feature_type) for job in jobs]

	# # FOR PROJECTING FILE e.g. projecting file is oqmd
	# project_file = "{0}/result/local_str_map/Sm2-Fe23-Ga1|{1}.csv".format(maindir, feature_type)
	
	project_file = None

	# # FOR MIXED JOBS, NOT RECOMMENDED
	# jobs = ["oqmd_Sm-Fe-Mo|with|{0}|{1}".format(job, feature_type) for job in jobs]

	is_local2file = True
	is_processing_map = True
	is_search_params = False
	is_local_map_stats = False
	 # # counting number of local structures in 11, 12, etc 
	concern_center = "B" # all, Fe, Nd, B, Sm
	print "concern_center", concern_center

	for job in jobs:
	  struct_dir_file = "{0}/input/rep_dirs/{1}.txt".format(maindir, job)
	  local_struct_file = "{0}/result/local_str_map/{1}.csv".format(maindir, job)
	  if is_local2file:
		  # # save local vector representation to file
		  get_local_vector(struct_dir_file=struct_dir_file, 
			  save_csv=local_struct_file, n_neighbors=5)
	  save_localdir = "{0}/result/local_str_map/{1}/concern_region_{2}".format(maindir, job, concern_center)
	  if is_processing_map:
		  # # create map of local structures
		  savefig = "{0}/result/local_str_map/{1}/{2}_{3}.pdf".format(maindir, job, job, concern_center)
		  # # X
		  X, full_local_idx = df2X(filename=local_struct_file, 
			  concern_center=concern_center, pv=pv)

		  # # X project
		  if project_file is not None:
			  dfp = pd.read_csv(project_file, index_col=0)
			  dfp = dfp[dfp["center"] == concern_center]
			  dfp.drop("center", 1, inplace=True)
			  Xp = dfp.values
		  else:
			  Xp = None

		  if os.path.exists(save_localdir) and os.path.isdir(save_localdir):
			  shutil.rmtree(save_localdir)
  
		  annot_name = np.array([get_basename(k) for k in full_local_idx])
		  # color_array = ["orange"] * X.shape[0]

		  # with open(struct_dir_file) as file:
		  #     struct_reprs = file.read().splitlines()
		  # count = 1
		  # with open("/media/nguyen/work/SmFe12_screening/input/rep_dirs/Nd-Fe-B|ofm1_no_d|stable.txt", "w") as fds:
		  #   for idx in struct_reprs:
		  #     for Nd_stab in NdFeB_stab_ids:
		  #       if Nd_stab in idx:
		  #         fds.write(str(idx) +"\n")
		  #         print count
		  #         count += 1

		  color_array = []
		  marker_array = []

		  for idx in full_local_idx:

			is_color = False
			for Nd_stab in NdFeB_stab_ids:
			  if Nd_stab in idx:
				  color_array.append("green") # red
				  marker_array.append("v")
				  is_color = True
				  break
			for oqmd_stab in oqmd_NdFeB_stab:
			  if oqmd_stab in idx:
				color_array.append("blue")
				marker_array.append("v")
				is_color = True
				break
			for oqmd_stab in oqmd_NdFeB_unstab:
			  if oqmd_stab in idx:
				color_array.append("blue")
				marker_array.append("o")
				is_color = True
				break

			if not is_color:
			  color_array.append("orange") # blue
			  marker_array.append("o")

		  # blue_idx = np.where((np.array(color_array)=="blue"))[0]
		  config_isomap = dict({"n_neighbors": 700, "n_components": 2, # # 50, 100, 500
				  "eigen_solver": "auto", "tol": 0, "max_iter": None, 
				  "path_method": "auto", "neighbors_algorithm": "auto", 
				  "n_jobs": None})
		  config_tsne = dict({"n_components":2, "perplexity":100.0,  # same meaning as n_neighbors
				  "early_exaggeration":50.0, # same meaning as n_cluster
				  "learning_rate":1000.0, "n_iter":1000,
				   "n_iter_without_progress":300, "min_grad_norm":1e-07, "metric":'euclidean', "init":'random',
				   "verbose":0, "random_state":None, "method":'barnes_hut', "angle":0.5, "n_jobs":None})
		  # # for search params
		  full_results = []
		  if is_search_params:
			tsne_search_params(config_tsne=config_tsne, job=job, 
			  full_local_idx=full_local_idx,
			  concern_center=concern_center, 
			  marker_array=marker_array,
			  annot_name=annot_name, color_array=color_array)

		  # # for best params
		  if not is_search_params:
			# Nd: 451-451, Fe: 451-451, B:201-251, B stable: 20-30
			config_tsne["perplexity"] = 201 # in manuscript Nd: 400-200, Fe: 800-200, B:201-100, B stable: 20-30
			config_tsne["early_exaggeration"] = 100 #
			title = "locals from {}".format(get_basename(struct_dir_file))
			ranking_results = processing_map(X=X, title=title, saveat=savefig,
				full_local_idx=full_local_idx, config=config_tsne,
				annot_name=annot_name, color_array=color_array,
				marker_array=marker_array,
				save_localdir=save_localdir, Xp=Xp)
			fig = plt.figure(figsize=(12, 8))

			for key, value in ranking_results.items():
			  plt.plot(value,  marker="o", linestyle="-.", 
				# color=color,
				alpha=1.0, label=key, linewidth=1, markersize=2, mfc='none')

			saveat = "{0}/result/local_str_map/{1}/{2}_ranking.pdf".format(maindir, 
					job, concern_center)
			plt.savefig(saveat)

	  
	  if is_local_map_stats:
		  # # create map of local structures
		  saveat = "{0}/result/local_str_map/{1}/region_stats.csv".format(maindir, job)
		  local_map_stats(struct_dir_file=struct_dir_file, save_localdir=save_localdir,
			  saveat=saveat)


	if False:
		# # for testing the get local_xyz only
		local_adress = "/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_4*j8_17*i8.ofm1_no_d_4"
		saveat = "/media/nguyen/work/SmFe12_screening/input/locals/check_Mo_2-22-2/ofm1_no_d/Mo_4*j8_17*i8.ofm1_no_d_4"
	   
		concern_locals = ['/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_9*j8_23*f8.ofm1_no_d_14',
				'/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_9*j8_23*f8.ofm1_no_d_16',
				'/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_9*j8_23*f8.ofm1_no_d_18',
				'/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_9*j8_23*f8.ofm1_no_d_19',
				'/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_9*j8_23*f8.ofm1_no_d_20',
				'/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_9*j8_23*f8.ofm1_no_d_22',
				'/media/nguyen/work/SmFe12_screening/input/feature/check_Mo_2-22-2/ofm1_no_d/Mo_9*j8_23*f8.ofm1_no_d_24']
		
		local_xyz_saves = [k.replace("/feature/", "/locals/") for k in concern_locals]
		
		for local_adress, saveat in zip(concern_locals, local_xyz_saves):
			get_local_xyz(local_adress=local_adress,
				n_neighbors=5, saveat=saveat)





