import os, glob, ntpath, pickle, functools, copy, sys, re, gc, math, shutil, random
from matplotlib.colors import Normalize, LogNorm
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
import matplotlib as mpl

import pandas as pd 
import numpy as np
import datetime
try:
	from tensorflow.io import gfile
	from utils import utils # # ignore in call create_data
	from embedding_space import EmbeddingSpace
except:
	pass
from sklearn.preprocessing import MinMaxScaler
from ase.io import read, write
from pymatgen.core.structure import Structure
from graphviz import Digraph


tv = "energy_substance_pa"
numGenerations = 30
SmFe12_ene = 0.07

n_symm_seeding = 32 
max_symm_gen = 30 
min_symm_eps = 0.1 
max_symm_eps = 0.5 
n_symm_levels = 1 
max_one_attempt = 50

feature_type = "ofm1_with_d" # ofm1_with_d ofm1_no_d
estimator_method = "u_gp"
embedding_method = "MLKR"

symprec = 0.1
ith_trial = 1
eps_precision = 2 # 

main_dir = "/Users/nguyennguyenduong/Dropbox/My_code/USPEX_run/Structure Prompting" 
code_dir =  main_dir + "/code"
input_dir =  main_dir + "/input"
result_dir = copy.copy(main_dir)
jobs = [sys.argv[1]]

def get_basename(filename):
	head, tail = ntpath.split(filename)
	basename = os.path.splitext(tail)[0]
	return tail

# # for fix composition
pool_name = "skeleton/{0}".format(get_basename(jobs[0]))

	

is_inter = False
dt = datetime.datetime.today()
daytime = "-".join([str(dt.year), str(dt.month), str(dt.day)])

vmin_plt = dict({"energy_substance_pa":  0.05, "magmom_pa":1.5, "magmom_pv":0.0})
vmax_plt = dict({"energy_substance_pa":  0.2, "magmom_pa":2.0, "magmom_pv":1.2 })

subs_eles = ["Al", "Mo", "Zn", "Co", "Cu", "Ti", "Ga", "B", "Zr"] # # 
cdicts = dict({"Ga":"lightcoral", "Mo":"darkorange", "Zn":"yellow", 
		"Co":"lightskyblue", "Cu":"navy", "Ti":"maroon", "Al":"green",
		"B":"gray", "Zr":"pink"
		})


family_cdicts = dict({"Sm1Fe12":"#e30022",
					"Sm1Fe7":"#195905",
					"Sm1Fe5":"lime", # "#32cd32"  lime chartreuse
					"Sm3Fe9":"#fada5e",  
					"Sm2Fe4":"#bfff00",
					"Sm1Fe3":"blue",
					"Sm1Fe4":"pink"
					})

family_mdicts = dict({"Sm1Fe12":"o",
					"Sm1Fe7":"x",
					"Sm1Fe5":"+",
					"Sm3Fe9":"s",  
					"Sm2Fe4":"^"})

org_comp = ["Sm", "Fe"]

initialPopSize = 30
populationSize = 30 


toten_simple = {
	"Sm": -4.709, # http://oqmd.org/materials/entry/13648, R-3m, checked
	"Fe": -8.284, # http://oqmd.org/materials/entry/8568, Im-3m, checked
	# "Nd": -4.762, # http://oqmd.org/materials/entry/8137, P63/mmc
	"Mo": -10.848, # http://oqmd.org/materials/entry/676469, Im-3m
	"Zn": -1.240, # http://oqmd.org/materials/entry/9225, P63/mmc
	# "Co": -8.099, # http://oqmd.org/materials/entry/594128, Fm-3m

	"Co": -7.045, # http://oqmd.org/materials/entry/592136, Fm-3m
	"N": -8.207, # https://oqmd.org/analysis/calculation/37196, Pa-3
	"C": -9.193, #  https://oqmd.org/analysis/calculation/1258431, R-3m
	"B": -6.678 , # https://oqmd.org/materials/entry/620379. R-3m


	"Cu": -3.681, # http://oqmd.org/materials/entry/592441, Fm-3m
	"Ti": -7.775, # http://oqmd.org/materials/entry/9315, P6/mmm
	"Al": -3.744, # http://oqmd.org/materials/entry/8100, R-3m
	"Ga": -3.009, # http://oqmd.org/materials/entry/107684, Cmca
	# "V": -8.939, # http://oqmd.org/materials/entry/8082,  Im-3m
	"Zr": -8.546, # https://oqmd.org/materials/entry/8196, P63/mmc


	"P": -5.399, # https://oqmd.org/analysis/calculation/41443, P1
	"As": -4.652, # https://oqmd.org/analysis/calculation/1232946, R-3m

	}
init_mag_dict = dict({
	"Sm": 7.0, 
	"Fe": 5.0, 
	"Mo": 5.0, 
	"Zn": 0.0, 

	"Co": 5.0, 

	"Cu": 0.0, 
	"Ti": 5.0, 
	"Al": 0.0, 
	"Ga": 0.0, 
	"B":  0.0, 
	"C":  0.0, 	
	"Zr": 5.0, 
	"N":  0.0, 
	"P":  0.0, 
	"As":  0.0, 
	})

LANTHANIDE_MAG = {"La": 0.0, "Ce":2.143, "Pr": 3.200, "Nd":3.273,
 "Sm":0.741, "Eu":0.0, "Gd":-7.00, "Tb":-9.0, "Dy":-10.0, "Ho":-10.0,
 "Er":-9.0, "Tm":-7.0, "Yb":-4.0, "Lu":-0.0}
 

J4f = {"Ce":2.5, "Pr":4.0 , "Nd":4.5, "Pm":4,
	 "Sm":2.5, "Eu":0, "Gd":3.5, "Tb":6, "Dy":7.5, "Ho":8,
	 "Er":7.5, "Tm":6, "Yb":3.5}

gJ4f = {"Ce":6/float(7), "Pr":0.8 , "Nd":8/float(11) , "Pm":3/float(5),
	 "Sm":2/float(7), "Eu":0.0, "Gd":2.0, "Tb":3/float(2), "Dy":4/float(3), 
	 "Ho":5/float(4), "Er":6/float(5), "Tm":7/float(6), "Yb":8/float(7)}


eval_color = dict({"Al":"#83A340", "Mo":"#DDD1E6", "Zn":"#642455", "Co":"#FFCB56", 
				"Cu":"#1389A5", "Ti":"#9CDEEE", "Ga":"#C9EA3D"
	})


prop_lim = dict({
			"Volume": dict({"Sm1Fe12": [155, 180],
				"Sm1Fe11Co1": [155, 175], "Sm1Fe11Zr1": [170, 180], "Sm1Fe11Ti1": [165, 175], 
				"Sm1Fe11B1": [152, 170]}), 

			"magmom_pa": dict({"Sm1Fe12": [1.8, 2.2], "Sm1Fe11Zr1": [1.4, 1.8], 
				"Sm1Fe11Ti1": [1.4, 1.8], "Sm1Fe11Co1": [1.8, 2.2], "Sm1Fe11B1": [1.4, 1.8]}),

			"Density": dict({"Sm1Fe12": [7.5, 9.5], 
				"Sm3Fe9": [7.5, 9.5], 	"Sm1Fe5": [7.5, 9.5], "Sm2Fe4": [7.5, 9.5]})
		})


all_symm = list(range(1, 231, 1))
def release_mem(fig):
	fig.clf()
	plt.close()
	gc.collect()

def ax_setting():
	plt.style.use('default')
	plt.tick_params(axis='x', which='major', labelsize=13)
	plt.tick_params(axis='y', which='major', labelsize=13)

def makedirs(file):
	if not os.path.isdir(os.path.dirname(file)):
		os.makedirs(os.path.dirname(file))


def get_subdirs(sdir, tail=None):
	subdirs = glob.glob(sdir+"/*")
	if tail is not None:
		subdirs = [_ for _ in glob.glob(sdir+"/*") if  tail in _]
	return subdirs


def merge_two_dicts(x, y):
		"""Given two dictionaries, merge them into a new dict as a shallow copy."""
		z = x.copy()
		z.update(y)
		return z

def filter_array(x, y):
	# for i in y:
	# 	if type(i) != np.float64:
	# 		print (i, type(i))
	nan_x = np.isnan(x)
	nan_y = np.isnan(np.array(y).astype(np.float64))
	nan_id = nan_x + nan_y
	return x[~nan_id], y[~nan_id]


def load_pickle(filename):
	if not gfile.exists(filename):
		raise NameError("ERROR: the following data not available \n" + filename)
	data = pickle.load(gfile.GFile(filename, "rb"))
	return data

def load_flags(filename):
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument(filename, type=argparse.FileType('r'))
	args = parser.parse_args()

	return args


def dump_pickle(data, filename):
	pickle.dump(data, gfile.GFile(filename, 'wb'))

def dump_flags(data, filename):
	FLAGS.append_flags_into_file(filename)

def get_uspex_dir(job, uspex_file):
	idv = result_dir + "{0}/{1}".format(job, uspex_file)
	return idv

def get_micro_config(job):
	filename = result_dir + "{0}/config.yaml".format(job)
	return filename

def get_ofm_csv_dir(job):
	filename = result_dir + "{0}/ML/feature/ofm_csv/{0}___{1}.csv".format(job, feature_type)
	return filename

def get_ofm_ene_csv_dir(job):
	filename = result_dir + "{0}/ML/feature/ofm_ene_csv/{0}___{1}.csv".format(job, feature_type)
	return filename

def get_phase_diag_dir(job, ith_gen):
	filename = result_dir + "{0}/ML/Individuals_magmom_ene/Gen{1}/{0}___{2}.csv".format(job, ith_gen, job)
	return filename


def get_phase_corner(input_file):
	with open(input_file, 'r') as f:
		lines = f.readlines()
		for ith, line in enumerate(lines):
			if "numSpecies" in line:
				corner1 = np.array(lines[ith+1].split(),dtype=float)
				corner2 = np.array(lines[ith+2].split(),dtype=float)
				corner3 = np.array(lines[ith+3].split(),dtype=float)
				phase_corners = [corner1, corner2, corner3]
	return phase_corners

def get_input_block(input_file):
	with open(input_file, 'r') as f:
		lines = f.readlines()
		for ith, line in enumerate(lines):
			if "numSpecies" in line:
				block_id = 1
				is_all_blocks = False
				blocks = []
				while not is_all_blocks:
					if "EndNumSpecies" not in lines[ith+block_id]:
						block = np.array(lines[ith+block_id].split(),dtype=float)
						blocks.append(block)
						block_id += 1
					else:
						is_all_blocks = True

			if "atomType" in line:
				atomType = np.array(lines[ith+1].split(), dtype=str)

	return atomType, blocks


def Individuals2csv(filename):
	# # # # # # # # # # # # 
	# # convert Individuals uspex to csv
	# # # # # # # # # # # # 
	tocsv_file = filename + ".csv" 
	tof = open(tocsv_file, 'w')

	input_file = copy.copy(filename)
	input_file = input_file.replace("results1/Individuals", "INPUT.txt")

	phase_corners = []
	with open(input_file, 'r') as f:
		lines = f.readlines()
		for ith, line in enumerate(lines):
			if "atomType" in line:
				comps = lines[ith+1]
			if "numSpecies" in line and  "EndNumSpecies" in lines[ith+4]:
				corner1 = np.array(lines[ith+1].split(),dtype=float)
				corner2 = np.array(lines[ith+2].split(),dtype=float)
				corner3 = np.array(lines[ith+3].split(),dtype=float)
				phase_corners = [corner1, corner2, corner3]


	comps_array = comps.split()
	comps_list = []
	composition_order = ""
	for ith, comp in enumerate(comps_array):
		composition_order += "Composition{},".format(comp)
		comps_list.append("Composition{}".format(comp))

	with open(filename, 'r') as f:
		lines = f.readlines()
		for ith, line in enumerate(lines):
			line = line.replace("[","").replace("]", "")
			newline = ','.join(line.split())
			newline += '\n'
			if ith==0:
				newline = newline.replace("Composition,",composition_order)
				newline = newline.replace("KPOINTS,","KPOINTS1,KPOINTS2,KPOINTS3,")
			elif ith==1:
				continue
			tof.write(newline)
				# idx = "structure{}".format(ith -1)
	tof.close()
	return tocsv_file, comps_array, comps_list, phase_corners

def enthalpie2csv(filename):
	# # # # # # # # # # # # 
	# # convert enthalpies_complete.dat uspex to csv
	# # # # # # # # # # # # 
	tocsv_file = filename + ".csv" 
	tof = open(tocsv_file, 'w')

	lines = [i.strip().split() for i in open(filename).readlines()]
	# write it as a new CSV file
	with open(tocsv_file, 'wb') as f:
		for ith, line in enumerate(lines):
			newline = "{},".format(ith+1)
			newline += ','.join(line)
			newline += '\n'
			if ith==0:
				n_vr = len(line)
				headline = "ID,"
				headline += ','.join(['vr{}'.format(k+1) for k in range(n_vr)])
				headline += '\n'
				tof.write(headline)
			elif ith==1:
				continue
			tof.write(newline)
				# idx = "structure{}".format(ith -1)
	tof.close()
	return tocsv_file

def origin2csv(filename):
	# # # # # # # # # # # # 
	# # convert enthalpies_complete.dat uspex to csv
	# # # # # # # # # # # # 
	tocsv_file = filename + ".csv" 
	tof = open(tocsv_file, 'w')

	with open(filename, "r") as f:
		lines = f.readlines()

	# write it as a new CSV file
	with open(tocsv_file, 'wb') as f:
		for ith, line in enumerate(lines):
			parent_id = line[line.find("["):line.find("]")]
			parent_id = parent_id.replace("[", "").replace("]", "")

			parent_name = '+'.join(parent_id.split())
			print (parent_id, parent_name)
			line = line.replace(parent_id, parent_name)


			line = line.replace("[", "").replace("]", "")
			newline = ','.join(line.split())
			newline += '\n'
			tof.write(newline)
	tof.close()
	return tocsv_file


def get_magmom_v(filename):
	# f = open(filename,'r')
	# magmoms = f.readlines()[2:]
	# f.close()
	magmoms = np.loadtxt(filename, skiprows=2)
	return np.array(magmoms)


def get_SmFeX_family(X):
	tmp = copy.copy(org_comp)
	tmp.append(X)
	cmp_name = "".join(tmp)
	return tmp, cmp_name

def add_element(job):
	tmp = copy.copy(org_comp)
	job = job[job.find("/"):]
	for e in subs_eles:
		if e in job:
			tmp.append(e)
	return tmp


def color_code(index_all):
	list_cdict = np.array([dict({"black":"full"})] * len(index_all))
	marker_array = np.array(["o"] * len(index_all))
	alphas = np.array([0.1] * len(index_all))

	return list_cdict, marker_array, alphas

def getname_compound(df, idx, comps_array):
	comps = df.loc[idx, ["Composition{}".format(k) for k in comps_array]]
	str_comp = ""
	for ele, c in zip(comps_array, comps):
		if c!=0:
			str_comp += "{0}{1}".format(ele, int(c))
	return str_comp



def get_normcolor(c_value, v_range=None, islog=False, cmap=cm.RdYlGn):
	if v_range is None:
		vmin = float(min(c_value))
		vmax = float(max(c_value))
		# cmap = cm.Oranges # jet

	else:
		vmin, vmax = v_range[0], v_range[1]
		
	 # seismic, bwr

	if islog:
		normc = LogNorm(vmin=vmin, vmax=vmax)
	else:
		normc = Normalize(vmin=vmin, vmax=vmax)
	colors = []
	n_point = len(c_value)
	for i in range(n_point):
		value = normc(float(c_value[i]))
		colors.append(cmap(value))
	mappable = cm.ScalarMappable(norm=normc, cmap=cmap)
	return colors, mappable

def write_seed_file(write_to, seed_files):
	if os.path.isdir(write_to):
		os.remove(write_to)
	with open(write_to, 'w') as result:
		for ith, seed_file in enumerate(seed_files):
			# struct = qmpy.io.poscar.read(seed_file)
			# qmpy.io.poscar.write(struct, seed_file)
			struct = read(seed_file)
			write(filename=seed_file, images=struct)

			with open(seed_file, "r") as f:
				lines = f.readlines()

				for ith, line in enumerate(lines):
					line = line.replace("Cartesian", "Direct")
					tmp = line.split()
					
					if ith==1:
						if "1.0000000000000000" in line:
							newline = "1.0"
					elif ith in [2, 3, 4] or ith > 7:
						newline = "  {0}   {1}   {2}".format(tmp[0], tmp[1], tmp[2])
					else:
						newline =" ".join(tmp)

					newline = newline.replace(" -", "-")
					# minus_pos = [i for i, a in enumerate(newline) if a == "-"]
					# for p in minus_pos:
					# 	newline[]

					newline += "\n"
					result.write(newline)


def write_seed_file_for_subs(write_to, seed_files):
	if os.path.isdir(write_to):
		os.remove(write_to)
	with open(write_to, 'w') as result:
		for ith, seed_file in enumerate(seed_files):
			struct = read(seed_file)
			# write(filename=seed_file, images=struct)
			with open(seed_file, "r") as f:
				lines = f.readlines()
				for ith, line in enumerate(lines):
					line = line.replace("direct", "Direct")
					tmp = line.split()
					if ith in [2, 3, 4] or ith > 7:
						newline = "  {0}   {1}   {2}".format(tmp[0], tmp[1], tmp[2])
					else:
						newline =" ".join(tmp)

					newline = newline.replace(" -", "-")
					newline += "\n"
					result.write(newline)


def norm_uspex_id(full_path, rmvs=None):
	full_path = full_path.replace(result_dir, "")


	if rmvs is not None:
		for rm in rmvs:
			full_path = full_path.replace(rm, "")
	idx = full_path[full_path.find("EA")+len("EA"):full_path.find(".poscar")]

	# # to fix for SmBFeTi_00001/1020 10.142
	space = idx.find(" ") 
	print (idx, space)
	if space != -1:
		idx = idx[:space]

	return int(idx)


def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def get_angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


def copyfiles(org_fld, dest_fld, files, prefixes):
	for ith, file in enumerate(files):
		src = org_fld + "/" + file
		if prefixes is None:
			dst = dest_fld + "/" + file
		else:
			dst = dest_fld + "/{0}_{1}".format(prefixes[ith], file)	

		if os.path.exists(src):
			shutil.copyfile(src, dst)



def get_hidx(x, y, indexes, n_bin=10):
	
	hx = []
	hy = []
	hidx = []
	bins = np.linspace(np.nanmin(x), np.nanmax(x), int(n_bin))
	for ith in range(n_bin-1):
		idx = np.where((x>bins[ith]) & (x<bins[ith+1]))[0]
		if len(idx) != 0:
			xb = x[idx]
			yb = y[idx]
			idxb = indexes[idx]
			min_idx = np.argmin(yb)

			hx.append(xb[min_idx])
			hy.append(yb[min_idx])
			hidx.append(idxb[min_idx])
	return hx, hy, hidx


def get_color_marker(symm, bounds=None, is_periodic=False):

	if not is_periodic:
		cnorm = mpl.colors.Normalize(vmin=12, vmax=191)
		symm_marker = dict({ 
			191: "_",
			166: "^",
			139: "s",
			123: "2",
			115: "D", 
			
			65: "p",
			12: "+",



			# # # for ternary SmFeB
			160: "o", 
			146: "h", 

			})

		if symm in symm_marker.keys():
			m = symm_marker[symm]
		else:
			m = "." # 2
	else:
		# # # norm by % of symm
		cnorm = mpl.colors.Normalize(vmin=1, vmax=230)
		symm_marker = dict({0:"o", 1: "D", 2:"x", 3:"1", 4:"s", 5:"+", 
					6:7, 7:"^", 8:"p", 9:"2",  })
		m = symm_marker[int(symm % 10)]


	# rainbow, RdYlGn_r hsv tab20 YlOrRd, seismic
	cmap = plt.get_cmap('seismic')

	if bounds is not None:
		## Prepare bins for the normalizer
		norm_bins = np.sort(bounds) + 0.5
		norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

		# # https://colorkit.co/palette/c7522a-e5c185-fbf2c4-74a892-008585/
		discret_color = [
				# "lightsteelblue", "royalblue", "navy",  
				"#F0FFFF",   "#DCDCDC", "#696969", # # # very light white
				"#74a892", "#008585",  "#003f5c", # green
				 "#fbf2c4", "#e5c185",  "#c7522a", # red

				]


		cmap = mpl.colors.ListedColormap(discret_color)
		cnorm = mpl.colors.BoundaryNorm(norm_bins, cmap.N)
	return cnorm, cmap, m

def get_ratio(index, element):
	# # e.g. mix-_-Sm-Fe10-Al1-Ga1-_-Ga_9___Al_5
	start = index.find("mix-_-") + len("mix-_-")
	end = index.rfind("-_-")
	short_index = index[start:end]
	pos = short_index.find(element)
	r = int(short_index[pos+2:pos+3])

	return r


def get_color_112(index):
	# c = "black"
	colors = dict()
	ratios = []

	state_subs = 0

	cdicts = dict({"Ga":"lightcoral", "Mo":"darkorange", "Zn":"yellow", 
		"Co":"lightskyblue", "Cu":"navy", "Ti":"maroon", "Al":"green"})

	index = index.replace("CuAlZnTi_", "")
	if "mix" in index:
		for element, color in cdicts.items():
			if element in index:
				ratio = get_ratio(index=index, element=element)
				colors[color] = ratio
	else:
		for element, color in cdicts.items():
			if element in index and "CuAlZnTi" not in index:
				colors[color] = "full"
	return colors


def read_poscar(filename):
	with open(filename, 'r') as f:
		inp = f.read()
		struct = Structure.from_str(input_string=inp, fmt='poscar', sort=True,
			 merge_tol=0.1)
	return struct



def get_gen_struct(job, struct_id):
	sdir = get_uspex_dir(job=job, uspex_file="ML/symstruct/gatheredPOSCARS")
	poscar_file = "{0}/EA{1}.poscar".format(sdir, struct_id) 
	structure = read_poscar(filename=poscar_file)
	symm = structure.get_space_group_info(symprec=0.1, angle_tolerance=5.0)[1]

	return structure, symm

def normalize(x_train, x_test):
	scaler = MinMaxScaler().fit(x_train)
	x_train = scaler.transform(x_train)
	x_test = scaler.transform(x_test)
	return x_train, x_test

def get_common_idx(a, b):
	a = set(a)
	b = set(b)
	common_idx = list(a.intersection(b))
	# assert len(common_idx) < 0.9 * len(b)

	return common_idx

def plot_ene_vol_sym(df, v1, v2, ax, symm, label, s, org_struct_dir, save_hull_structure, job):
	facecolors = None
	cdicts = dict({"Sm1Fe12":"red", "Sm3Fe9":"cyan", 
		"Sm1Fe5":"green", "Sm2Fe4":"orange", "Sm1Fe7":"white"})
	for e in cdicts.keys():
		if e in job:
			facecolors = cdicts[e]

	if v1 in df.columns and v2 in df.columns:
		y = df[v1].values
		x = df[v2].values
		mag = df["magmom_pa"].values

		structure_indexes = df.index
		symm_label = []
		plt_label = []
		hidx = []
		if symm == "Seeds":
			seed_idx = df[df["Origin"]=="Seeds"].index
			hy = df.loc[seed_idx, v1].values 
			hx = df.loc[seed_idx, v2].values
			hmag = df.loc[seed_idx, "magmom_pa"].values

			hidx = copy.copy(seed_idx)
			for _, seed_structure in enumerate(seed_idx):
				structure, this_symm = get_gen_struct(job=job, struct_id=seed_structure)
				# if this_symm > 1:
				cnorm, cmap, m = get_color_marker(symm=this_symm)
				ax.scatter(hx[_], hy[_], s=s, alpha=0.8, marker=m, 
							c=this_symm, label="Seeds", norm=cnorm, 
							cmap=cmap, edgecolor="black"
							)
				if facecolors is not None:
					ax.scatter(hx[_], hy[_], s=10, marker="o", 
					facecolors=facecolors,
					edgecolor="white", linewidth=1
							)
				if this_symm not in symm_label and hy[_] < 0.1:
					plt_label.append([hx[_], hy[_], int(this_symm) ])
					symm_label.append(this_symm)

		# elif symm == 1:
		# 	flt_idx = np.where(y < 0.07)[0]
		# 	hx, hy, hidx = x[flt_idx], y[flt_idx], structure_indexes[flt_idx]
		# 	hmag = df.loc[hidx, "magmom_pa"].values

		# 	# ax.plot(hx, hy, linestyle='-', color="gray")
		# 	# ax.scatter(hx, hy, s=1, alpha=0.8, marker=".", c="black")
		# 	if facecolors is not None:
		# 		ax.scatter(hx, hy, s=10, marker="o", 
		# 		facecolors=facecolors)

		# 	if len(hy) != 0:
		# 		lbl_x, lbl_y = hx[np.argmin(hy)], hy[np.argmin(hy)]
		# 		if lbl_y < 0.1:
		# 			ax.text(lbl_x, lbl_y, "{} \n {}".format(int(symm), get_basename(job)))
		else:
			cnorm, cmap, m = get_color_marker(symm=symm)


			# hx, hy, hidx = get_hidx(x=x, y=y, indexes=structure_indexes, n_bin=20)
			# hmag = df.loc[hidx, "magmom_pa"].values

			hx = copy.copy(x)
			hy = copy.copy(y)
			hmag = copy.copy(mag)
			hidx = copy.copy(structure_indexes)
			
			c = [symm] * len(hx)
			ax.scatter(hx, hy, s=s, alpha=0.8, marker=m, 
						c=c, norm=cnorm, 
						cmap=cmap, label=label,
						edgecolor="black"
						# edgecolor=edgecolor, linewidth=3
						)
			if facecolors is not None:
				ax.scatter(hx, hy, s=10, marker="o", 
						facecolors=facecolors,
						# edgecolor="white", linewidth=1
						)
			# lbl_x, lbl_y = hx[np.argmin(hy)], hy[np.argmin(hy)]
			# ax.text(lbl_x, lbl_y, "{} \n {}".format(int(symm), get_basename(job)))

		for lbl in plt_label:
			ax.text(lbl[0], lbl[1], lbl[2])

		# # save hull structures to file
		if save_hull_structure is not None and len(hidx) != 0:
			files = ["EA{0}.poscar".format(_) for _ in hidx]
			prefixes = []
			min_ene_id = np.argmin(hy)
			for ith, (_m, _x, _y) in enumerate(zip(hmag, hx, hy)):
				# if symm == "Seeds":
				# 	text = "s_{0}_{1}_{2}_{3}".format(symm, round(_m, 3), round(_x, 3), round(_y, 3)) 
				# elif _y < 0.12 and symm > 10:
				# 	text = "+{0}_{1}_{2}_{3}".format(symm, round(_m, 3), round(_x, 3), round(_y, 3)) 
				# else:
				text = "{0}_{1}_{2}_{3}".format(symm, round(_m, 3), round(_x, 3), round(_y, 3))

				if ith == min_ene_id:
					text = "min_" + text
				prefixes.append(text)

			print (symm, files)
			makedirs(save_hull_structure+"/tmp.txt")
			copyfiles(org_fld=org_struct_dir, dest_fld=save_hull_structure, 
				files=files, prefixes=prefixes)
	return df.loc[hidx, :]








