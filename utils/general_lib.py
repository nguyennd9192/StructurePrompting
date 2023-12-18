

import os, glob, ntpath, pickle, functools, copy, sys, random
from params import *

import pandas as pd
import numpy as np
from tensorflow.io import gfile

try:
	from utils import utils # # ignore in call create_data
except:
	pass
	

from embedding_space import EmbeddingSpace
from sklearn.preprocessing import MinMaxScaler


vmin_plt = dict({"fe":-0.15, "magmom_pa":1.2})
vmax_plt = dict({"fe":0.20, "magmom_pa":2.2})

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


def get_subdirs(sdir):
	subdirs = glob.glob(sdir+"/*")
	return subdirs

def get_basename(filename):
		head, tail = ntpath.split(filename)
		basename = os.path.splitext(tail)[0]
		return tail

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



def get_savedir(ith_trial):

	# # for running local, for visualize only
	s_dir = str(os.path.join(
		FLAGS.save_dir,
		"/".join([FLAGS.data_init, FLAGS.sampling_method,
					FLAGS.score_method, FLAGS.embedding_method, 
					FLAGS.mae_update_threshold, str(FLAGS.active_p),
					str(FLAGS.is_search_params)
					])))

	# # for running server, with different initial percentage
	# s_dir = str(os.path.join(
	# 		FLAGS.save_dir,
	# 		"/".join([FLAGS.data_init, str(FLAGS.init_pct), FLAGS.sampling_method,
	# 					FLAGS.score_method, FLAGS.embedding_method, 
	# 					FLAGS.mae_update_threshold, str(FLAGS.active_p),
	# 					str(FLAGS.is_search_params)
	# 					])))

	result_file = s_dir + "/trial_{}.pkl".format(ith_trial)
	unlbl_dir = result_file.replace(".pkl","/")

	if FLAGS.score_method == "u_gp_mt":
		mt_kernel = 1.0 # 1.0, 0.001
		fix_update_coeff = 1
		unlbl_dir += "_mt{}".format(mt_kernel)

	# # to mark whether update estimator by DQ only or DQ vs RND
	# # "" is update all
	# estimator_update_by = str(FLAGS.estimator_update_by).split('_') # , "RND", "OS"
	# if len(estimator_update_by) < 3:
	# 	for k in estimator_update_by:
	unlbl_dir += FLAGS.estimator_update_by

	return unlbl_dir


def get_qrindex(df, qid):
	dq = "query2update_DQ_{}".format(qid)
	os = "query_outstanding_{}".format(qid)
	rnd = "query_random_{}".format(qid)



	update_DQ_str = df.loc[df[dq]==dq, "unlbl_index"].to_list()
	outstand_str = df.loc[df[os]==os, "unlbl_index"].to_list()
	random_str = df.loc[df[rnd]==rnd, "unlbl_index"].to_list()
	return update_DQ_str, outstand_str, random_str

def get_color_feature(v):
	# # Nguyen

	feature_dict = dict({
			"of":"black",
			"s1":"darkblue", "s2":"royalblue",
			"p1":"moccasin", "f6":"orange",
			
			"d2":"azure", "d5":"turquoise", "d6":"powderblue", 
			"d7":"teal", "d10":"darkslategray",
			  
			})
	c = "black"
	for term in feature_dict.keys():
		tmp = term+"-"
		if tmp in v:
			c = feature_dict[term]

	return dict({c:"full"})

def get_color_112(index):
	# c = "black"
	colors = dict()
	ratios = []

	state_subs = 0
	# cdicts = dict({"Ga":"purple", "Mo":"red", "Zn":"orange", 
	# 	"Co":"brown", "Cu":"blue", "Ti":"cyan", "Al":"green"})

	# cdicts = dict({"Ga":"#C9EA3D", "Mo":"#DDD1E6", "Zn":"#642455", 
	# 	"Co":"#FFCB56", "Cu":"#1389A5", "Ti":"#9CDEEE", "Al":"#83A340"})

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

	#normalize item number values to colormap
	# norm = matplotlib.colors.Normalize(vmin=0, vmax=1000)

	#colormap possible values = viridis, jet, spectral
	# rgba_color = cm.jet(norm(400),bytes=True) 
	return colors

def get_marker_112(index):
	m = "|"
	if "1-11-1" in index:
		m = "P"
	elif "1-10-2" in index:
		m = "+" 
	elif "1-9-3" in index:
		m = "x" # 
	return m

def get_family(index):
	if "mix" in index:
		if "Sm-Fe9" in index:
			f = "1-9-3"
		elif "Sm-Fe10" in index:
			f = "1-10-2"
		elif "Sm-Fe11" in index:
			f = "1-11-1"
		elif "Sm2-Fe23" in index:
			f = "2-23-1"
		elif "Sm2-Fe22" in index:
			f = "2-22-2"
		elif "Sm2-Fe21" in index:
			f = "2-21-3"
		elif "2-22-2" in index:
			f = "2-22-2"
	else:
		family_match = dict({1:"1-11-1", 2:"1-10-2", 3:"1-9-3"})
		n_subs = str(index).count("__")
		f = family_match[n_subs]
	return f

def get_ratio(index, element):
	# # e.g. mix-_-Sm-Fe10-Al1-Ga1-_-Ga_9___Al_5
	start = index.find("mix-_-") + len("mix-_-")
	end = index.rfind("-_-")
	short_index = index[start:end]
	pos = short_index.find(element)
	r = int(short_index[pos+2:pos+3])

	return r

def	get_scatter_config(unlbl_index, index_train, selected_inds):
	mix_index = ["mix__"+k for k in unlbl_index]
	plot_index = np.concatenate((mix_index, index_train), axis=0)
	family = [get_family(k) for k in plot_index]

	list_cdict = np.array([get_color_112(k) for k in plot_index])
	marker_array = np.array([get_marker_112(k) for k in family])
	alphas = np.array([0.3] * len(plot_index))
	alphas[selected_inds] = 1.0 
	alphas[len(unlbl_index):] = 1.0
	return list_cdict, marker_array, alphas


def est_alpha_updated(X_train, y_train, 
				X_test, y_test, selected_inds, 
				estimator):
	# # 1. filter 
	# # hold or dimiss
	model = "empty"
	if selected_inds != []:
		# # update X_train, y_train by selected_inds
		if selected_inds is not None:
			X_train = np.concatenate((X_train, X_test[selected_inds]), axis=0)
			y_train = np.concatenate((y_train, y_test[selected_inds]), axis=0)
			assert X_test[selected_inds].all() != None

		# if FLAGS.mae_update_threshold != "update_all":
		# 	estimator.fit(X_train, y_train)
		# 	errors = estimator.predict(X_train)
		# 	for ith, err in zip(errors):
		# 		if err > FLAGS.mae_update_threshold:
		# 			X_train.remove(ith)
		# 			y_train.remove(ith)


	# # normalize
	scaler = MinMaxScaler().fit(X_train)
	X_train = scaler.transform(X_train)
	X_test = scaler.transform(X_test)

	# # transform to embedding or not
	if FLAGS.embedding_method != "org_space":
		model = EmbeddingSpace(embedding_method=FLAGS.embedding_method)
		if FLAGS.embedding_method == "LMNN":
			_y_train = np.round(y_train, 1)#.reshape(len(y_train), 1) #.astype(str)
			print ("_y_train", _y_train.shape)
			# print ("_y_train", len(set(_y_train)))

		else:
			_y_train = copy.copy(y_train)

		model.fit(X_train=X_train, y_train=_y_train)

		X_train = model.transform(X_val=X_train, get_min_dist=False)
		X_test = model.transform(X_val=X_test, get_min_dist=False)

	return X_train, y_train, X_test, model


def load_data(init_train_percent=None):
	# # read_load train data
	csvname = os.path.join(FLAGS.data_dir, FLAGS.data_init + ".csv")
	df = pd.read_csv(csvname, index_col=0)
	pv = list(df.columns)
	pv.remove(FLAGS.tv)


	X_train, y_train, index_train = utils.get_mldata(FLAGS.data_dir, FLAGS.data_init)
	X_test, y_test, index_test = utils.get_mldata(FLAGS.data_dir, FLAGS.data_target)

	if init_train_percent is not None:
		X_all = np.concatenate((X_train, X_test), axis=0)
		y_all = np.concatenate((y_train, y_test), axis=0)
		index_all = np.concatenate((index_train, index_test), axis=0)


		total_data = len(index_train) + len(index_test)
		n_train = int(init_train_percent * total_data)
		
		pos_idx = list(range(total_data))
		random.shuffle(pos_idx)

		pos_idx_train = pos_idx[:n_train]
		pos_idx_test = pos_idx[n_train:]

		X_train = X_all[pos_idx_train]
		y_train = y_all[pos_idx_train]
		index_train = index_all[pos_idx_train]

		X_test = X_all[pos_idx_test]
		y_test = y_all[pos_idx_test]
		index_test = index_all[pos_idx_test]

	return X_train, y_train, index_train, X_test, y_test, index_test, pv

def norm_id(id_qr):
	rmvs = [
			"/glusterfs/damlabfs/vasp_data/",
			"/Volumes/Nguyen_6TB/work/SmFe12_screening/input/feature/",
			"/coarse_relax", "/fine_relax", "/standard",
			".ofm1_no_d", "ofm1_no_d/", 
			# mix/Sm-Fe9-Ti1-Mo2/Mo_2-9___Ti_7
		# "6103324a5e1c4c889b5b74cbff538098/single_193_162jobs/",
		# "0fb43913fa7540e9aafd8d49465c7123/single_1102_77jobs/",
		# "0cfb92a839764d70a894766286cefcd2/single_1111_21jobs/",
		# "0f0ce41cc6a146a0a3c6595725634176/supp_0/",
		# "ce7eeeea21fb4b8784344b0ed536f8dd/supp_1/",
		# "bc336d85dca14791a99a35e1cac3ae59/supp_2/",
		# "d6719d4ecdb0405b81d7ad7ff554b693/supp_3/",
		# "9c40cfc311e2430e84718bf93f651cd4/supp_4/",
		# "0277b5462e8a4efc984b17562e80b226/supp_5/",
		# "6daead0ebb7044cf9cae51ced4267854/supp_6/",
		# "c11dc70398594be48aa0b71d335e0c1a/supp_7/",
		# "944b6d0b64cb45aeadda4dd71ab0a667/supp_8/",
		# # "437d235204e34f23953f71c12aa6724a/supp_9/",
		# # "53f9ec52625044a28bd51b0e17444cd8/supp_10/",
		# "b2d63e38b5a84a6199fb2ffaa2c63733/supp_11/"

		# # /home/nguyen/vasp_data/standard/mix/supp_2/mix-_-Sm-Fe9-Al2-Co1-_-Co_7___Al_9-11
		# # 
		"/home/nguyen/vasp_data",
		"single/single_193_162jobs/",	"single/single_1102_77jobs/",	"single/single_1111_21jobs/",
		"/mix/supp_0/",	"mix/supp_1/",	"/mix/supp_2/",	"/mix/supp_3/",
		"/mix/supp_4/",	"mix/supp_5/",	"/mix/supp_6/",	"/mix/supp_7/",
		"/mix/supp_8/",	"mix/supp_11/"
			]
	# # single/Sm-Fe9-Co3/ofm1_no_d/Co_2__j8_4__i8_9__i8.ofm1_no_d
	# # mix/Sm-Fe10-Al1-Ga1/ofm1_no_d/Ga_1___Al_8.ofm1_no_d
	is_query_from_feature = False
	if ".ofm1_no_d" in id_qr:
		is_query_from_feature = True

	for rm in rmvs:
		id_qr = id_qr.replace(rm, "")

	if "mix" in id_qr and is_query_from_feature:
		id_qr = id_qr.replace("/", "-_-")

	id_qr = get_basename(id_qr)


	# dot = "-_-"
	# if dot in id_qr:
	# 	last_point = id_qr.rfind(dot) #+ len("mix-_-")
	# 	id_qr = id_qr[last_point:]
	# 	print ("last_point", last_point)

	return id_qr
 
def id_qr_to_database(id_qr, std_results, 
		coarse_results=None, fine_results=None, tv="energy_substance_pa"):
	# # Still use in create_data
	# id_qr = arg[0]
	id_qr_cvt = norm_id(id_qr)
	if id_qr_cvt in std_results.index:
		target_y = std_results.loc[id_qr_cvt, tv]
	elif id_qr_cvt in fine_results.index:
		target_y = fine_results.loc[id_qr_cvt, tv]
		# print ("Add fine_relax results", target_y)
	elif id_qr_cvt in coarse_results.index:
		target_y = coarse_results.loc[id_qr_cvt, tv]
		# print ("Add coarse_relax results", target_y)
	else:
		target_y = None

	if "mix" in id_qr and type(target_y) is not float and target_y is not None:
		print ("=======")
		print ("old:", id_qr)
		print ("id_qr:", id_qr_cvt, target_y)
		# print ("None index:", id_qr_cvt, len(std_results.index))
	return id_qr, id_qr_cvt, target_y


def id_qr(qr, y, index):
	if qr in index:
		idx = np.where(index==qr)[0][0]
		value = y[idx]
	else:
		value = None
	return (qr, value)



def get_queried_data(qids, queried_files, unlbl_X, unlbl_y, unlbl_index,
			 embedding_model):
	"""
	database_results: *.csv of all vasp calculated data, normally in the standard step
	queried_files: all queried files
	unlbl_X: original ublbl data, before any querying
	unlbl_index: original index of ublbl data, before any querying
	"""

	dqs = []
	oss = []
	rnds = []

	for qid, qf in zip(qids, queried_files):
		this_df = pd.read_csv(qf, index_col=0)
		dq, os, rnd = get_qrindex(df=this_df, qid=qid)

		dq_cvt = map(functools.partial(id_qr, y=unlbl_y, index=unlbl_index), dq)
		os_cvt = map(functools.partial(id_qr, y=unlbl_y, index=unlbl_index), os)
		rnd_cvt = map(functools.partial(id_qr, y=unlbl_y, index=unlbl_index), rnd)

		dqs.extend(dq_cvt)
		oss.extend(os_cvt)
		rnds.extend(rnd_cvt)

	DQs, OSs, RNDs = np.array(dqs), np.array(oss), np.array(rnds)

	valid_Xyid = []

	for data in [DQs, OSs, RNDs]:
		_y =  data[:, -1] 
		valid_id = [i for i, val in enumerate(_y) if val != None] 
		_idx = data[valid_id, 0]
		valid_Xyid.append(_idx)
	
	return valid_Xyid


def get_all_unlb_y(database_results, unlbl_X, unlbl_index,
			coarse_db_rst, fine_db_rst):
	"""
	database_results: *.csv of all vasp calculated data, normally in the standard step
	queried_files: all queried files
	unlbl_X: original ublbl data, before any querying
	unlbl_index: original index of ublbl data, before any querying
	"""

	dqs = []

	frames = [pd.read_csv(k, index_col=0) for k in database_results]
	db_results = pd.concat(frames)
	index_reduce = [get_basename(k) for k in db_results.index]
	db_results["index_reduce"] = index_reduce
	db_results.set_index('index_reduce', inplace=True)

	# # coarse, fine db
	crs_frames = [pd.read_csv(k, index_col=0) for k in coarse_db_rst]
	crs_db_results = pd.concat(crs_frames)
	crs_db_results = crs_db_results.dropna()
	index_reduce = [get_basename(k) for k in crs_db_results.index]
	crs_db_results["index_reduce"] = index_reduce
	crs_db_results.set_index('index_reduce', inplace=True)

	fine_frames = [pd.read_csv(k, index_col=0) for k in fine_db_rst]
	fine_db_results = pd.concat(fine_frames)
	fine_db_results = fine_db_results.dropna()
	index_reduce = [get_basename(k) for k in fine_db_results.index]
	fine_db_results["index_reduce"] = index_reduce
	fine_db_results.set_index('index_reduce', inplace=True)


	data_cvt = map(functools.partial(id_qr_to_database, db_results=db_results,
		crs_db_results=crs_db_results, fine_db_results=fine_db_results), unlbl_index)
	data_cvt = np.array(list(data_cvt))

	valid_Xyid = []

	_y =  data_cvt[:, -1] 
	# # last column as predicted variable
	# # data in the format of: [id_qr, id_qr_cvt, target]
	valid_id = [i for i, val in enumerate(_y) if val is not None and val is not np.nan] 
	_y = np.array([float(k) for k in _y[valid_id]]) 

	# # _idx is index in queried file, unlbl_X
	_idx = np.array(data_cvt[valid_id, 0])

	_X = np.array(unlbl_X[[np.where(unlbl_index==k)[0][0] for k in _idx], :])
	valid_Xyid = (_X, _y, _idx)

	return valid_Xyid


def vasp_lbl2mix(unlbl_file, database_results, coarse_db_rst, fine_db_rst):
	"""
	database_results: *.csv of all vasp calculated data, normally in the standard step
	queried_files: all queried files
	unlbl_X: original ublbl data, before any querying
	unlbl_index: original index of ublbl data, before any querying
	"""
	unlbl_df = pd.read_csv(unlbl_file, index_col=0)
	unlbl_index = unlbl_df.index.to_list()

	frames = [pd.read_csv(k, index_col=0) for k in database_results]
	db_results = pd.concat(frames)
	index_reduce = [get_basename(k) for k in db_results.index]
	db_results["index_reduce"] = index_reduce
	db_results.set_index('index_reduce', inplace=True)

	# # coarse, fine db
	crs_frames = [pd.read_csv(k, index_col=0) for k in coarse_db_rst]
	crs_db_results = pd.concat(crs_frames)
	crs_db_results = crs_db_results.dropna()
	index_reduce = [get_basename(k) for k in crs_db_results.index]
	crs_db_results["index_reduce"] = index_reduce
	crs_db_results.set_index('index_reduce', inplace=True)

	fine_frames = [pd.read_csv(k, index_col=0) for k in fine_db_rst]
	fine_db_results = pd.concat(fine_frames)
	fine_db_results = fine_db_results.dropna()
	index_reduce = [get_basename(k) for k in fine_db_results.index]
	fine_db_results["index_reduce"] = index_reduce
	fine_db_results.set_index('index_reduce', inplace=True)

	y_index_cvt = map(functools.partial(id_qr_to_database, db_results=db_results,
		crs_db_results=crs_db_results, fine_db_results=fine_db_results), unlbl_index)
	y_index_cvt = np.array(list(y_index_cvt))
	unlbl_df["y_obs"] = None
	for a in y_index_cvt:
		id_qr, id_qr_cvt, target_y = a[0], a[1], a[2]
		if target_y != None:
			unlbl_df.loc[id_qr, "y_obs"] = target_y

	unlbl_df.to_csv(unlbl_file.replace(".csv","_with_lbl.csv",))
	



