from lib import *
import pymatgen as pm
from pymatgen.io.cif import CifWriter
import shutil, graphviz, pydot
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from acquisition import acquisition
from symmetry_generation import generate_subgroup_path
from descriptor import *
from process_features import ProcessFeature
from utils import get_model
from symmetry_transformation import symmetry_structure_generator

def get_train_test_data(train_file, test_file):
	df = pd.read_csv(train_file, index_col="ID")
	columns = list(df.columns)
	index_all = df.index
	pv = []
	for c in columns:
		if "-" in c and "-Type" not in c:
			pv.append(c)
	x_test, y_test, idx_test = None, None, None

	if test_file is not None:
		test_df = pd.read_csv(test_file, index_col="ID")
		test_columns = set(test_df.columns)
		idx_test = list(test_df.index)
		pv = list(test_columns.intersection(set(pv)))
		x_test = test_df[pv].values
		if tv in test_df.columns:
			y_test = test_df[tv].values


	x_train = df[pv].values
	y_train = df[tv].values
	scaler = MinMaxScaler().fit(x_train)
	x_train = scaler.transform(x_train)

	if test_file is not None:
		x_test = scaler.transform(x_test)

	return x_train, y_train, x_test, y_test, pv, scaler, idx_test

def get_predictor(estimator, x_train, y_train, x_test, y_test):
	if embedding_method=="MLKR":
		# # method as MLKR then transform x_train, x_test
		embedding = EmbeddingSpace(embedding_method="MLKR")
		embedding.fit(X_train=x_train, y_train=y_train)

		x_train_trans = embedding.transform(X_val=x_train, get_min_dist=False)
		x_test_trans = embedding.transform(X_val=x_test, get_min_dist=False)

		pca = PCA(n_components=2)
		x_train = pca.fit_transform(x_train_trans)
		x_test = pca.transform(x_test_trans)

	if estimator is None:
		estimator = get_model(
			method=estimator_method, seed=ith_trial, is_search_params=True, 
			n_shuffle=10000, mt_kernel=None)
		estimator.fit(x_train, y_train)
	y_test_pred = estimator.predict(x_test)
	variance = estimator.predict_proba(X=x_test, is_norm=True)


	return embedding, estimator, x_train, x_test, y_test_pred, variance




def symmetry_candidate_generator(gen, initial_candidates):
	eps = round(random.uniform(min_symm_eps, max_symm_eps), eps_precision)
	gen_dir = get_uspex_dir(job=job, uspex_file="ML/symm_gen/gen{0}/eps_{1}/".format(gen, eps))
	gen_ofm_dir = get_uspex_dir(job=job, uspex_file="ML/symm_gen_ft/gen{0}/eps_{1}/".format(gen, eps))

	gen_dir += str(pool_name)
	if not os.path.exists(gen_dir):
		for input_struct_file in initial_candidates:
			generate_subgroup_path(savedir=gen_dir+"/"+get_basename(input_struct_file), eps=eps, 
								is_primcell=True, input_struct_file=input_struct_file)

			symmetry_structure_generator(gen_dir=gen_dir+"/"+get_basename(input_struct_file), 
					group=None,	species=None, numIons=None,
					initial_candidate=input_struct_file, n_gen=1
					)

	candidates = [os.path.join(path, name) for path, subdirs, files in os.walk(gen_dir) for name in files]

	for cand_poscar in candidates:
		config = dict({"filename": cand_poscar, 
				"ft_type": feature_type,
				"is_ofm1":True, "is_including_d": "with_d" in feature_type })
		saveat = cand_poscar.replace("/symm_gen/", "/symm_gen_ft/").replace(".poscar", ".{0}".format(feature_type))
		if not os.path.isfile(saveat):
			c = OFMFeature(name=cand_poscar, saveat=saveat, config=config)
			feature = c.create_feature()
			feature_name = c.get_feature_name()
			c.saveto()
		else:
			print ("File existed.")
	ofm_dirs = [path.replace("/symm_gen/", "/symm_gen_ft/").replace(".poscar", ".{0}".format(feature_type)) 
				for path, subdirs, files in os.walk(gen_dir) if len(subdirs) == 0]

	for ofm_dir in ofm_dirs:
		csv_saveat = ofm_dir + ".csv"
		if not os.path.isfile(csv_saveat):
			process = ProcessFeature()
			process.is_rm_const_feature = False
			process.process_file_paths(structure_dir=ofm_dir)
			process.csv_export(more_info=["name", "saveat"], 
				csv_saveat=csv_saveat, from_feature_names="all")

	all_gen_dir = get_uspex_dir(job=job, uspex_file="ML/symm_gen_ft/gen{0}".format(gen))
	all_gen_files = [os.path.join(path, name) for path, subdirs, files in os.walk(all_gen_dir) for name in files if ".csv" in name]


	df_concat = pd.concat([pd.read_csv(f, index_col="name") for f in all_gen_files ], ignore_index=False)
	gen_ft_csv = all_gen_dir + ".csv"
	df_concat.rename(columns = {'name':'ID'}, inplace = True)
	df_concat.index.names = ['ID']
	df_concat.to_csv(gen_ft_csv)
	return gen_dir, gen_ofm_dir, gen_ft_csv


def symmetry_seeding(job, uspex_file, idv_file, save_fig_at):
	df = pd.read_csv(uspex_file, index_col="ID")
	idv_df = pd.read_csv(idv_file, index_col="ID")
	this_gen = np.nanmax(idv_df["Gen"])
	next_gen = int(this_gen+1)
	# # # create symmetry generator file

	for ith_gen in range(next_gen):
		lastgen_struct = get_uspex_dir(job=job, uspex_file="ML/symm_gen/gen{0}".format(ith_gen))
		lastgen_ft = get_uspex_dir(job=job, uspex_file="ML/symm_gen_ft/gen{0}".format(ith_gen))
		
		for d in [lastgen_struct, lastgen_ft]:
			if os.path.exists(d):
				shutil.rmtree(d) 

		lastgen_pkl = get_uspex_dir(job=job, uspex_file="ML/symm_gen/gen{0}.pkl".format(ith_gen))
		if os.path.exists(lastgen_pkl):
			os.remove(lastgen_pkl) 

	
	save_estimator = get_uspex_dir(job=job, uspex_file="ML/symm_gen/gen{0}.pkl".format(next_gen))
	save_seed_at = get_uspex_dir(job=job, uspex_file="Seeds/POSCARS_{0}".format(next_gen))

	selected_symm = []
	
	initial_candidates = []

	# # get last two gens
	init_gen_df = copy.copy(idv_df) # idv_df[idv_df["Gen"] >= this_gen-2]

	micro_prompt_ranking_file = get_uspex_dir(job=job, uspex_file="ML/micro-prompt/ranking.csv")
	if os.path.isfile(micro_prompt_ranking_file):
		mp_df = pd.read_csv(micro_prompt_ranking_file, index_col=0)
		mp_df_sorted = mp_df.sort_values(["acquisition_scores"], ascending=False).head(1)
		prompted_symm_set = mp_df_sorted["SYMM"].values
	else:
		prompted_symm_set = sorted(list(set(init_gen_df["SYMM"].values)), reverse=True)

	for n_symm, this_symm in enumerate(prompted_symm_set):
		# this_df = idv_df[idv_df["SYMM"] == this_symm]
		this_df = init_gen_df[init_gen_df["SYMM"] == this_symm]

		min_eneEA_index = this_df[tv].idxmin()
		this_min_ene = this_df.loc[min_eneEA_index, tv]

		comp_cols = [k for k in idv_df.columns if "Composition" in k]
		this_comps =  this_df.loc[min_eneEA_index, comp_cols].values

		if 0 in this_comps:
			continue
		max_magEA_index = None

		if "magmom_pa" in this_df.columns:
			max_magEA_index = this_df["magmom_pa"].idxmax()
		
		# if this_min_ene < 0.1:
		min_ene_EA = get_uspex_dir(job=job, 
			uspex_file="ML/symstruct/gatheredPOSCARS/EA{}.poscar".format(min_eneEA_index))
		max_mag_EA = get_uspex_dir(job=job, 
			uspex_file="ML/symstruct/gatheredPOSCARS/EA{}.poscar".format(max_magEA_index))

		for this_EA in [min_ene_EA]: # max_mag_EA
			if os.path.isfile(this_EA):
				initial_candidates.append(this_EA)
				selected_symm.append(this_symm)

		if len(initial_candidates) == 20:
		# if this_symm < 65: # n_symm > 20 or 
			break

	gen_dir, gen_ofm_dir, gen_ft_csv = symmetry_candidate_generator(
		gen=next_gen, initial_candidates=initial_candidates)
	
	one_batch = int(n_symm_seeding)

	x_train, y_train, x_test, y_test, pv, scaler, idx_test = get_train_test_data(
		train_file=uspex_file, test_file=gen_ft_csv)

	estimator = None
	if os.path.isfile(save_estimator):
		estimator = load_pickle(save_estimator)
	
	embedding, estimator, x_train_ebd, x_test_ebd, y_test_pred, variance = get_predictor(
		estimator, x_train, y_train, x_test, y_test)

	if not os.path.isfile(save_estimator):
		dump_pickle(data=estimator, filename=save_estimator)

	symm_test = []
	for idt in idx_test:
		if "_trial_" in idt:
			s = int(idt[idt.find(".poscar/") + 8: idt.find("_trial_")])
		else:
			idt = get_basename(idt)
			s = int(idt[:idt.find("_")])
		symm_test.append(s)

	symm_sort = np.argsort(symm_test)[::-1] 
	y_pred_sort = np.argsort(y_test_pred)
	variance_sort = np.argsort(variance)# [::-1]

	# # get 1/3 as min_ene + 1/3 as max_var +1/3 as random
	min_ene_pick_ids = [idx_test[_] for _ in y_pred_sort[:one_batch]]
	max_var_pick_ids = []
	max_symm_pick_ids = []
	random_id = []

	for i in variance_sort:
		this_id = idx_test[i]
		if this_id not in min_ene_pick_ids and len(max_var_pick_ids) < one_batch:
			max_var_pick_ids.append(this_id)

	for i in symm_sort:
		this_id = idx_test[i]
		if this_id not in min_ene_pick_ids and len(max_symm_pick_ids) < one_batch:
			max_symm_pick_ids.append(this_id)

	rest_id = list(set(idx_test) - set(min_ene_pick_ids) - set(max_var_pick_ids) - set(max_symm_pick_ids)) # 
	random_id = random.sample(rest_id, one_batch)



	test_df = pd.read_csv(gen_ft_csv, index_col="ID")
	test_df.loc[idx_test, "{}_pred".format(tv)] = y_test_pred
	test_df.loc[idx_test, "{}_var".format(tv)] = variance

	test_df["pickup"] = None
	test_df["pickupID"] = None


	test_df.loc[min_ene_pick_ids, "pickup"] = "selected_by_min_{}_pred".format(tv)
	# test_df.loc[max_var_pick_ids, "pickup"] = "selected_by_min_var_{}_pred".format(tv)
	test_df.loc[max_symm_pick_ids, "pickup"] = "selected_by_max_symm"
	test_df.loc[random_id, "pickup"] = "selected_by_random"


	seed_files = list(set(np.concatenate((min_ene_pick_ids, max_var_pick_ids, random_id, max_symm_pick_ids))))
	test_df.to_csv(gen_ft_csv.replace(".csv", "_pred.csv"))
	# save seed to next gen
	if os.path.isfile(save_seed_at):
		os.remove(save_seed_at)

	makedirs(save_seed_at)
	write_seed_file(write_to=save_seed_at, seed_files=seed_files)


if __name__ == '__main__':
	for job in jobs:
		ofm_ene_file = get_ofm_ene_csv_dir(job=job)
		idv_file = get_uspex_dir(job=job, uspex_file="results1/Individuals_magmom_ene.csv")
		save_fig_at = get_uspex_dir(job=job, uspex_file="ML/uspex_population")
		symmetry_seeding(job=job, uspex_file=ofm_ene_file, 
				idv_file=idv_file, save_fig_at=save_fig_at)



















