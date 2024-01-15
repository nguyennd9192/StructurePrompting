import os, sys, pickle
import numpy as np
import pandas as pd
from lib import *
from descriptor import OFMFeature, XRDFeature

class ProcessFeature(): 
	def __init__(self):
		self.export_features = []
		self.export_features_name = []
		self.is_rm_const_feature = False
	def create_export_feature(self, feature):
		list_features = [i for i in dir(feature) if not i.startswith('_')]

	def process_file_paths(self, structure_dir):
		struct_files = get_subdirs(sdir=structure_dir)
		data = []
		for struct_file in struct_files:
			if " " not in struct_file:
				with open(struct_file, "rb") as tmp_f:
					struct_repr_obj = pickle.load(tmp_f) #  encoding='latin1'  
					data.append(struct_repr_obj)

		self.data_path = data
		names = np.array([file.__getattribute__('name')
							 for file in self.data_path])
		saveats = np.array([file.__getattribute__('saveat')
							 for file in self.data_path])
		self.names = names
		self.saveats = saveats

		features = np.array([file.__getattribute__('feature')
							 for file in self.data_path])
		feature_names = np.array([file.__getattribute__('feature_name')
								 for file in self.data_path])

		if self.is_rm_const_feature:
			non_zeros = np.where(features.any(axis=0))[0]
			self.features = features[:, non_zeros]
			self.feature_names = feature_names[:, non_zeros][0]
		else:
			self.features = features
			self.feature_names = feature_names[0]

	def csv_export(self, more_info, csv_saveat, from_feature_names="all" ):
		# # #
		more_info_feature = np.array([[file.__getattribute__(f) for f in more_info]
							for file in self.data_path])
		# # take subset or take all self.features to csv
		if from_feature_names == "all":
			self.export_features_name = self.feature_names
		else:
			# # export only subset of original feature_names
			self.export_features_name = from_feature_names
		export_features_idxs = np.array([self.feature_names.tolist().index(i) for i in self.export_features_name])
		n_ft = len(export_features_idxs)
		
		
		# # except some structure cannot export feature
		exp_org_feature = []
		for ft in self.features:
			if ft is not None:
				tmp = ft[export_features_idxs]
			else:
				tmp = np.zeros(n_ft) 
			exp_org_feature.append(tmp)

		csv_data = np.concatenate([more_info_feature, exp_org_feature], axis=1)
		columns  = np.concatenate([more_info, self.export_features_name])
		# indexes = [get_basename(filename=k) for k in saveat]
		df = pd.DataFrame(csv_data, columns=columns)
		# df = df.loc[:, (df != 0).any(axis=0)]
		makedirs(csv_saveat)
		df.to_csv(csv_saveat)

	def get_local_structure(self, structure_dir, is_atom="all"):
		struct_files = get_subdirs(sdir=structure_dir)
		data = []
		# # load all struct_repr_obj
		for struct_file in struct_files:
			with open(struct_file, "rb") as tmp_f:
				struct_repr_obj = pickle.load(tmp_f) #  encoding='latin1'  
				data.append(struct_repr_obj)
		all_locals = np.array([file.__getattribute__('locals')
							 for file in data])
		all_atoms = np.array([file.__getattribute__('atoms')
							 for file in data])
		all_local_xyz = np.array([file.__getattribute__('local_xyz')
							 for file in data])

		local_idxes = []
		local_repr = []
		local_xyz = []

		if is_atom == "all":
			for atoms, struct_file in zip(all_atoms, struct_files):
				idxes = ["{0}/{1}".format(struct_file, a) for a in atoms]
				local_idxes.append(idxes)
			local_repr = copy.copy(all_locals)
		else:
			for atoms, this_locals, this_local_xyz, struct_file in zip(all_atoms, all_locals, all_local_xyz, struct_files):
				if atoms is not None:
					for a, c, c_x in zip(atoms, this_locals, this_local_xyz):
						idxes = []
						lcs = []
						lc_xyz = []
						if a == is_atom:
							idxes.append("{0}/{1}".format(struct_file, a))
							lcs.append(c)
							lc_xyz.append(c_x)
						if len(lcs) != 0:
							local_idxes.append(idxes)
							local_repr.append(np.array(lcs))
							local_xyz.append(np.array(lc_xyz))

		return local_repr, local_idxes, local_xyz

if __name__ == '__main__':
	for job in jobs:
		ft_dir = "{0}/{1}/ML/feature/{2}/".format(result_dir, 
					job, feature_type)
		csv_saveat = get_ofm_csv_dir(job=job)
		process = ProcessFeature()
		process.is_rm_const_feature = False
		process.process_file_paths(structure_dir=ft_dir)
		process.csv_export(more_info=["name", "saveat"], 
			csv_saveat=csv_saveat, from_feature_names="all")






