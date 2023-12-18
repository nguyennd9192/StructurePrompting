import re 
import pymatgen as pm
import numpy as np
from pymatgen.analysis.local_env import VoronoiNN

"""
python poscar2ofm.py number_cores is_ofm1 is_including_d dir_1 dir_2 ...
"""
 
def get_element_representation(name='Si'):

	"""
	generate one-hot representation for a element, e.g, si = [0.0, 1.0, 0.0, 0.0, ...]
	:param name: element symbol
	"""
	element = pm.Element(name)
	general_element_electronic = {'s1': 0.0, 's2': 0.0, 'p1': 0.0, 'p2': 0.0, 'p3': 0.0,
																'p4': 0.0, 'p5': 0.0, 'p6': 0.0,
																'd1': 0.0, 'd2': 0.0, 'd3': 0.0, 'd4': 0.0, 'd5': 0.0,
																'd6': 0.0, 'd7': 0.0, 'd8': 0.0, 'd9': 0.0, 'd10': 0.0,
																'f1': 0.0, 'f2': 0.0, 'f3': 0.0, 'f4': 0.0, 'f5': 0.0, 'f6': 0.0, 'f7': 0.0,
																'f8': 0.0, 'f9': 0.0, 'f10': 0.0, 'f11': 0.0, 'f12': 0.0, 'f13': 0.0, 'f14': 0.0}

	general_electron_subshells = ['s1', 's2', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6',
																'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10',
																'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7',
																'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'f14']

	if name == 'H':
			element_electronic_structure = ['s1']
	elif name == 'He':
			element_electronic_structure = ['s2']
	else:
			element_electronic_structure = [''.join(pair) for pair in re.findall("\.\d(\w+)<sup>(\d+)</sup>",
																																			 element.electronic_structure)]
	for eletron_subshell in element_electronic_structure:
			general_element_electronic[eletron_subshell] = 1.0
	return np.array([general_element_electronic[key] for key in general_electron_subshells])

def ofm(struct, is_ofm1=True, is_including_d=True):
	"""
	Generate OFM descriptor from a poscar file
	:param struct: pymatgen Structure object
	"""
	
	atoms = np.array([site.species_string for site in struct])
	local_xyz = []
	local_orbital_field_matrices = []
	for i_atom, atom in enumerate(atoms):

		coordinator_finder = VoronoiNN(cutoff=10.0)
		neighbors = coordinator_finder.get_nn_info(structure=struct, n=i_atom)

		site = struct[i_atom]
		center_vector = get_element_representation(atom)
		env_vector = np.zeros(32)
		
		atom_xyz = [atom]
		coords_xyz = [site.coords]
		for nn in neighbors:
			site_x = nn['site']
			w = nn['weight']
			site_x_label = site_x.species_string
			atom_xyz += [site_x_label]
			coords_xyz += [site_x.coords]
			neigh_vector = get_element_representation(site_x_label)
			d = np.sqrt(np.sum((site.coords - site_x.coords)**2))
			if is_including_d:
				env_vector += neigh_vector * w  / d
			else:
				env_vector += neigh_vector * w
		
		if is_ofm1:
			env_vector = np.concatenate(([1.0], env_vector))

		local_matrix = center_vector[None, :] * env_vector[:, None]
		local_matrix = np.ravel(local_matrix) # ravel to make 2024-Dimensional vector
		local_orbital_field_matrices.append(local_matrix)
		local_xyz.append({"atoms": np.array(atom_xyz), "coords": np.array(coords_xyz)})  

	local_orbital_field_matrices = np.array(local_orbital_field_matrices)
	material_descriptor = np.mean(local_orbital_field_matrices, axis=0)

	return {'mean': material_descriptor, 
					'locals': local_orbital_field_matrices, 
					'atoms': atoms,
					"local_xyz": local_xyz,
					'pm_struct': struct}



def get_ofm1_name():
	centers = ['s1', 's2', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6',
		'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10',
		'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7',
		'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'f14']
	envs = ['ofcenter', 's1', 's2', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6',
		'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10',
		'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7',
		'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'f14']


	# "s1-s1", "s1-s2", "s1-d5", "s1-f3", "s2-f14", "p1-center", "p1-d3", "p1-f1", "p1-f13", "p1-f14", "p2-d2", "p2-d10", "d6-f2", "d6-f3", "d6-f14", "d7-p5", "f4-d4", "f4-d5", "f4-f6", "f4-f14"
	# # the right format should be 
	# [center_size, 1] * [1, envs_size] = [center_size, envs_size]
	ofm_name = []
	for e in envs:
		for c in centers:
			ofm_name.append("{0}-{1}".format(c, e))
	return ofm_name


if __name__ == '__main__':

	# Testing
	# data = np.load("data.npy", encoding = 'latin1')
	# a = data[0]
	filename = "test.cif"
	with open(filename, 'r') as f:
		inp = f.read()
	struct = pm.Structure.from_str(inp, fmt="cif")
	test = ofm(struct=struct, is_including_d=True)

	# print(np.sum(test["mean"] - a["mean"]))