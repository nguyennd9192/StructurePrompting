
from general_lib_37 import * 
from pymatgen.transformations.advanced_transformations import *
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Specie, DummySpecie
from pyxtal.symmetry import Group
from pyxtal.viz import display_crystals
from pyxtal import pyxtal 
from pymatgen.vis.structure_vtk import StructureVis, make_movie, MultiStructuresVis
from pymatgen.vis.structure_chemview import quick_view
import subprocess
# from pymatgen.symmetry.site_symmetries import get_site_symmetries
from pymatgen.symmetry.structure import *
from pymatgen.transformations.advanced_transformations import CubicSupercellTransformation



def generate_descendant(parent_struct, H, savedir, ith_trial, eps, is_primcell, n_site_prime, init_symm):
	pmg_struct, new_structure, pmg_symm = None, None, None
	try:
		new_structure = parent_struct.subgroup_once(eps=eps, H=H, perms=None, group_type='k', 
					# # t, k--please no use t+k
					max_cell=2, min_cell=1, mut_lat=True, ignore_special=False)
		# lattice.mutate(degree=eps, frozen=True)
		pmg_struct = new_structure.to_pymatgen()

		pmg_symm = pmg_struct.get_space_group_info(symprec=symprec, angle_tolerance=5.0)[1]
		if is_primcell:
			prim_cell = pmg_struct.get_primitive_structure()
			prim_cell_symm_info = prim_cell.get_space_group_info(symprec=symprec, angle_tolerance=5.0)
			# print (prim_cell_symm_info)
			prim_cell_symm = prim_cell_symm_info[1]
			# assert pmg_symm == prim_cell_symm
			assert len(prim_cell.species) == n_site_prime
			pmg_struct = copy.copy(prim_cell)
		
		saveat = "{0}/{1}_trial_{2}.poscar".format(savedir, pmg_symm, ith_trial) 
		makedirs(saveat)

		pmg_struct.to(fmt="poscar", filename=saveat)
	except Exception as e:
		print ("Error: ", e)
		pass 

	return pmg_struct, new_structure, pmg_symm

	
 
def get_sub_groups(symm):
	g = Group(symm)
	subgroups = g.get_max_t_subgroup()['subgroup']
	set_subgroups = set(subgroups)
	return subgroups

def trans_by_path(ref_struct):
	structs = ref_struct._get_subgroup_ids(H=None, group_type="t", idx=None, max_cell=1, min_cell=1)

	# structs = p.get_transition_by_path(ref_struct, path, d_tol, d_tol2=0.5, N_images=2, both=False)
	# structs = ref_struct.get_transition(ref_struc=ref_struct, d_tol=1.0, d_tol2=0.3, 
	# 	N_images=2, max_path=30, both=False)
	return structs


def generate_subgroup_path(savedir, eps, is_primcell, input_struct_file):
	# # read from input 
	input_struct = read_poscar(filename=input_struct_file)
	input_symm = input_struct.get_space_group_info(symprec=symprec, angle_tolerance=5.0)

	n_site_prime = len(input_struct.species)
	# # derive subgroups
	inp_symms = get_sub_groups(symm=input_symm[1])

	
	input_struct_pyxtal = pyxtal()
	input_struct_pyxtal.from_seed(seed=input_struct, tol=1e-2, a_tol=5.0) # style='spglib'

	gen_symm = []
	all_gene = []
	xrds = []

	n_candidates = 0
	n_attempt = 0

	while n_candidates < max_symm_gen:
		H = None
		parent_struct = copy.copy(input_struct_pyxtal)
		sub_savedir = "{0}/".format(savedir)
		for ith_symm in range(n_symm_levels):
			child_pmg_struct, child_xtal_structure, child_symm = generate_descendant(
				parent_struct=parent_struct, H=None, savedir=sub_savedir,
				ith_trial=n_candidates, eps=eps, is_primcell=is_primcell,
				n_site_prime=n_site_prime, init_symm=input_symm[1],
				)

			if child_symm is not None:
				H = None # copy.copy(child_symm)
				parent_struct = copy.copy(child_xtal_structure)
				sub_savedir += "/{0}".format(child_symm)
				
		candidates = [os.path.join(path, name) for path, subdirs, files in os.walk(savedir) for name in files]
		n_candidates = len(candidates)
		n_attempt += 1
		if n_attempt > max_one_attempt:
			break

