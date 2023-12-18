
from general_lib_37 import * 
from pymatgen.core.structure import Structure
from pyxtal import pyxtal 
from pyxtal.lattice import Lattice

import nevergrad as ng
from robocrys import StructureCondenser, StructureDescriber
import qmpy

from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from timeout_decorator import timeout
from symmetry_generation import make_supercell

from pymatgen.transformations.advanced_transformations import SuperTransformation, MonteCarloRattleTransformation
import re


@timeout(20)

def structure_training(org_structure, 
		Sm_frac_coords, Fe_frac_coords, X_frac_coords,
		species, numIons, 
		Sm_labels, Fe_labels, X_labels, 
		symm_group, is_return_structure=False):
	structure = read_poscar(filename=org_structure)
	

	lattice = structure.lattice
	sg_analyzer = SpacegroupAnalyzer(structure)
	lattice_type = sg_analyzer.get_lattice_type()
	pyxtal_cell = Lattice.from_para(
		a=lattice.a, b=lattice.b, c=lattice.c, 
		alpha=lattice.alpha, beta=lattice.beta, gamma=lattice.gamma, 
		ltype=lattice_type)



	Sm_sites = dict({})
	wyckoff_letters = []

	for label, coords in zip(Sm_labels, Sm_frac_coords):
		# sites.extend(lattice.get_wyckoff_positions(label, coords, species=["Sm"]*len(coords) ))
		# wyckoff_site = {"species": "Sm", "abc": coords, "label": label}
		# wyckoff_letters.append(label)
		rename_label = "{}{}".format(label[1:], label[:1])
		Sm_sites[rename_label] = coords

	Fe_sites = dict({})
	for label, coords in zip(Fe_labels, Fe_frac_coords):
		rename_label = "{}{}".format(label[1:], label[:1])
		Fe_sites[rename_label] = coords

	sites = [Sm_sites, Fe_sites]

	if X_labels is not None:
		X_sites = dict({})
		for label, coords in zip(X_sites, X_frac_coords):
			rename_label = "{}{}".format(label[1:], label[:1])
			X_sites[rename_label] = coords
			sites.append(X_sites)

	# elements = ['Sm', 'Fe']
	# composition = [2, 24]

	# sites = [{'2b': [0. , 0. , 0.5]}, # # original: "b2"
	# 		{'8h': [0.230331, 0.230331, 0.      ], 
	# 		'16k': [0.337969, 0.162031, 0.25    ]}]

	pert_structure = None
	s = pyxtal()
	try:
		s.from_random(3, symm_group, species, numIons, 
			max_count=20,
			lattice=pyxtal_cell, sites=sites)

		pert_structure = s.to_pymatgen()
		if is_return_structure:
			return pert_structure
		# rattle_std = 2  # Adjust as needed
		# min_distance = 2.6
		# trans = MonteCarloRattleTransformation(rattle_std=rattle_std, min_distance=min_distance)
		# # Step 5: Create a super transformation
		# # This is needed to apply multiple transformations in sequence
		# # trans = SuperTransformation(transformations)

		# # Step 6: Apply the transformation to the structure
		# pert_structure = trans.apply_transformation(pert_structure)

		analyzer = SpacegroupAnalyzer(structure=pert_structure, symprec=symprec)
		# pert_structure =  analyzer.get_symmetrized_structure()

		pert_structure = analyzer.get_conventional_standard_structure(
			international_monoclinic=True)
		# a += 1 

		if is_return_structure:
			return pert_structure

		condenser = StructureCondenser(
				use_conventional_cell=False, # True with symmetry site, False with other
				symprec=symprec)
		describer = StructureDescriber()
		condensed_structure = condenser.condense_structure(pert_structure)
		description = describer.describe(condensed_structure)

	except Exception as e:
		print (e)
		print ("=========")
		return 100.0






	# structure = Structure(lattice, sites)

	# spg_ops = SpacegroupOperations(
 #            spacegroup.get_space_group_symbol(),
 #            spacegroup.get_space_group_number(),
 #            spacegroup.get_symmetry_operations(),
 #        )
	# print (spg_ops)
	# structure = SymmetrizedStructure(structure, 
	# 	spacegroup=spg_ops, 
	# 	equivalent_positions=sites, wyckoff_letters=wyckoff_letters)
	# print (structure)


	symm_sites = describer._da.sites
	score = 0



	for symm_site in symm_sites:
		element = symm_sites[symm_site]["element"]
		site_type = symm_sites[symm_site]["geometry"]["type"]
		likeness = round(symm_sites[symm_site]["geometry"]["likeness"], 2)
		sym_labels = symm_sites[symm_site]["sym_labels"][0]

		print (element, site_type)
		print ("========")
		if "coordinate" in site_type:
			n_coord = float(site_type[:site_type.find("-")])
			score -= n_coord 
			if element == "Sm":
				if n_coord >= 16 and n_coord <= 20:
					score -= n_coord
				else:
					score += n_coord
			if element == "Fe":
				if n_coord >= 9 and n_coord <= 14:
					score -= n_coord
				else:
					score += n_coord



		if element == "Sm":
			if site_type == "square co-planar":
				score += 10
			if site_type == "body-centered cubic":
				score += 10

			if site_type == "cuboctahedral":
				score -= 10
			
		elif element == "Fe":
			if site_type == "body-centered cubic":
				score -= 10 
			if site_type == "cuboctahedral":
				score -= 10 
			else:
				score += 10 
	print ("score: ", score)
	return score


def optimize_query(filename, 
				species, numIons,
				symm_group, saveat):
	
	qm_struct = qmpy.io.poscar.read(filename)
	qm_struct.group_atoms_by_symmetry()
	qm_struct.symmetrize(tol=symprec, angle_tol=5)

	# init_frac_coords = np.array(qm_struct.frac_coords)
	tmp = list(set(species) - set(["Fe", "Sm"]))
	if len(tmp) == 1:
		X_element = tmp[0]
	else:
		X_element = None


	wyckoff_info =dict({
		"Fe": dict({"labels": [], "frac_coords":[]}),
		"Sm": dict({"labels": [], "frac_coords":[]}),
		})

	if X_element is not None:
		wyckoff_info[X_element] = dict({"labels": [], "frac_coords":[]})

	
	for i, qm_site in enumerate(qm_struct.sites):
		element = str(qm_site.label)
		wyckoff_id = str(qm_site.wyckoff)
		qm_f_coord = list(qm_site.coord)

		if wyckoff_id not in wyckoff_info[element]["labels"]:
			wyckoff_info[element]["frac_coords"].append(qm_f_coord)
			wyckoff_info[element]["labels"].append(wyckoff_id)

	Sm_sites_shape = np.array(wyckoff_info["Sm"]["frac_coords"]).shape
	Sm_frac_coords_search_space = ng.p.Array(
		init=wyckoff_info["Sm"]["frac_coords"], 
				# shape=init_frac_coords.shape,
				lower=np.full(Sm_sites_shape, 0.0), upper=np.full(Sm_sites_shape, 1.0),
					)
	Fe_sites_shape = np.array(wyckoff_info["Fe"]["frac_coords"]).shape
	Fe_frac_coords_search_space = ng.p.Array(
		init=wyckoff_info["Fe"]["frac_coords"], 
						# shape=init_frac_coords.shape,
						lower=np.full(Fe_sites_shape, 0.0), upper=np.full(Fe_sites_shape, 1.0),
							)

	if X_element is not None:
		X_sites_shape = np.array(wyckoff_info[X_element]["frac_coords"]).shape
		X_frac_coords_search_space = ng.p.Array(
			init=wyckoff_info[X_element]["frac_coords"], 
							lower=np.full(X_sites_shape, 0.0), upper=np.full(X_sites_shape, 1.0),
								)
		X_labels = wyckoff_info[X_element]["labels"]
	else:
		X_frac_coords_search_space = None
		X_labels = None


	parametrization = ng.p.Instrumentation(
		org_structure=filename,
		Sm_frac_coords=Sm_frac_coords_search_space,
		Fe_frac_coords=Fe_frac_coords_search_space,
		X_frac_coords=X_frac_coords_search_space,


		symm_group=symm_group,
		Sm_labels=wyckoff_info["Sm"]["labels"], 
		Fe_labels=wyckoff_info["Fe"]["labels"],
		X_labels=X_labels,


		is_return_structure=False,
		species=species, numIons=numIons,

				)


	# optimizer = ng.optimizers.NGOpt15(parametrization=parametrization, budget=20)

	# optimizer = ng.optimizers.BayesOptim(init_budget=10, pca=True, n_components=2)
	# print (dir(optimizer))
	optimizer = ng.optimizers.NGOpt39(parametrization=parametrization, budget=10)

	recommendation = optimizer.minimize(structure_training)

	best_kwargs = recommendation.kwargs
	print (best_kwargs)

	opt_structure = structure_training(
		org_structure=filename,
		Sm_frac_coords=best_kwargs["Sm_frac_coords"],
		Fe_frac_coords=best_kwargs["Fe_frac_coords"],
		X_frac_coords=best_kwargs["X_frac_coords"],


		symm_group=symm_group,
		Sm_labels=best_kwargs["Sm_labels"], 
		Fe_labels=best_kwargs["Fe_labels"],
		X_labels=best_kwargs["X_labels"],


		is_return_structure=True,
		species=species, numIons=numIons,


		)

	makedirs(saveat)

	if opt_structure is not None and not isinstance(opt_structure, float):
		opt_structure.to(fmt="poscar", filename=saveat)
		print ("===================")


def symmetry_structure_generator(gen_dir, initial_candidate, 
			group, species, numIons, n_gen):
	factor = 0.8
	input_struct = read_poscar(filename=initial_candidate)
	# input_struct = make_supercell(struct=input_struct, scale=2)
	lattice = input_struct.lattice
	sg_analyzer = SpacegroupAnalyzer(input_struct)
	lattice_type = sg_analyzer.get_lattice_type()
	

	if species is None or numIons is None:
		compositions = str(input_struct.composition).split(" ")
		r = re.compile("([a-zA-Z]+)([0-9]+)")
		tmp = np.array([r.match(composition).groups() for composition in compositions])
		species = [str(_) for _ in tmp[:, 0]]
		numIons = [int(_) for _ in tmp[:, 1]]

	if group is None:
		group = input_struct.get_space_group_info(symprec=symprec, angle_tolerance=5.0)[1]


	parent_struct = pyxtal()
	parent_struct.from_seed(seed=input_struct, tol=1e-2, a_tol=5.0) # style='spglib'

	random_structure = pyxtal()
	for trial in range(n_gen):
		# pyxtal_cell = Lattice.from_para(
		# 	# a=2*lattice.a, b=lattice.b, c=lattice.c, 
		# 	# alpha=lattice.alpha, beta=lattice.beta, gamma=lattice.gamma, 
		# 	a=8.4, b=8.4, c=4.8, 
		# 	alpha=90, beta=90, gamma=90, 

		# 	ltype=lattice_type)
		try:
			random_structure.from_random(
					dim=3,
					group=group, species=species, numIons=numIons,
					factor=factor, thickness=None,
					area=None, lattice=None, # pyxtal_cell
					sites=None,
					conventional=True, t_factor=1.0,
					max_count=3000,	torsions=None,
					force_pass=False,	block=None,
					num_block=None, seed=None,
					tm = None, use_hall=False,
				) #
			print (random_structure)
			Fe_sites = random_structure._get_coords_and_species()
		
			pmg_struct = random_structure.to_pymatgen()

			pmg_symm = pmg_struct.get_space_group_info(symprec=symprec, angle_tolerance=5.0)[1]
			saveat = gen_dir +"/{0}_{1}.poscar".format(group, trial) # Sm2Fe24N4 Sm2Fe24_30-230
			makedirs(saveat)
			pmg_struct.to(fmt="poscar", filename=saveat)

			init_dir = gen_dir +"/Sm2Fe24_139_opt/{0}_{1}.poscar".format(group, trial) # Sm2Fe24N4_opt
			opt_saveat = saveat.replace(".poscar", "_opt.poscar")
			optimize_query(filename=saveat, 
				species=species, numIons=numIons,
				symm_group=group, saveat=saveat)
			print ("Optimized. Symmetry generated saveat: ", saveat)
		# prim_cell_symm_info = prim_cell.get_space_group_info(symprec=symprec, angle_tolerance=5.0)
		# # print (prim_cell_symm_info)
		# prim_cell_symm = prim_cell_symm_info[1]
		except Exception as e:
			print (e)
			pass







