
import os, re, yaml, sys, glob, ntpath
import numpy as np
from sys import argv
import pymatgen as pm
import re
from pymatgen.io.cif import CifParser
from general_lib_37 import *


if __name__ == '__main__': 
	uspex_output = [
			# "BESTgatheredPOSCARS", "BESTgatheredPOSCARS_order", #"goodStructures_POSCARS",
			# "gatheredPOSCARS_unrelaxed", "gatheredPOSCARS_order", 
			"gatheredPOSCARS"
			]

	for ea_type in uspex_output:
		for job in jobs:
			symstruct_file = get_uspex_dir(job=job, uspex_file="results1/{0}".format(ea_type))

			with open(symstruct_file, 'r') as f:
				inp = f.read()
			iter_ =  list(re.finditer("EA", inp))
			indices = [m.start(0) for m in iter_]
			symstruct_savedir = get_uspex_dir(job=job, uspex_file="ML/symstruct/{0}".format(ea_type))
			for ith, indice in enumerate(indices):
				if ith != len(indices)-1:
					one_structure = inp[indice:indices[ith+1]]
				else:
					one_structure = inp[indice:]

				# data = cb.from_string(one_structure).get_structures()
				
				final = one_structure.find(" ")
				struct_name = one_structure[:final]
				savedir = symstruct_savedir + "/{0}.poscar".format(struct_name)# cif
				makedirs(savedir)

				with open(savedir, "w") as save_f:
					save_f.write(one_structure)















