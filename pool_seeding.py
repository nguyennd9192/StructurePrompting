import scipy.io
import numpy as np
import os, copy, glob, random
import numpy as np
import pandas as pd
from general_lib_37 import *

from symmetry_generation import *


def get_subdirs(sdir):
	subdirs = glob.glob(sdir+"/*")
	return subdirs

def Individuals2csv(filename):
	# # # # # # # # # # # # 
	# # convert Individuals uspex to csv
	# # # # # # # # # # # # 
	tocsv_file = filename + ".csv" 
	tof = open(tocsv_file, 'w')

	input_file = copy.copy(filename)
	input_file = input_file.replace("results1/Individuals", "INPUT.txt")

	with open(input_file, 'r') as f:
		lines = f.readlines()
		for ith, line in enumerate(lines):
			if "atomType" in line:
				comps = lines[ith+1]
				break
		print ("comps", comps.split())
	
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
	return tocsv_file, comps_array, comps_list

def gen_from_pool():
	indi_file = get_uspex_dir(job=job, uspex_file="results1/Individuals")

	if os.path.isfile(indi_file):
		# # link to the predefined storage by getting name of the task
		indi_csv, comps_array, comps_list = Individuals2csv(filename=indi_file)
		the_pool_dir = "/home/nguyen/work/USPEX/application/archive/SmFe12/the_pool/oqmd/"

		comps_array.remove("Fe")
		comps_array.remove("Sm")

		for e in comps_array:
			# # from composition to storage dir
			pool_dir = "{0}/{1}".format(the_pool_dir, "".join(["Sm", "Fe", e]))
			struc_in_pool = get_subdirs(sdir=pool_dir)
			
			indi_df = pd.read_csv(indi_csv, index_col="ID")
			n_individuals = len(indi_df)
			this_gen = np.nanmax(indi_df["Gen"])
			print ("pool_dir:", pool_dir)

			print ("n in pool:", len(struc_in_pool))
			print ("this gen:", this_gen)

			the_seed_file = pwd+"/Seeds/POSCARS_{}".format(int(this_gen + 1))
			if not os.path.isfile(the_seed_file):
				n_seeds = 5
				seed_files = random.sample(struc_in_pool, n_seeds)
				with open(the_seed_file, 'w') as result:
					for ith, seed_file in enumerate(seed_files):
						print ("File number:", ith)
						with open(seed_file, "r") as f:
							result.write(f.read())


def gen_sub_symm():
	# job = "base/1e-4/Sm1Fe12"
	for job in jobs:
		indi_file = get_uspex_dir(job=job, uspex_file="results1/Individuals")
		seed_dir = get_uspex_dir(job=job, uspex_file="Seeds")

		if os.path.isfile(indi_file):
			# # 
			indi_csv, comps_array, comps_list = Individuals2csv(filename=indi_file)
			indi_df = pd.read_csv(indi_csv, index_col="ID")
			n_individuals = len(indi_df)
			this_gen = np.nanmax(indi_df["Gen"])
			the_seed_file = seed_dir+"/POSCARS_{}".format(int(this_gen + 1))



			pool_dir = "{0}/{1}/ML/pool_dir/gen{2}".format(result_dir, job,  int(this_gen + 1))
			print (pool_dir)
			# generate_subgroup_path(savedir=pool_dir)
			# generate_subgroup_path(savedir=gen_dir, eps=0.2)


			
			makedirs(the_seed_file)
			seed_files = get_subdirs(pool_dir)
			write_seed_file(write_to=the_seed_file, seed_files=seed_files)



if __name__ == '__main__':
	gen_from_pool()
	# gen_sub_symm()


















