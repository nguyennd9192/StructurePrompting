import os, sys, abc, pickle
# sc_dir = "/media/nguyen/work/SmFe12_screening/source_code/"
# for ld, subdirs, files in os.walk(sc_dir):
# 	if os.path.isdir(ld) and ld not in sys.path:
# 		sys.path.append(ld)

import numpy as np
import pymatgen as pm
from abc import ABCMeta, abstractmethod

from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis import transition_state
from pymatgen.io.vasp import inputs, outputs
from pymatgen.util.plotting import pretty_plot


from ofm import ofm, get_ofm1_name
from general_lib import *
from plot import release_mem, ax_setting

class Feature(object):
	def __init__(self, name, saveat, config):
		self.name = name
		self.config = config
		self.saveat = saveat

		with open(config["filename"], 'r') as f:
			inp = f.read()
		struct = pm.Structure.from_str(input_string=inp, fmt='poscar') # poscar
		self.struct = struct

	@abstractmethod
	def reader_data(self): 
		pass 
	
	@abstractmethod
	def create_feature(self):
		pass

	@abstractmethod
	def saveto(self):
		if self.saveat is not None:
			tmp = self.create_feature()
			makedirs(self.saveat)
			with open(self.saveat, "wb") as f:
				pickle.dump(self, f)

		return

class OFMFeature(Feature):
	def __init__(self, name, saveat, config):
		# super(Feature, self).__init__(name=name, saveat=saveat, config=config)
		# super(Feature, self).__init__()
		Feature.__init__(self, name, saveat, config)
		# with open(config["filename"], 'r') as f:
		#         inp = f.read()
		# struct = pm.Structure.from_str(input_string=inp, fmt='poscar')    
		# self.struct = struct

	def reader_data(self):
		ofm_calc = ofm(struct=self.struct, 
				is_ofm1=self.config["is_ofm1"], 
				is_including_d=self.config["is_including_d"])

		self.feature = ofm_calc["mean"]
		self.locals = ofm_calc["locals"]
		self.atoms = ofm_calc["atoms"]
		self.natoms = len(self.atoms)
		self.local_xyz = ofm_calc["local_xyz"]
		self.pm_struct = ofm_calc["pm_struct"]

		self.feature_name = get_ofm1_name()
		# 'mean': material_descriptor, 
		# 'locals': local_orbital_field_matrices, 
		# 'atoms': atoms,
		# 'local_xyz': local_xyz,
		# 'pm_struct': struct

		# wyckoff_uniq_sites, wyckoff_all_sites = get_all_wickoff_idx(name=self.name)

		# self.wyckoff_uniq_sites = wyckoff_uniq_sites
		# self.wyckoff_all_sites = [wyckoff_all_sites[k] for k in range(self.natoms)]
		return 


	def create_feature(self):
		self.reader_data()
		return self.feature

	def get_feature_name(self):
		return get_ofm1_name()

class XRDFeature(Feature):
	def __init__(self, name, saveat, config):

		Feature.__init__(self, name, saveat, config)

		# with open(config["filename"], 'r') as f:
		#         inp = f.read()
		# struct = pm.Structure.from_str(input_string=inp, fmt='poscar')
		# self.struct = struct
		print (name)

	def reader_data(self):
		_struct = self.struct
		# sym = bandstructure.HighSymmKpath(structure=_struct)
		# sym.get_kpath_plot()

		# ana = analyzer.SpacegroupAnalyzer
		#test = structure_analyzer.VoronoiConnectivity(structure=structure)
		#print "get connection", test.get_connections()
		#test = GWvaspinputsets.SingleVaspGWWork(structure=structure, job= 'G0W0', spec=1.0)
		#print "create input = ", test.create_input()
		#print "=============================="
		#test = input(out[0])
		#print "velocities = ", test.natoms()
		#test = outputs.Outcar(out[0])
		#print test.charge
		#test = transition_state.NEBAnalysis(structures=structure)
		#test.get_plot()


		cal = pm.analysis.diffraction.xrd.XRDCalculator(wavelength="CuKa")
		pattern = cal.get_xrd_pattern(structure=_struct, 
			scaled=True, two_theta_range=(0, 90)).__dict__
		
		# pattern_keys = ['d_hkls', 'ydim', 'hkls', '_kwargs', '_args', 'y', 'x']
		self.d_hkls = pattern["d_hkls"]
		self.ydim = pattern["ydim"]
		self.hkls = pattern["hkls"]
		self._kwargs = pattern["_kwargs"]
		self._args = pattern["_args"]
		self.feature = pattern["y"] # # feature is intensity
		self.feature_name = pattern["x"] # # feature_name is two_theta


		# # x: 2_theta angle
		# # y: intensity
		this_plt = cal.get_xrd_plot(structure=_struct,  two_theta_range=(0, 90),
			fontsize=12,
			# title="XRD of {0}".format(get_basename(self.config["filename"])), 
			# savefig=self.config["savefig_at"], 
			# tight_layout=True
			)
		# this_plt.show()
		this_plt.tight_layout()
		makedirs(self.config["savefig_at"])
		this_plt.savefig(self.config["savefig_at"])
		release_mem(fig=this_plt)
		# fig = cal.plot_structures(structure=[_struct],  two_theta_range=(0, 90)
		#     # title="XRD of {0}".format(get_basename(self.config["filename"])), 
		#     # show=True, savefig=self.config["savefig_at"], tight_layout=True
		#     )

		return self.feature

	def create_feature(self):
		self.reader_data()
		return self.feature


class XenonPyFeature(Feature):
	def __init__(self, name, config):
		Feature.__init__(self, name, config)


	def reader_data(self):
		
		return ''

	def create_feature(self):
		self.reader_data()
		return ''


def load_feature(ft_file, prop):

	# # # prop: "feature", "locals", "atoms", "config"
	with open(ft_file, "rb") as f:
		load_feature = pickle.load(f)
		return load_feature.__dict__[prop]

def cfg_ofm(input_dir, subs_dir, current_job, is_including_d):
	basename = get_basename(subs_dir)

	if is_including_d:
		ft_type = "ofm1_with_d"
	else:
		ft_type = "ofm1_no_d"

	config = dict({"filename": subs_dir, 
		"ft_type": ft_type,
		"is_ofm1":True, "is_including_d": is_including_d})

	if "poscar" in basename:
		bn_ext = basename.replace("poscar", ft_type)
	else:
		bn_ext = "{0}.{1}".format(basename, ft_type)

	saveat = "{0}/feature/{1}/{2}/{3}".format(input_dir, 
			current_job, ft_type, bn_ext)

	return config, saveat



def cfg_xrd(input_dir, subs_dir, current_job, result_dir):
	ft_type = "xrd"
	basename = get_basename(subs_dir)

	config = dict({"filename": subs_dir, 
		"ft_type":ft_type})

	saveat = "{0}/feature/{1}/{2}/{3}".format(input_dir, 
			current_job, ft_type, basename.replace("poscar", ft_type))


	if "poscar" in basename:
		bn_ext = basename.replace("poscar", ft_type)
	else:
		bn_ext = "{0}.{1}".format(basename, ft_type)

	saveat = "{0}/feature/{1}/{2}/{3}".format(input_dir, 
			current_job, ft_type, bn_ext)

	savefig = "{0}/feature/{1}/{2}/{3}.pdf".format(result_dir, 
		current_job, ft_type, bn_ext)
	config["savefig_at"] = savefig

	return config, saveat

def process_ofm_xrd(input_dir, result_dir, current_job):
	subs_dirs = get_subdirs(sdir="{0}/origin_struct/{1}".format(input_dir, current_job))
	for subs_dir in subs_dirs:
		for is_including_d in [True, False]:
			config, saveat = cfg_ofm(input_dir=input_dir, 
					subs_dir=subs_dir,
					current_job=current_job,
					is_including_d=is_including_d)
			c = OFMFeature(name=subs_dir, saveat=saveat, config=config)
			feature = c.create_feature()
			feature_name = c.get_feature_name()
			c.saveto()

		config, saveat = cfg_xrd(input_dir=input_dir,
					subs_dir=subs_dir, 
					current_job=current_job, 
					result_dir=result_dir)
		c = XRDFeature(name=subs_dir, saveat=saveat, config=config)

		feature = c.create_feature()
		c.saveto()



def get_std_representation(input_dir, subs_dir, current_job, is_including_d):
	config, _ = cfg_ofm(input_dir=input_dir, 
						subs_dir=subs_dir,
						current_job=current_job,
						is_including_d=is_including_d)
	saveat = _.replace("/feature/", "/feature/standard/")

	c = OFMFeature(name=subs_dir, saveat=saveat, config=config)
	feature = c.create_feature()
	feature_name = c.get_feature_name()
	# c.saveto()
	return feature, feature_name

if __name__ == '__main__':
	is_get_ofm = True
	is_xrd = True

	# # for mix only
	mix_job = input_dir+"/origin_struct/mix_2-24" # mix
	jobs = ["mix_2-24/"+get_basename(k) for k in get_subdirs(mix_job)]
	print (jobs)

	# jobs = ["Sm2-Fe21-Ga3_std6000"] # CuAlZnTi_Sm2-Fe21-M3_std6000
	# # FOR INITIAL SUBSTITUTION STRUCTURES
	for current_job in jobs:
		
		input_dir = "/media/nguyen/work/SmFe12_screening/input"
		subs_dirs = get_subdirs(sdir="{0}/origin_struct/{1}".format(input_dir, current_job))

		if is_get_ofm:
			for is_including_d in [False, True]:
				for subs_dir in subs_dirs:
					config, saveat = cfg_ofm(input_dir=input_dir, 
							subs_dir=subs_dir,
							current_job=current_job,
							is_including_d=is_including_d)
					c = OFMFeature(name=subs_dir, saveat=saveat, config=config)
					feature = c.create_feature()
					feature_name = c.get_feature_name()
					c.saveto()
	
		if is_xrd:
			for subs_dir in subs_dirs:

				config, saveat = cfg_xrd(input_dir=input_dir,
							subs_dir=subs_dir, 
							current_job=current_job, 
							result_dir=result_dir)
				c = XRDFeature(name=subs_dir, saveat=saveat, config=config)

				feature = c.create_feature()
				c.saveto()


