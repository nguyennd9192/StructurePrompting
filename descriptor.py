import numpy as np
import pymatgen as pm
from abc import ABCMeta, abstractmethod
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from lib import *
from pymatgen.analysis import transition_state
from pymatgen.io.vasp import inputs, outputs
from pymatgen.util.plotting import pretty_plot


from ofm import ofm, get_ofm1_name

class Feature(object):
	def __init__(self, name, saveat, config):
		self.name = name
		self.config = config
		self.saveat = saveat

		with open(config["filename"], 'r') as f:
			inp = f.read()
		
		struct = pm.core.structure.Structure.from_str(input_string=inp, fmt='poscar') 

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
		Feature.__init__(self, name, saveat, config)

	def reader_data(self):
		try:
			ofm_calc = ofm(struct=self.struct, 
				is_ofm1=self.config["is_ofm1"], 
				is_including_d=self.config["is_including_d"])
			self.feature = ofm_calc["mean"]
			self.locals = ofm_calc["locals"]
			self.atoms = ofm_calc["atoms"]
			self.natoms = len(self.atoms)
			self.local_xyz = ofm_calc["local_xyz"]
			self.pm_struct = ofm_calc["pm_struct"]
		except Exception as e:
			self.feature = None
			self.locals = None
			self.atoms = None
			self.natoms = None
			self.local_xyz = None
			pass
		
		self.feature_name = get_ofm1_name()
		# 'mean': material_descriptor, 
		# 'locals': local_orbital_field_matrices, 
		# 'atoms': atoms,
		# 'local_xyz': local_xyz,
		# 'pm_struct': struct
		return 


	def create_feature(self):
		self.reader_data()
		return self.feature

	def get_feature_name(self):
		return get_ofm1_name()

class XRDFeature(Feature):
	def __init__(self, name, saveat, config):

		Feature.__init__(self, name, saveat, config)

	def reader_data(self):
		_struct = self.struct

		cal = pm.analysis.diffraction.xrd.XRDCalculator(wavelength="CuKa")
		pattern = cal.get_pattern(structure=_struct, 
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
		this_plt = cal.get_plot(structure=_struct,  two_theta_range=(0, 90),
			fontsize=12,
			)
		# this_plt.show()
		if self.config["xs"] is not None:
			this_plt.vlines(self.config["xs"], ymin=0.0, ymax=100, color="red")
			this_plt.xticks(self.config["xs"], rotation='vertical')
			this_plt.title(get_basename(self.config["savefig_at"]))
		this_plt.tight_layout()
		makedirs(self.config["savefig_at"])
		this_plt.savefig(self.config["savefig_at"])
 
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

def cfg_ofm(subs_dir, current_job, is_including_d):
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

	if input_dir in subs_dir:
		# # save for oqmd only
		saveat = "{0}/oqmd_feature/{1}/{2}/{3}".format(input_dir, current_job, ft_type, bn_ext)
	else:
		saveat = "{0}/{1}/ML/feature/{2}/{3}".format(result_dir, 
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
			config, saveat = cfg_ofm( 
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
	config, _ = cfg_ofm(subs_dir=subs_dir,
						current_job=current_job,
						is_including_d=is_including_d)
	saveat = _.replace("/feature/", "/feature/standard/")

	c = OFMFeature(name=subs_dir, saveat=saveat, config=config)
	feature = c.create_feature()
	feature_name = c.get_feature_name()
	# c.saveto()
	return feature, feature_name

if __name__ == '__main__':
	is_xrd = False
	for job in jobs:
		# # for uspex
		sdir="{0}/{1}/ML/symstruct/gatheredPOSCARS".format(result_dir, job)

		subs_dirs = get_subdirs(sdir=sdir)
		print (sdir)
		print (subs_dirs)

		for subs_dir in subs_dirs:
			if ".png" in subs_dir:
				continue
			config, saveat = cfg_ofm( 
					subs_dir=subs_dir,
					current_job=job,
					is_including_d=True)

			if not os.path.isfile(saveat):
				c = OFMFeature(name=subs_dir, saveat=saveat, config=config)
				feature = c.create_feature()
				print ("Summ all feature: ", np.sum(feature))
				feature_name = c.get_feature_name()
				c.saveto()
			else:
				print ("File existed.")
	
		if is_xrd:
			for subs_dir in subs_dirs:
				if ".png" in subs_dir:
					continue
				config, saveat = cfg_xrd(input_dir=input_dir,
							subs_dir=subs_dir,
							current_job=job, 
							result_dir=result_dir)
				config["xs"] = [30, 33, 36.0, 37.5, 42.5, 47.5, 71, 80]
				c = XRDFeature(name=subs_dir, saveat=saveat, config=config)

				feature = c.create_feature()
				c.saveto()







