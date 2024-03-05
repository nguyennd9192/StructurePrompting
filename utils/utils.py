from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import copy, os, pickle, sys, scipy

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.preprocessing import OneHotEncoder

from tensorflow.io import gfile

from utils.kernel_block_solver import BlockKernelSolver
from utils.small_cnn import SmallCNN
from utils.allconv import AllConv
from utils.uncertainty_regression import *
from utils.mixture_of_experts import NN_estimator

def load_pickle(filename):
  if not gfile.exists(filename):
    raise NameError("ERROR: the following data not available \n" + filename)
  data = pickle.load(gfile.GFile(filename, "rb"))
  return data

class Logger(object):
  """Logging object to write to file and stdout."""

  def __init__(self, filename):
    self.terminal = sys.stdout
    self.log = gfile.GFile(filename, "w")

  def write(self, message):
    self.terminal.write(message)
    self.log.write(message)

  def flush(self):
    self.terminal.flush()

  def flush_file(self):
    self.log.flush()

 
def create_checker_unbalanced(split, n, grid_size):
  """Creates a dataset with two classes that occupy one color of checkboard.

  Args:
    split: splits to use for class imbalance.
    n: number of datapoints to sample.
    grid_size: checkerboard size.
  Returns:
    X: 2d features.
    y: binary class.
  """
  y = np.zeros(0)
  X = np.zeros((0, 2))
  for i in range(grid_size):
    for j in range(grid_size):
      label = 0
      n_0 = int(n/(grid_size*grid_size) * split[0] * 2)
      if (i-j) % 2 == 0:
        label = 1
        n_0 = int(n/(grid_size*grid_size) * split[1] * 2)
      x_1 = np.random.uniform(i, i+1, n_0)
      x_2 = np.random.uniform(j, j+1, n_0)
      x = np.vstack((x_1, x_2))
      x = x.T
      X = np.concatenate((X, x))
      y_0 = label * np.ones(n_0)
      y = np.concatenate((y, y_0))
  return X, y


def flatten_X(X):
  shape = X.shape
  flat_X = X
  if len(shape) > 2:
    flat_X = np.reshape(X, (shape[0], np.product(shape[1:])))
  return flat_X 

def get_mldata(data_dir, name):
  """Loads data from data_dir.

  Looks for the file in data_dir.
  Assumes that data is in pickle format with dictionarmoeields data and target.


  Args:
    data_dir: directory to look in
    name: dataset name, assumes data is saved in the save_dir with filename
      <name>.pkl
  Returns:
    data and targets
  Raises:
    NameError: dataset not found in data folder.
  """
  dataname = name
  if dataname == "checkerboard":
    X, y = create_checker_unbalanced(split=[1./5, 4./5], n=10000, grid_size=4)
  else:
    filename = os.path.join(data_dir, dataname + ".pkl")
    data = load_pickle(filename)
    X = data["data"]
    y = data["target"]
    index=None
    if "index" in data:
      index = data["index"]
    if "keras" in dataname:
      X = X / 255
      y = y.flatten()
  return X, y, index

def parse_material_structure_data(
  df_ms_data, material_composition, list_of_elements,
  fe_threshold=0.070463
):
  """Parse material structure data into dataset of material families.

  Args:
    df_ms_data: Dataframe of calculated material structures
    material_composition: Material composition of calculated structures {Fe: 12, Sm:1, N:3}
    list_of_elements: List of elements are used to define the compositions, including symmetry index
    fe_threshold: Threshold of formation energy to determine which structures are stable

  Returns:
    Dataframe contains information of calculated material families
  """
  n_calculated_structures = len(df_ms_data)
  
  # Init onehot vector for material family
  df_ms_data[list_of_elements] = np.zeros((n_calculated_structures, len(list_of_elements)))
  df_ms_data[list_of_elements] = df_ms_data[list_of_elements].astype(int)
  
  # Update onehot vector for symmetry index
  enc = OneHotEncoder(handle_unknown='ignore')
  enc.fit((np.array(range(230))+1).reshape(-1,1))
  X_cat = enc.transform(df_ms_data["SYMM"].values.reshape(-1,1))
  df_symm = pd.DataFrame(X_cat.toarray(), columns=["SYMM_{}".format(i) for i in enc.categories_[0]], index=df_ms_data.index.values)
  df_ms_data = pd.concat([df_ms_data, df_symm], axis=1)
  
  # Update onehot vector for material composition
  set_name = []
  for element, n_elements in material_composition.items():
    df_ms_data[element] = np.ones(n_calculated_structures) * n_elements
    set_name.extend([element for i in range(n_elements)])
  
  set_names, elements, indices = [], []
  for ith, row in df_ms_data.iterrows():
      elements.append("|".join(set_name))
      set_names.append("|".join(set_name + ["SYMM_{}".format(row["SYMM"])]))
      indices.append("{}_{}_{}".format(element, row["Gen"], ith))
      
  df_ms_data["length"] = df_ms_data[list_of_elements].values.sum(axis=1)
  df_ms_data["Label"] = ["High" if i < fe_threshold else "Low" for i in df_ms_data["energy_substance_pa"].values]
  df_ms_data["set_name"] = set_names
  df_ms_data["element"] = elements
  df_ms_data["ID"] = indices
  
  df_ms_data = df_ms_data.set_index("ID")
  return df_ms_data

def filter_data(X, y, keep=None):
  """Filters data by class indicated in keep.

  Args:
    X: train data
    y: train targets
    keep: defaults to None which will keep everything, otherwise takes a list
      of classes to keep

  Returns:
    filtered data and targets
  """
  if keep is None:
    return X, y
  keep_ind = [i for i in range(len(y)) if y[i] in keep]
  return X[keep_ind], y[keep_ind]


def get_class_counts(y_full, y):
  """Gets the count of all classes in a sample.

  Args:
    y_full: full target vector containing all classes
    y: sample vector for which to perform the count
  Returns:
    count of classes for the sample vector y, the class order for count will
    be the same as long as same y_full is fed in
  """
  classes = np.unique(y_full)
  classes = np.sort(classes)
  unique, counts = np.unique(y, return_counts=True)
  complete_counts = []
  for c in classes:
    if c not in unique:
      complete_counts.append(0)
    else:
      index = np.where(unique == c)[0][0]
      complete_counts.append(counts[index])
  return np.array(complete_counts)


def flip_label(y, percent_random):
  """Flips a percentage of labels for one class to the other.

  Randomly sample a percent of points and randomly label the sampled points as
  one of the other classes.
  Does not introduce bias.

  Args:
    y: labels of all datapoints
    percent_random: percent of datapoints to corrupt the labels

  Returns:
    new labels with noisy labels for indicated percent of data
  """
  classes = np.unique(y)
  if len(classes) > 10:
    print("Warning!!! Reconsidering is_clf tag. The number of classes is larger than 10.")
  y_orig = copy.copy(y)
  indices = list(range(y_orig.shape[0]))
  np.random.shuffle(indices)
  sample = indices[0:int(len(indices) * 1.0 * percent_random)]
  fake_labels = []
  for s in sample:
    label = y[s]
    class_ind = np.where(classes == label)[0][0]
    other_classes = np.delete(classes, class_ind)
    np.random.shuffle(other_classes)
    fake_label = other_classes[0]
    assert fake_label != label
    fake_labels.append(fake_label)
  y[sample] = np.array(fake_labels)
  assert all(y[indices[len(sample):]] == y_orig[indices[len(sample):]])
  return y


def get_model(method, seed=13, is_search_params=True, 
  n_shuffle=10000, mt_kernel=None):
  """Construct sklearn model using either logistic regression or linear svm.

  Wraps grid search on regularization parameter over either logistic regression
  or svm, returns constructed model

  Args:
    method: string indicating scikit method to use, currently accepts logistic
      and linear svm.
    seed: int or rng to use for random state fed to scikit method

  Returns:
    scikit learn model
  """
  # TODO(lishal): extend to include any scikit model that implements
  #   a decision function.
  # TODO(lishal): for kernel methods, currently using default value for gamma
  # but should probably tune.
  # # for my building u_gp vs e_krr

  nn_libs = ["fully_connected", "moe", "LeNet"]
  if method=="u_gp":
    model = UncertainGaussianProcess(
      name=method, random_state=seed, 
      cv=5, search_param=is_search_params, 
      mt_kernel=None) 
    return model

  if method=="u_gp_mt" and mt_kernel is not None:
    model = UncertainGaussianProcess(name=method,
      random_state=seed, cv=10, 
      search_param=is_search_params, 
      mt_kernel=mt_kernel) 
    return model

  if method=="e_krr":
    model = UncertainEnsembleRegression(
      name=method, random_state=seed, 
      n_shuffle=n_shuffle, alpha=0.1, gamma=0.1,
      cv=10, score_method="kr", search_param=is_search_params,
      )
    return model

  if method=="u_knn":
    model = UncertainKNearestNeighbor(
      name=method, random_state=seed, 
      cv=10, search_param=is_search_params)
    return model


  if method in nn_libs:
    NN_kwargs = dict({"method":method, "n_inputs":85, "n_epoches":10,
      "batch_size":10, "lr": 0.01, "momentum":0.9})

    if method=="moe":
      NN_kwargs["num_experts"] = 5
      NN_kwargs["hidden_dim"] = 15


    model = NN_estimator(random_state=1, cv=3, n_times=3, 
      search_param=False,
      verbose=False, NN_kwargs=NN_kwargs)
    return model
    

  if method == "logistic":
    model = LogisticRegression(random_state=seed, multi_class="multinomial",
                               solver="lbfgs", max_iter=200)
    params = {"C": [10.0**(i) for i in range(-4, 5)]}
  elif method == "krr":
    # print ("method:", method)
    model = KernelRidge(kernel='rbf')
    params = {"alpha": [10.0**(i) for i in range(-3, 1)], 
      "gamma": [10.0**(i) for i in range(-2, 1)]}
  elif method == "logistic_ovr":
    model = LogisticRegression(random_state=seed)
    params = {"C": [10.0**(i) for i in range(-5, 4)]}
  elif method == "linear_svm":
    model = LinearSVC(random_state=seed)
    params = {"C": [10.0**(i) for i in range(-4, 5)]}
  elif method == "kernel_svm":
    model = SVC(random_state=seed)
    params = {"C": [10.0**(i) for i in range(-4, 5)]}
  elif method == "kernel_ls":
    model = BlockKernelSolver(random_state=seed)
    params = {"C": [10.0**(i) for i in range(-6, 1)]}
  elif method == "small_cnn":
    # Model does not work with weighted_expert or simulate_batch
    model = SmallCNN(random_state=seed)
    return model
  elif method == "allconv":
    # Model does not work with weighted_expert or simulate_batch
    model = AllConv(random_state=seed)
    return model
  else:
    raise NotImplementedError("ERROR: " + method + " not implemented")

  model = GridSearchCV(model, params, cv=3)
  return model


def calculate_entropy(batch_size, y_s):
  """Calculates KL div between training targets and targets selected by AL.

  Args:
    batch_size: batch size of datapoints selected by AL
    y_s: vector of datapoints selected by AL.  Assumes that the order of the
      data is the order in which points were labeled by AL.  Also assumes
      that in the offline setting y_s will eventually overlap completely with
      original training targets.
  Returns:
    entropy between actual distribution of classes and distribution of
    samples selected by AL
  """
  n_batches = int(np.ceil(len(y_s) * 1.0 / batch_size))
  counts = get_class_counts(y_s, y_s)
  true_dist = counts / (len(y_s) * 1.0)
  entropy = []
  for b in range(n_batches):
    sample = y_s[b * batch_size:(b + 1) * batch_size]
    counts = get_class_counts(y_s, sample)
    sample_dist = counts / (1.0 * len(sample))
    entropy.append(scipy.stats.entropy(true_dist, sample_dist))
  return entropy


def get_train_val_test_splits(X, y, max_points, seed, confusion, seed_batch,
                              split=(2./3, 1./6, 1./6), is_clf=None):
  """Return training, validation, and test splits for X and y.

  Args:
    X: features
    y: targets
    max_points: # of points to use when creating splits.
    seed: seed for shuffling.
    confusion: labeling noise to introduce.  0.1 means randomize 10% of labels.
    seed_batch: # of initial datapoints to ensure sufficient class membership.
    split: percent splits for train, val, and test.
    is_clf: None mean auto check target variable is numerical variable or categorical variable
  Returns:
    indices: shuffled indices to recreate splits given original input data X.
    y_noise: y with noise injected, needed to reproduce results outside of
      run_experiments using original data.
  """
  np.random.seed(seed)
  X_copy = copy.copy(X)
  y_copy = copy.copy(y)

  # Introduce labeling noise
    
  if is_clf:
    # # for classification tasks
    y_noise = flip_label(y_copy, confusion)
  else:
    # # for regression tasks
    y_noise = y_copy

  indices = np.arange(len(y))

  if max_points is None:
    max_points = len(y_noise)
  else:
    max_points = min(len(y_noise), max_points)
  train_split = int(max_points * split[0])
  val_split = train_split + int(max_points * split[1])
  assert seed_batch <= train_split

  # Do this to make sure that the initial batch has examples from all classes
  min_shuffle = 3
  n_shuffle = 0
  y_tmp = y_noise

  if is_clf:
    # # Need at least 4 obs of each class for 2 fold CV to work in grid search step, for clf
    while (any(get_class_counts(y_tmp, y_tmp[0:seed_batch]) < 4)
           or n_shuffle < min_shuffle):
      # # break if shuffle > 3times AND number of each class > 4
      np.random.shuffle(indices)
      y_tmp = y_noise[indices]
      n_shuffle += 1
  else:
    # # for regression
    while (n_shuffle < min_shuffle):
      np.random.shuffle(indices)
      n_shuffle += 1
  idx_train = indices[0:train_split]
  idx_val = indices[train_split:val_split]
  idx_test = indices[val_split:max_points]

  X_train = X_copy[idx_train]
  X_val = X_copy[idx_val]
  X_test = X_copy[idx_test]

  y_train = y_noise[idx_train]
  y_val = y_noise[idx_val]
  y_test = y_noise[idx_test]
  if is_clf:
    # Make sure that we have enough observations of each class for 2-fold cv
    assert all(get_class_counts(y_noise, y_train[0:seed_batch]) >= 4)
  # Make sure that returned shuffled indices are correct
  assert all(y_noise[indices[0:max_points]] ==
             np.concatenate((y_train, y_val, y_test), axis=0))
  return (indices[0:max_points], X_train, y_train,
          X_val, y_val, X_test, y_test, y_noise, idx_train, idx_val, idx_test)


def get_sept_train_val_test(X, y, X_test, y_test, 
          max_points, seed, confusion, seed_batch, split=[2./3, 1./3], is_clf=None):

  # # create y_noise, X_copy, y_copy (copy now for train+validation only)
  indices = np.arange(len(y))
  np.random.seed(seed)
  X_copy = copy.copy(X)
  y_copy = copy.copy(y)
  if is_clf:
    y_noise = flip_label(y_copy, confusion)
  else:
    y_noise = y_copy
  # # preparing max points used in training
  if max_points is None:
    max_points = len(y_noise)
  else:
    max_points = min(len(y_noise), max_points)

  # # shuffling indices
  y_tmp = y_noise
  min_shuffle = 3
  n_shuffle = 0
  if is_clf:
    while (any(get_class_counts(y_tmp, y_tmp[0:seed_batch]) < 4)
           or n_shuffle < min_shuffle):
      # # break if shuffle > 3times AND number of each class > 4
      np.random.shuffle(indices)
      y_tmp = y_noise[indices]
      n_shuffle += 1
  else:
    while (n_shuffle < min_shuffle):
      np.random.shuffle(indices)
      n_shuffle += 1

  # # count number of train, val. split now have only 2 terms
  train_split = int(max_points * split[0])
  val_split = train_split + int(max_points * split[1])

  idx_train = indices[0:train_split]
  idx_val = indices[train_split:val_split]
  idx_test = None

  X_train = X_copy[idx_train]
  X_val = X_copy[idx_val]

  y_train = y_noise[idx_train]
  y_val = y_noise[idx_val]


  # print(np.concatenate((y_train, y_val)))
  assert all(y_noise[indices[0:max_points]] == np.concatenate((y_train, y_val)))
  return (indices[0:max_points], X_train, y_train,
          X_val, y_val, X_test, y_test, y_noise, idx_train, idx_val, idx_test)














