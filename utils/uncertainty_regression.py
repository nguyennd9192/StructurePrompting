# Copyright 2017 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Block kernel lsqr solver for multi-class classification."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import copy
import math

import numpy as np
import scipy.linalg as linalg
from scipy.sparse.linalg import spsolve
from sklearn import metrics
from regression import RegressionFactory, CV_predict_score
from sklearn.preprocessing import MinMaxScaler
from combination_generator import CombinationGeneratorFactory
from scipy import stats



# from statsmodels.distributions.empirical_distribution import ECDF
 
class UncertainEnsembleRegression(object):
  def __init__(self, name, 
        random_state=1, 
        n_shuffle=10000,
        alpha=0.1, gamma=0.1,
        cv=3,
        score_method="kr", search_param=True, # # GaussianProcess
):
    self.alpha = alpha
    self.gamma = gamma
    self.name = name 

    self.search_param = search_param
    self.kernel = 'rbf'
    self.coef_ = None
    self.cv = cv
    self.n_shuffle = n_shuffle
    self.random_state = random_state
    self.estimator = None

  def fit(self, X_train, y_train, sample_weight=None):
    # # in fit function
    # # just return estimator with best param with X_train, y_train
    np.random.seed(self.random_state)
    n_features = X_train.shape[1]

    if self.gamma is None:
      self.gamma = 1./n_features

    self.X_train = X_train
    self.y_train = y_train

    # # currently, anytime we fit the estimator with X_train, y_train
    # # searching for parameter
    if self.estimator is None: # # for not always search parameters:
      estimator, GridSearchCV = RegressionFactory.get_regression(
          method=self.name, 
          kernel=self.kernel, alpha=self.alpha, gamma=self.gamma, 
          search_param=self.search_param, X=X_train, y=y_train,  
          cv=self.cv)

      self.GridSearchCV = GridSearchCV
      self.estimator = estimator
    self.estimator.fit(X_train, y_train)
    return self.estimator

  def predict(self, X_val, get_pred_vals=False):
    indices = list(range(len(self.y_train)))
    multiple_sub_indices = CombinationGeneratorFactory.get_generator(
        method="bagging", items=indices, 
        sample_size=0.6, n_shuffle=self.n_shuffle
    )

    X_train_copy = copy.copy(self.X_train)
    y_train_copy = copy.copy(self.y_train)
    estimator_copy = copy.copy(self.estimator)

    y_val_preds = []
    for sub_indices in multiple_sub_indices:
      X_context = X_train_copy[sub_indices, :]
      y_context = y_train_copy[sub_indices]
      
      #######################
      # Fit the context model
      estimator_copy.fit(X_context, y_context)
      y_val_pred = estimator_copy.predict(X_val)
      # print("len_test_set:", len(y_val_pred))
      y_val_preds.append(y_val_pred)
      # r2 = r2_score(y_context, y_context_pred)
      # mae = mean_absolute_error(y_context, y_context_pred)
        #######################
    if get_pred_vals:
      return y_val_preds
    else:
      return np.mean(y_val_preds, axis=0)

  def score(self, X_val, y_val):
    y_pred = self.predict(X_val, get_pred_vals=False)
    val_acc = metrics.r2_score(y_val, y_pred)
    return val_acc

  def predict_proba(self, X, is_norm=True):
    # # large variance -> probability to be observed small -> sorting descending take first
    # # small variance -> probability to be observed large 
    y_val_preds = self.predict(X, get_pred_vals=True)
    # print("y_val_preds.shape", np.array(y_val_preds).shape)
    # # normalize variance to 0-1

    # # canonical method
    var = np.var(np.array(y_val_preds), axis=0).reshape(-1, 1)

    # # fitting with mixture gaussian, find cummulative
    # ecdf = ECDF(sample)
    # var = []
    # y_val_preds_T = np.array(y_val_preds).T
    # for y_val_pred in y_val_preds_T:
    #   # print (y_val_pred, len(y_val_pred))

    #   kernel = stats.gaussian_kde(y_val_pred)
    #   yref = np.linspace(min(y_val_pred),max(y_val_pred),100)
    #   pdf = kernel(yref).T

    #   # # find peaks
    #   peaks, _ = find_peaks(pdf)
    #   pr = peak_prominences(pdf, peaks)[0]
    #   argmax1, argmax2 = np.argsort(pr)[-2:] # # two largest prominence
    #   y_pred_peak1, y_pred_peak2 = yref[argmax1], yref[argmax2]

    #   # # variance between two largest peaks
    #   v = y_pred_peak2 - y_pred_peak1
    #   var.append([v])
    # var = np.array(var)

    var_norm = MinMaxScaler().fit_transform(X=var).ravel()

    # var_norm = var.reshape(-1, 1)
    # prob = 1 / (var_norm)
    if is_norm:
      return var_norm
    else:
      return var

  def best_score_(self, X=None, y=None):
    estimator, GridSearchCV = RegressionFactory.get_regression(
        method=self.name, 
        kernel=self.kernel, alpha=self.alpha, gamma=self.gamma, 
        search_param=self.search_param, X=X, y=y,  
        cv=self.cv) 

    r2, r2_std, mae, mae_std = CV_predict_score(estimator, X, y, 
                n_folds=3, n_times=3, score_type='r2')
    return mae


class UncertainGaussianProcess(object):
  def __init__(self, name,
        random_state=1, kernel="rbf",
        cv=3, search_param=False,
        mt_kernel=None):
    self.name = name 
    self.search_param = search_param
    self.kernel = kernel
    self.cv = cv
    self.estimator = None
    self.random_state = random_state
    self.mt_kernel = mt_kernel


  def fit(self, X_train, y_train, sample_weight=None):
    # # in fit function
    # # just return estimator with best param with X_train, y_train
    np.random.seed(self.random_state)


    self.X_train = X_train
    self.y_train = y_train
    if self.estimator is None:
      estimator, GridSearchCV = RegressionFactory().get_regression(
          method=self.name, 
          kernel=self.kernel, alpha=None, gamma=None, # # rbf
          search_param=self.search_param, X=X_train, y=y_train,  
          cv=self.cv, mt_kernel=self.mt_kernel) # mt_kernel=self.mt_kernel
      self.estimator = estimator
      self.GridSearchCV = GridSearchCV
    self.estimator.fit(X_train, y_train)
    return self.estimator

  def predict(self, X_val, get_variance=False):
    y_val_pred, y_val_pred_std = self.estimator.predict(X_val, 
      return_std=True, return_cov=False)
    if get_variance:
      return y_val_pred, y_val_pred_std 
    else:
      return y_val_pred

  def score(self, X_val, y_val):
    y_pred = self.predict(X_val, get_variance=False)    
    val_acc = metrics.r2_score(y_val, y_pred)
    return val_acc

  def predict_proba(self, X, is_norm=True):
    # # large variance -> probability to be observed small -> sorting descending take first
    # # small variance -> probability to be observed large 
    y_val_preds, y_val_pred_std = self.predict(X, get_variance=True)

    # # normalize variance to 0-1
    var_norm = MinMaxScaler().fit_transform(X=y_val_pred_std.reshape(-1, 1))
    # var_norm = y_val_pred_std.reshape(-1, 1)
    # prob = 1 / var_norm
    if is_norm:
      return var_norm.ravel()
    else:
      return y_val_pred_std.reshape(-1, 1)


  def best_score_(self, X=None, y=None):
    estimator, GridSearchCV = RegressionFactory.get_regression(
        method=self.name, 
        kernel=self.kernel, alpha=None, gamma=None, # # rbf, cosine
        search_param=self.search_param, X=X, y=y,  
        cv=self.cv, mt_kernel=self.mt_kernel) 
    r2, r2_std, mae, mae_std = CV_predict_score(
            estimator, X, y, n_folds=3, n_times=3, score_type='r2')
    return mae




class UncertainKNearestNeighbor(object):
  def __init__(self, 
        name,
        random_state=1, 
        cv=3, search_param=False,):
    self.name = name
    self.search_param = search_param
    self.cv = cv
    self.estimator = None
    self.random_state = random_state


  def fit(self, X_train, y_train, sample_weight=None):
    # # in fit function
    # # just return estimator with best param with X_train, y_train
    np.random.seed(self.random_state)
    n_features = X_train.shape[1]

    self.X_train = X_train
    self.y_train = y_train

    if self.estimator is None: # # for not always search parameters:
    # if self.estimator is None or self.search_param: # # either self.estimator is None or search_param is True-> require search
      estimator, GridSearchCV = RegressionFactory.get_regression(method="u_knn", 
          search_param=self.search_param, X=X_train, y=y_train,  
          cv=self.cv) 

      self.estimator = estimator
      self.GridSearchCV = GridSearchCV
    self.estimator.fit(X_train, y_train)
    return self.estimator

  def predict(self, X_val, get_variance=False):
    y_val_pred = self.estimator.predict(X_val)
    if get_variance:

      best_nb =  self.estimator.get_params()["n_neighbors"]
      nb_varlist = np.random.normal(loc=best_nb, scale=20, size=10)
      nb_varlist += abs(min(nb_varlist)) + 1

      y_val_preds = []
      for nb in nb_varlist:
        self.estimator.n_neighbors = int(nb)
        self.estimator.fit(self.X_train, self.y_train)
        y_val_preds.append(y_val_pred)

      return y_val_pred, np.mean(y_val_preds, axis=0)
    else:
      return y_val_pred

  def score(self, X_val, y_val):
    y_pred = self.predict(X_val, get_variance=False)    
    val_acc = metrics.r2_score(y_val, y_pred)
    return val_acc

  def predict_proba(self, X, is_norm=True):
    # # large variance -> probability to be observed small -> sorting descending take first
    # # small variance -> probability to be observed large 
    y_val_preds, y_val_pred_std = self.predict(X, get_variance=True)

    # # normalize variance to 0-1
    var_norm = MinMaxScaler().fit_transform(X=y_val_pred_std.reshape(-1, 1))
    if is_norm:
      return var_norm.ravel()
    else:
      return y_val_pred_std.reshape(-1, 1)


  def best_score_(self, X=None, y=None):
    estimator, GridSearchCV = RegressionFactory.get_regression(
        method="u_knn", alpha=None, gamma=None, # # rbf, cosine
        search_param=self.search_param, X=X, y=y,  
        cv=self.cv) 
    r2, r2_std, mae, mae_std = CV_predict_score(
            estimator, X, y, n_folds=3, n_times=3, score_type='r2')
    return mae





