import pandas as pd
from sklearn import manifold
import numpy as np


import numpy as np
from sklearn import svm
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.metrics import pairwise_distances

from coranking import coranking_matrix
from coranking.metrics import trustworthiness, continuity, LCMC
from nose import tools as nose

class Preprocessing():
    def __init__(self, similarity_matrix='', ticklabels=''):
        self.similarity_matrix = similarity_matrix
        self.ticklabels = ticklabels

    def iso_map(self, n_neighbors=5, n_components=2, eigen_solver='auto', tol=0, max_iter=None, path_method='auto',
                neighbors_algorithm='auto', n_jobs=None, metric="minkowski", p=2):
        # Y = manifold.Isomap(n_neighbors, n_components, eigen_solver, tol, max_iter, path_method,
        #                     neighbors_algorithm, n_jobs).fit_transform(self.similarity_matrix)
        # return Y, self.ticklabels

        dim_reduc_method = manifold.Isomap(n_neighbors, n_components, eigen_solver, tol, max_iter, path_method,
                            neighbors_algorithm, n_jobs).fit(self.similarity_matrix)
        Y = dim_reduc_method.transform(self.similarity_matrix)

        return Y, self.ticklabels, dim_reduc_method

    def standard(self, n_neighbors=5, n_components=2, reg=0.001, eigen_solver='auto', tol=1e-06, max_iter=100,
                 method='standard', hessian_tol=0.0001, modified_tol=1e-12,
                 neighbors_algorithm='auto', random_state=None, n_jobs=None):
        Y = manifold.LocallyLinearEmbedding(n_neighbors, n_components, reg, eigen_solver, tol, max_iter, 'standard',
                                            hessian_tol, modified_tol, neighbors_algorithm,
                                            random_state, n_jobs).fit_transform(self.similarity_matrix)

        return Y, self.ticklabels

    def locallyLinearEmbedding(self, n_neighbors=5, n_components=2, reg=0.001, eigen_solver='auto', tol=1e-06, max_iter=100,
                               method='modified', hessian_tol=0.0001, modified_tol=1e-12,
                               neighbors_algorithm='auto', random_state=None, n_jobs=None):
        Y = manifold.LocallyLinearEmbedding(n_neighbors, n_components, reg, eigen_solver, tol, max_iter, 'modified',
                                            hessian_tol, modified_tol, neighbors_algorithm,
                                            random_state, n_jobs).fit_transform(self.similarity_matrix)
        return Y, self.ticklabels

    def hessianEigenmapping(self, n_neighbors=5, n_components=2, reg=0.001, eigen_solver='auto', tol=1e-06, max_iter=100,
                            method='hessian', hessian_tol=0.0001, modified_tol=1e-12,
                            neighbors_algorithm='auto', random_state=None, n_jobs=None):
        Y = manifold.LocallyLinearEmbedding(n_neighbors, n_components, reg, eigen_solver, tol, max_iter, 'hessian',
                                            hessian_tol, modified_tol, neighbors_algorithm,
                                            random_state, n_jobs).fit_transform(self.similarity_matrix)
        return Y, self.ticklabels

    def tsne(self, n_components=2, perplexity=30.0, early_exaggeration=12.0, learning_rate=200.0, n_iter=1000,
             n_iter_without_progress=300, min_grad_norm=1e-07, metric='euclidean', init='random',
             verbose=0, random_state=None, method='barnes_hut', angle=0.5, n_jobs=None):

        dim_reduc_method = manifold.TSNE(n_components=n_components, perplexity=perplexity, early_exaggeration=early_exaggeration,
                          learning_rate=learning_rate, n_iter=n_iter, n_iter_without_progress=n_iter_without_progress,
                          min_grad_norm=min_grad_norm, metric=metric, init=init,
                          verbose=verbose, random_state=random_state, method=method, angle=angle)
        Y = dim_reduc_method.fit_transform(self.similarity_matrix)


        high_data = pairwise_distances(X=self.similarity_matrix, metric=metric)
        low_data = pairwise_distances(X=Y, metric=metric)

        results = get_ranking(high_data, low_data)

        # return Y, self.ticklabels
        return Y, self.ticklabels, dim_reduc_method, results


    def mds(self, n_components=2, metric=True, n_init=4, max_iter=300, verbose=0,
            eps=0.001, n_jobs=None, random_state=None, dissimilarity='euclidean'):
        Y = manifold.MDS(n_components, metric, n_init, max_iter, verbose, eps,
                         n_jobs, random_state, dissimilarity).fit_transform(self.similarity_matrix)
        return Y, self.ticklabels

    def get_all_preprocess(self):

        return ['iso_map', 'standard', 'locallyLinearEmbedding',
                'hessianEigenmapping', 'tsne', 'mds']


def get_ranking(high_data, low_data):
  # # https://www.sciencedirect.com/science/article/pii/S0925231209000101
  Q = coranking_matrix(high_data, low_data)

  results = dict({})
  c = continuity(Q)
  results["continuity"] = c

  low_data = LCMC(Q)
  results["LCMC"] = c

  t = trustworthiness(Q)
  results["trustworthiness"] = t

  return results
  # print nose.assert_equal(c, 1.)
  # t = trustworthiness(Q.astype(np.int64), min_k=2)






class OutlierRanking():
  def __init__(self, algo, configs):
    self.algo = algo
    self.configs = configs
    if self.algo == "Elliptic-Envelope":
      tmp = self.EllipticEnvelope(**self.configs)

    elif self.algo == "One-Class-SVM": 
      tmp = self.OneClassSVM(**self.configs)

    elif self.algo == "Local-Outlier-Factor":
      tmp = self.LocalOutlierFactore(**self.configs)

    elif self.algo == "Isolation-Forest":
      tmp = self.IsolationForest(**self.configs)

    self.algorithm = tmp

  def EllipticEnvelope(self, outliers_fraction=0.1):
    return EllipticEnvelope(contamination=outliers_fraction)#.fit(self.X).predict(self.X)

  def OneClassSVM(self, outliers_fraction=0.1, kernel="rbf", gamma=0.1):
    return svm.OneClassSVM(nu=outliers_fraction,
      kernel=kernel, gamma=gamma)#.fit(self.X).predict(self.X)

  def LocalOutlierFactore(self, outliers_fraction=0.1, n_neighbors=10):
    return LocalOutlierFactor(n_neighbors=35, contamination=outliers_fraction)#.fit_predict(self.X)

  def IsolationForest(self, outliers_fraction=0.1, random_state=42):
    return IsolationForest(contamination=outliers_fraction, 
      random_state=random_state)#.fit(self.X).predict(self.X)   

  def rank(self, X):

    self.fit(X_train=X)
    if self.algo == "Local-Outlier-Factor":
      y_pred = self.algorithm.fit_predict(X)
    else:
      y_pred = self.algorithm.fit(X).predict(X)
    
    return y_pred

  def fit(self, X_train):
    self.algorithm.fit(X_train)

  def predict(self, X_test):
    return self.algorithm.predict(X_test)

