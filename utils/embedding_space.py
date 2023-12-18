
from sklearn import neighbors
import metric_learn as mkl
import numpy as np
from sklearn.metrics import pairwise_distances
from plot import scatter_plot_2
import umap

class EmbeddingSpace(object):
	def __init__(self, embedding_method):
		self.embedding_method = embedding_method
		if embedding_method == "MLKR":
			learn_metric = mkl.MLKR(n_components=2, init="auto")
		elif embedding_method == "LFDA":
			learn_metric = mkl.LFDA(n_components=2, 
		 		k=20, embedding_type="plain") 
		 		# weighted, orthonormalized, plain
		elif embedding_method == "LMNN":
			learn_metric = mkl.LMNN(k=10, learn_rate=0.1,
				n_components=3) 
			# k: number of neighbors
		self.learn_metric = learn_metric
	# learn_metric = mkl.ITML_Supervised()
	# learn_metric = mkl.SDML_Supervised(sparsity_param=0.1, balance_param=0.0015,
	#           prior='covariance')

	def fit(self, X_train, y_train):
		n_features = X_train.shape[1]

		self.learn_metric.fit(X_train, y_train)
		X_train_embedded = self.learn_metric.transform(X_train)

		self.X_train = X_train
		self.y_train = y_train
		self.X_train_embedded = X_train_embedded

	def transform(self, X_val, get_min_dist=False):
		X_val_transform = self.learn_metric.transform(X_val)
		if get_min_dist:
			distances = pairwise_distances(X_val, self.X_train)
			max_distances = np.min(distances, axis=1)
			return X_val_transform, max_distances
		return X_val_transform

	# def predict_proba(self, X, is_norm=True):
	# 	# # large variance -> probability to be observed small -> sorting descending take first
	# 	# # small variance -> probability to be observed large 
	# 	X_val_transform, variance = self.transform(X, get_min_dist=True)

	# 	# # normalize variance to 0-1
	# 	var_norm = MinMaxScaler().fit_transform(X=variance.reshape(-1, 1))
	# 	if is_norm:
	# 	  return var_norm.ravel()
	# 	else:
	# 	  return variance.reshape(-1, 1)


class InversableEmbeddingSpace(object):
	def __init__(self, invs_emb_method="umap"):
		self.invs_emb_method = invs_emb_method
		if invs_emb_method == "umap":
			mapper = umap.UMAP(random_state=42)
		self.mapper = mapper


	def fit(self, X_train, y_train, save_at=None):

		self.mapper.fit(X_train)
		X_train_embedded = self.mapper.transform(X_train)

		self.X_train_embedded = X_train_embedded

	def transform(self, X_val, get_min_dist=False):
		X_val_transform = self.mapper.transform(X_val)
		if get_min_dist:
			distances = pairwise_distances(X_val, self.X_train)
			max_distances = np.min(distances, axis=1)
			return X_val_transform, max_distances
		return X_val_transform

	def invs_transform(self, X_test):
		inv_transformed_points = mapper.inverse_transform(X_test)














