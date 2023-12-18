
import torch
from torch import nn
from mixture_of_experts import MoE
import numpy as np
from pandas import read_csv
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torch.utils.data import random_split
from torch import Tensor
from torch.nn import Linear, Dropout
from torch.nn import ReLU
from torch.nn import Sigmoid
from torch.nn import Module

from torch.nn import MSELoss
from torch.optim import SGD
import torch
from numpy import vstack
from sklearn.metrics import r2_score, mean_absolute_error

from torch.nn.init import kaiming_uniform_
from torch.nn.init import xavier_uniform_
import random, copy
from sklearn.preprocessing import MinMaxScaler


from smt.applications import MOE
from smt.problems import LpNorm

# dataset definition
class CSVDataset(Dataset):
	# load the dataset
	def __init__(self, X, y):
		# store the inputs and outputs
		self.X = X
		self.y = y

	# number of rows in the dataset
	def __len__(self):
		return len(self.X)

	# get a row at an index
	def __getitem__(self, idx):
		return [self.X[idx], self.y[idx]]

	# get indexes for train and test rows
	def get_splits(self, n_test=0.33):
		# determine sizes
		test_size = round(n_test * len(self.X))
		train_size = len(self.X) - test_size
		# calculate the split
		return random_split(self, [train_size, test_size])



def enable_dropout(m):
	for each_module in m.modules():
		if each_module.__class__.__name__.startswith('Dropout'):
			# print ("Enable module:", each_module.__class__.__name__)
			each_module.train()
	return m


 

class NN_estimator():
	def __init__(self, 
		verbose=False, NN_kwargs=None):

		self.estimator = None
		self.NN_kwargs = NN_kwargs

		# calculate split
		# train, test = dataset.get_splits(n_test=0.33)


	def fit(self, X_train, y_train):
		# train the model
		# define the optimization
		train = CSVDataset(X_train.astype(float), y_train.astype(float))
		train_dl = DataLoader(train, 
				batch_size=3, shuffle=True)
		
		if self.estimator == None:
			if self.NN_kwargs["method"] == "fully_connected":
				estimator = MLP(n_inputs=X_train.shape[1])
			elif self.NN_kwargs["method"] == "moe":
				estimator = MixtureOfExperts(
					n_inputs=X_train.shape[1], kwargs=self.NN_kwargs)
			elif self.NN_kwargs["method"] == "LeNet":
				estimator = LeNet(
					n_inputs=X_train.shape[1], kwargs=self.NN_kwargs)


			self.estimator = estimator

		model = copy.copy(self.estimator)
		criterion = MSELoss()
		optimizer = SGD(model.parameters(), 
			lr=self.NN_kwargs["lr"], momentum=self.NN_kwargs["momentum"])
		# enumerate epochs
		for epoch in range(self.NN_kwargs["n_epoches"]):
			# enumerate mini batches
			for i, (inputs, targets) in enumerate(train_dl):
				# clear the gradients
				optimizer.zero_grad()
				# compute the model output
				yhat = model(inputs.float())

				# calculate loss
				yhat = torch.reshape(yhat, targets.shape)

				loss = criterion(yhat, targets.float())
				# credit assignment
				loss.backward()
				# update model weights
				optimizer.step()
		self.estimator = copy.copy(model)


	# evaluate the model
	def evaluate_model(self, test_dl, model):
		predictions, actuals = list(), list()
		for i, (inputs, targets) in enumerate(test_dl):
			# evaluate the model on the test set
			yhat = model(inputs)
			# retrieve numpy array
			yhat = yhat.detach().numpy()
			actual = targets.numpy()
			actual = actual.reshape((len(actual), 1))
			# round to class values
			yhat = yhat.round()
			# store
			predictions.append(yhat)
			actuals.append(actual)
		predictions, actuals = vstack(predictions).ravel(), vstack(actuals).ravel()
		# print ("predictions.shape:", predictions.shape)
		error = mean_absolute_error(predictions, actuals)
		# calculate accuracy
		# acc = accuracy_score(actuals, predictions)
		return error, predictions, actuals

	def predict(self, X_test, get_variance=False):
		test = CSVDataset(X_test.astype(float), np.array([0.0]*X_test.shape[0]))
		test_dl = DataLoader(test, batch_size=32, shuffle=True)

		model = copy.copy(self.estimator)
		predictions = list()
		# print ("Here")
		for i, (inputs, targets) in enumerate(test_dl):
			# evaluate the model on the test set
			yhat = model(inputs.float())
			# retrieve numpy array
			yhat = yhat.detach().numpy().ravel()
			# print  ("Batch", i, yhat)

			# round to class values
			# yhat = yhat.round()
			# store
			predictions.append(yhat)
		return np.concatenate(predictions)

	def predict_proba(self, X_test):
		T = 10 
		# enable_dropout(model)
		all_preds = []
		all_errors = []

		# predict stochastic dropout model T times
		# model = self.estimator
		self.estimator.is_drop = True

		for t in range(T):
			# model = enable_dropout(model) # STILL NOT WORKING WITH DROPOUT AT TEST-TIME
			preds = self.predict(X_test)
			# errors = np.abs(preds - actuals)

			all_preds.append(preds) 
			# print (preds)
			# all_errors.append(errors)
		   
		# mean prediction
		var = np.var(all_preds, axis=0) # 
		mean = np.mean(all_preds, axis=0)

		var_norm = MinMaxScaler().fit_transform(X=var.reshape(-1, 1))

		return var_norm.ravel()

	def score(self, X_val, y_val):
		y_pred = self.predict(X_val, get_variance=False)    
		val_acc = r2_score(y_val, y_pred)
		return val_acc

	def best_score_(self):
		# # has not implemented yet
		return None

	def get_params(self):
		# # has not implemented yet
		return self.NN_kwargs




# model definition
class MLP(Module):
	# define model elements
	def __init__(self, n_inputs):
		super(MLP, self).__init__()
		# input to first hidden layer
		hyperparams = {'l1_out': 5,  'l2_out': 5,
				  'l1_drop': 0.1,  'l2_drop': 0.1,
				  'batch_size': 32,  'epochs': 10}
		self.hidden1 = Linear(n_inputs, 10)
		kaiming_uniform_(self.hidden1.weight, nonlinearity='relu')
		self.act1 = ReLU()
		self.dropout1 = Dropout(hyperparams['l1_drop'])

		# second hidden layer
		self.hidden2 = Linear(10, 8)
		kaiming_uniform_(self.hidden2.weight, nonlinearity='relu')
		self.act2 = ReLU()
		self.dropout2 = Dropout(hyperparams['l2_drop'])

		# third hidden layer and output
		self.hidden3 = Linear(8, 4)
		xavier_uniform_(self.hidden3.weight)
		self.act3 = Linear(4, 1)
		self.is_drop = False

 
	# forward propagate input
	def forward(self, X):
		# input to first hidden layer
		if self.is_drop:
			X = self.hidden1(X.float())
			X = self.act1(X)

			self.dropout1 = Dropout(random.uniform(0.5, 1))
			X = self.dropout1(X)

			 # second hidden layer
			X = self.hidden2(X)
			X = self.act2(X)
			self.dropout2 = Dropout(random.uniform(0.5, 1))
			X = self.dropout2(X)

			# third hidden layer and output
			X = self.hidden3(X)
			X = self.act3(X)
			self.is_drop = False

			return X
		else:
			X = self.hidden1(X.float())
			X = self.act1(X)

			# second hidden layer
			X = self.hidden2(X)
			X = self.act2(X)

			# third hidden layer and output
			X = self.hidden3(X)
			X = self.act3(X)
			self.is_drop = False

			return X


class Flatten(nn.Module):
	def __init__(self):
		super(Flatten, self).__init__()

	def forward(self, x):
		x = x.view(x.size(0), -1)
		return x


class LeNet(nn.Module):
	def __init__(self, n_inputs, kwargs, droprate=0.5):
		super(LeNet, self).__init__()
		self.model = nn.Sequential()
		self.model.add_module('conv1', nn.Conv1d(n_inputs, 20, kernel_size=5, padding=2))
		self.model.add_module('dropout1', nn.Dropout(p=droprate))
		self.model.add_module('maxpool1', nn.MaxPool1d(20, stride=2))
		# self.model.add_module('conv2', nn.Conv1d(20, 50, kernel_size=5, padding=2))
		# self.model.add_module('dropout2', nn.Dropout(p=droprate))
		# self.model.add_module('maxpool2', nn.MaxPool1d(2, stride=2))
		# self.model.add_module('flatten', Flatten())
		# self.model.add_module('dense3', nn.Linear(50*7*7, 500))
		# self.model.add_module('relu3', nn.ReLU())
		# self.model.add_module('dropout3', nn.Dropout(p=droprate))
		self.model.add_module('final', nn.Linear(50, 1))
		
	def forward(self, x):
		return self.model(x)



class MixtureOfExperts(Module):
	def __init__(self,	n_inputs=1, kwargs=None):
		super(MixtureOfExperts, self).__init__()

		self.kwargs = kwargs
		if self.kwargs["hidden_dim"] == "default":
			self.hidden_dim = kwargs["num_experts"] * 4
		else:
			self.hidden_dim = kwargs["hidden_dim"]

		self.linear1 = Linear(n_inputs, n_inputs)
		self.activation = Sigmoid()
		
		self.moe2 = MoE(
				dim = n_inputs,
				num_experts = self.kwargs["num_experts"],               # increase the experts (# parameters) of your model without increasing computation
				hidden_dim = self.hidden_dim,           # size of hidden dimension in each expert, defaults to 4 * dimension
				activation = nn.LeakyReLU,      # use your preferred activation, will default to GELU
				second_policy_train = 'random', # in top_2 gating, policy for whether to use a second-place expert
				second_policy_eval = 'random',  # all (always) | none (never) | threshold (if gate value > the given threshold) | random (if gate value > threshold * random_uniform(0, 1))
				second_threshold_train = 0.2,
				second_threshold_eval = 0.2,
				capacity_factor_train = 1.25,   # experts have fixed capacity per batch. we need some extra capacity in case gating is not perfectly balanced.
				capacity_factor_eval = 2.,      # capacity_factor_* should be set to a value >=1
				loss_coef = 1e-2                # multiplier on the auxiliary expert balancing auxiliary loss
			)
		
		# third hidden layer and output
		self.hidden3 = Linear(n_inputs, 4)
		xavier_uniform_(self.hidden3.weight)
		self.act3 = Linear(4, 1)

		self.is_drop = False


	def forward(self, X):
		# input to first hidden layer
		if self.is_drop:
			X = self.linear1(X.float())

			self.dropout1 = Dropout(random.uniform(0.5, 1))
			X = self.dropout1(X)

			X, loss = self.moe2(X)
			self.dropout2 = Dropout(random.uniform(0.5, 1))
			X = self.dropout2(X)

			# third hidden layer and output
			X = self.hidden3(X)
			X = self.act3(X)
			return X
		else:
			X = self.linear1(X.float())

			X, loss = self.moe2(X)

			# third hidden layer and output
			X = self.hidden3(X)
			X = self.act3(X)
			return X
	
class MixtureOfExperts_mpt(object):
	def __init__(self, 
		random_state=1, cv=3, n_times=3, 
		search_param=False, verbose=False):

		self.search_param = search_param
		self.kernel = 'rbf'
		self.verbose = verbose
		self.cv = cv
		self.n_times = n_times
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
			moe = MOE(smooth_recombination=True, n_clusters=2)
			moe.experts.remove("RMTC")
			print (moe.experts)
			moe.set_training_values(X_train, y_train)
			moe.train()

			self.prob = LpNorm(ndim=X_train.shape[1])

		 
			self.estimator = moe

		return self.estimator

	def predict(self, X_val, get_variance=False):
		y_val_pred = self.estimator.predict_values(X_val)
		y_val_pred_std = self.prob(xe)

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


	def best_score_(self):
		# # some conflict meaning between best_score_ for GridSearchCV object and this attribute:
		# # GridSearchCV.best_score_ returns cv score of best parameter
		# # this UncertainGaussianProcess.best_score_returns cv score of given params
		if self.GridSearchCV is None:
		  r2, r2_std, mae, mae_std = CV_predict_score(self.estimator, self.X_train, self.y_train, 
					n_folds=3, n_times=3, score_type='r2')
		  result = r2
		else:
		  result = self.GridSearchCV.best_score_
		return result



