import pandas as pd
import numpy as np
import warnings
from joblib import parallel_backend
from sklearn.metrics import r2_score, mean_absolute_error, accuracy_score, precision_recall_fscore_support
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import MinMaxScaler, MaxAbsScaler, StandardScaler, scale
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel, ExpSineSquared
from sklearn.neighbors import KNeighborsRegressor
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning

class RegressionFactory(object): 
    
  @staticmethod
  @ignore_warnings(category=ConvergenceWarning)
  def get_regression(method, kernel='rbf', alpha=1, gamma=1, 
      search_param=False, X=None, y=None, cv=3, 
      mt_kernel=None):
    method = method.strip().lower()
    if method == "e_krr":
        if search_param:
          # alpha, gamma, scores_mean, scores_std = RegressionFactory.kernel_ridge_parameter_search(
          #       X=X, y_obs=y, kernel=kernel, n_folds=cv, n_times=n_times)
          model, md_selection = RegressionFactory.kernel_ridge_cv(X=X, y_obs=y, 
                      kernel=kernel, cv=10)
        else:
          model = KernelRidge(kernel=kernel, alpha=alpha, gamma=gamma)
        return model, md_selection

    elif method == "u_gp":
        if search_param:
          model, md_selection = RegressionFactory.gaussian_process_cv_with_noise(
              X=X, y_obs=y, cv=cv, mt_kernel=mt_kernel)
        else:          
          default_kernel = RegressionFactory.gp_kernel(c=1.0, l=0.5, n=0.2)
          model = GaussianProcessRegressor(alpha=0.1, kernel=default_kernel)
          md_selection = None
        model.fit(X, y)

    elif method == "u_knn":
        if search_param:
          model, md_selection = RegressionFactory.knn_cv(X=X, y_obs=y, cv=cv)
        else:
          model = KNeighborsRegressor(n_neighbors=10, metric="minkowski")
          md_selection = None
        model.fit(X, y)


    return model, md_selection

  @staticmethod
  def CV_predict(model, X, y, n_folds=3, n_times=3, is_gp=False):
  
      if (n_folds <= 0) or (n_folds > len(y)):
          n_folds = len(y)
          n_times = 1

      y_predicts = []
      for i in range(n_times):
          indexes = np.random.permutation(range(len(y)))

          kf = KFold(n_splits=n_folds)

          y_cv_predict = []
          cv_test_indexes = []
          cv_train_indexes = []
          for train, test in kf.split(indexes):
              # cv_train_indexes += list(indexes[train])
              # print(train, test)
              cv_test_indexes += list(indexes[test])

              X_train, X_test = X[indexes[train]], X[indexes[test]]
              y_train, Y_test = y[indexes[train]], y[indexes[test]]

              model.fit(X_train, y_train)

              # y_train_predict = model.predict(X_train)
              if is_gp:
                  y_test_predict, y_test_pred_prob = model.predict(X_test)
              else:
                  y_test_predict = model.predict(X_test)
              y_cv_predict += list(y_test_predict)

          cv_test_indexes = np.array(cv_test_indexes)
          # print(cv_test_indexes)
          rev_indexes = np.argsort(cv_test_indexes)

          y_cv_predict = np.array(y_cv_predict)

          y_predicts += [y_cv_predict[rev_indexes]]

      y_predicts = np.array(y_predicts)

      return y_predicts

  @staticmethod
  def kernel_ridge_parameter_search(X, y_obs, kernel='rbf',
                                n_folds=3, n_times=3):
      # parameter initialize
    gamma_log_lb = -2.0
    gamma_log_ub = 2.0
    alpha_log_lb = -4.0
    alpha_log_ub = 1.0
    n_steps = 10
    n_rounds = 4
    alpha = 1
    gamma = 1
    lb = 0.8
    ub = 1.2
    n_instance = len(y_obs)

    if (n_folds <= 0) or (n_folds > n_instance):
        n_folds = n_instance
        n_times = 1 

    # Start
    for i in range(n_rounds):
        scores_mean = []
        scores_std = []
        gammas = np.logspace(gamma_log_lb, gamma_log_ub, num=n_steps)
        for gamma in gammas:
            k_ridge = KernelRidge(alpha=alpha, gamma=gamma, kernel=kernel)
            y_predicts = RegressionFactory.CV_predict(
                k_ridge, X, y_obs, n_folds=n_folds, n_times=n_times)
            cv_scores = list(map(lambda y_predict: r2_score(y_obs, y_predict), y_predicts))

            scores_mean += [np.mean(cv_scores)]
            scores_std += [np.std(cv_scores)]

        best_index = np.argmax(scores_mean)
        gamma = gammas[best_index]
        gamma_log_lb = np.log10(gamma * lb)
        gamma_log_ub = np.log10(gamma * ub)
        scores_mean = []
        scores_std = []
        alphas = np.logspace(alpha_log_lb, alpha_log_ub, num=n_steps)
        for alpha in alphas:
            k_ridge = KernelRidge(alpha=alpha, gamma=gamma, kernel=kernel)
            y_predicts = RegressionFactory.CV_predict(
                k_ridge, X, y_obs, n_folds=n_folds, n_times=n_times)
            cv_scores = list(map(lambda y_predict: r2_score(y_obs, y_predict), y_predicts))

            scores_mean += [np.mean(cv_scores)]
            scores_std += [np.std(cv_scores)]

        best_index = np.argmax(scores_mean)
        alpha = alphas[best_index]
        alpha_log_lb = np.log10(alpha * lb)
        alpha_log_ub = np.log10(alpha * ub)

    return alpha, gamma, scores_mean[best_index], scores_std[best_index]

  def kernel_ridge_cv(X, y_obs, kernel, cv=10):
    n_steps = 10

    alpha_lb = -2
    alpha_ub = 2
    alphas = np.logspace(alpha_lb, alpha_ub, n_steps)

    gamma_lb = -2
    gamma_ub = 2
    gammas = np.logspace(gamma_lb, gamma_ub, n_steps)
    param_grid = {"alpha": alphas, "gamma": gammas}

    GridSearch = GridSearchCV(KernelRidge(kernel=kernel), param_grid=param_grid,
          cv=cv, n_jobs=-1, scoring="neg_mean_absolute_error") # # scoring
    GridSearch.fit(X, y_obs)
    best_model = GridSearch.best_estimator_
    return best_model, GridSearch


  def gp_kernel(c, l, n=0.05):
    tmp = ConstantKernel(constant_value=c)*RBF(length_scale=l, length_scale_bounds=(1e-5, 1e-1)) + WhiteKernel(noise_level=n)
    # tmp = ExpSineSquared(length_scale=l, periodicity=c)
    return tmp

  @staticmethod
  @ignore_warnings(category=ConvergenceWarning)
  def gaussian_process_cv_with_noise(X, y_obs, cv=10, mt_kernel=None):
    n_steps = 5
    rbf_length_lb = 1e-3# 1e-2 # -2
    rbf_length_ub = 1e1# 3 # 5 1
    # rbf_lengths = np.logspace(rbf_length_lb, rbf_length_ub, n_steps)
    rbf_lengths = np.linspace(rbf_length_lb, rbf_length_ub, n_steps)

    # const_lb = -1
    # const_ub = 2
    # # consts = np.logspace(const_lb, const_ub, n_steps)
    # consts = [1.0]

    noise_lb = 1e-5
    noise_ub = 1e-1
    # noises = np.logspace(noise_lb, noise_ub, n_steps)
    noises = np.linspace(noise_lb, noise_ub, n_steps)


    # noises = [0.05]

    alpha_lb = 1e-6 # -2 
    alpha_ub = 1e-3 # 1
    # alphas = np.logspace(alpha_lb, alpha_ub, n_steps)
    alphas = np.linspace(alpha_lb, alpha_ub, n_steps)

    if mt_kernel is None:
      # # grid search with kernel length and noise
      param_grid = {"alpha": alphas,
                    "kernel": [RegressionFactory.gp_kernel(l, n) 
                    for l in rbf_lengths
                    for n in noises]}
    else:
      # # grid search with kernel length only
      param_grid = {"alpha": alphas,
                    "kernel": [RegressionFactory.gp_kernel(1.0, l, mt_kernel)
                    for l in rbf_lengths]}
    GridSearch = GridSearchCV(GaussianProcessRegressor(),param_grid=param_grid,
                cv=cv, n_jobs=1, scoring="neg_mean_absolute_error") # # scoring
    # with parallel_backend('threading'):
    GridSearch.fit(X, y_obs)
    best_model = GridSearch.best_estimator_
    print (best_model.get_params())
    return best_model, GridSearch


  def knn_cv(X, y_obs, cv=10):

    metric_list = ["minkowski", "euclidean", "chebyshev"]

    n_steps = 10

    n_inst = X.shape[0]
    neighbor_lb = 1
    neighbor_ub = int(n_inst / 2)
    n_neighbors_list = range(neighbor_lb, neighbor_ub, int(neighbor_ub / n_steps))

    param_grid = {"n_neighbors": n_neighbors_list,
        "metric": metric_list}

    model = KNeighborsRegressor()

    GridSearch = GridSearchCV(model,param_grid=param_grid,
                cv=cv, n_jobs=-1, scoring="neg_mean_absolute_error")

    with parallel_backend('threading'):
      GridSearch.fit(X, y_obs)

    best_model = GridSearch.best_estimator_
    return best_model, GridSearch




def CV_predict_score(model, X, y, n_folds=3, n_times=3, score_type='r2'):

    if (n_folds <= 0) or (n_folds > len(y)):
        n_folds = len(y)
        n_times = 1

    y_predicts = []
    scores = []
    errors = []
    for i in range(n_times):
        y_predict = RegressionFactory.CV_predict(model=model, 
          X=X, y=y, n_folds=n_folds, n_times=1)
        # n_times = 1 then the result has only 1 y_pred array
        y_predict = y_predict[0]

        if score_type == "r2":
            this_score = r2_score(y_true=y, y_pred=y_predict)
            this_err = mean_absolute_error(y_true=y, y_pred=y_predict)
            errors.append(this_err)
            scores.append(this_score)

        if score_type == "clf-score":
            this_score = precision_recall_fscore_support(y_true=y, y_pred=y_predict, 
                average='macro')
            scores.append(this_score)
    if score_type == "r2":
        return np.mean(scores), np.std(scores), np.mean(errors), np.std(errors)
    
    if score_type == "clf-score":
        scores = np.array(scores)
        
        precisions = scores[:, 0]
        recalls = scores[:, 1]
        f1_scores = scores[:, 2]
        support = scores[0, 3]

        return_result = [np.mean(precisions), np.std(precisions), 
                        np.mean(recalls), np.std(recalls), 
                        np.mean(f1_scores), np.std(f1_scores), support]

        return  return_result


    