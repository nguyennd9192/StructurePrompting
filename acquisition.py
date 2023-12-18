from general_lib_37 import *
from sklearn.svm import SVC
def acquisition(X, y): 

	# # negative as "-1" label, positive as "+1" label
	labels = np.where(y<0, -1, 1)
	svm = SVC(C=1.0, kernel='rbf', gamma='scale', probability=True)
	svm.fit(X, labels)
	proba = svm.predict_proba(X)

	distance = np.abs(proba[:, 0] - proba[:, 1])
	boundaries = np.argsort(distance)
	n_queries = int(0.1*len(y))
	query_indexes = boundaries[:n_queries]

	outstanding_indexes = np.argsort(y)[:n_queries]
	return query_indexes, outstanding_indexes, svm

if __name__ == '__main__':

	X = np.array([[0.35461603,0.64538397],
			 [0.10610817, 0.89389183],
			 [0.65146049, 0.34853951],
			 [0.1, 0.9],
			 [0.2, 0.8],
			 [0.3, 0.7],
			 [0.4, 0.6],
			 [0.5, 0.5],
			 [0.6, 0.4],


			 ])
	distance = np.abs(X[:, 0] - X[:, 1] )
	print (distance)

	boundaries = np.argsort(distance)[:10]

	print (boundaries)