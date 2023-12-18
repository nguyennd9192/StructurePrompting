import random
from itertools import combinations
from scipy.special import comb 

class CombinationGeneratorFactory(object):

    @staticmethod
    def get_generator(method, items,
            sample_size, n_shuffle,
            reverse=True, extened=False):
        
        method = method.strip().lower()

        if method == "bagging":
            return BaggingCombinationGenerator(
                items=items, sample_size=sample_size, n_shuffle=n_shuffle
            )
        elif method == "k-combination":
            return KCombinationGenerator(
                items=items, sample_size=sample_size, n_shuffle=n_shuffle
            )
            # n_combinations = int(comb(len(items), sample_size))
            # if n_combinations > n_shuffle:
            #     return KCombinationGenerator(
            #         items=items, sample_size=sample_size, n_shuffle=n_shuffle
            #     )
            # else:
            #     return FullKCombinationGenerator(
            #         items=items, sample_size=sample_size
            #     )
        elif method == "all":
            return FCGAll(
                items=items, sample_size=sample_size
            )

class BaggingCombinationGenerator(object):

    def __init__(self, items, sample_size, n_shuffle):
        self.__items = items
        self.__sample_size = sample_size
        self.__n_shuffle = n_shuffle

    def __iter__(self):
        bag_size = int(self.__sample_size * len(self.__items))
        for i in range(self.__n_shuffle):
            sample = random.sample(self.__items, bag_size)
            yield sample

    def __len__(self):
        return self.__n_shuffle
    
class KCombinationGenerator(object):
    
    def __init__(self, items, sample_size, n_shuffle):
        self.__items = items
        self.__sample_size = sample_size
        self.__n_shuffle = n_shuffle
    
    def __iter__(self):
        for i in range(self.__n_shuffle):
            sample = random.sample(self.__items, self.__sample_size)
            yield sample
    
    def __len__(self):
        return self.__n_shuffle

class FullKCombinationGenerator(object):
    
    def __init__(self, items, sample_size=10):
        self.__items = items
        self.__sample_size = sample_size
    
    def __iter__(self):
        for sample in combinations(self.__items, self.__sample_size):
            yield list(sample)
    
    def __len__(self):
        return int(comb(len(self.__items), self.__sample_size))

class FCGAll(object):
    """
    Generate all possible combinations of features. It means 2^N combinations.
    """

    def __init__(self, items, sample_size=None):
        self.__items = items
        self.__sample_size = sample_size

    def __iter__(self):
        if self.__sample_size is None:
            i_range = range(1, len(self.__items) + 1)

            for i in i_range:
                for c in combinations(self.__items, i):
                    yield c
        else:
            i_range = range(1, self.__sample_size + 1)

            for i in i_range:
                for c in combinations(self.__items, i):
                    yield c

    def __len__(self):
        return 2**len(self.__items) - 1
