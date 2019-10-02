#import pathos
#import pathos.helpers as hp
#import ctypes
#import numpy as np


#def shared_array(shape):   
#    shared_array_base = hp.mp.Array(ctypes.py_object, shape[0]*shape[1], lock=False)
#    shared_array = np.ctypeslib.as_array(shared_array_base)
#    shared_array = shared_array.reshape(*shape)
#    return shared_array


#num_lines = (15, 1)
#array = shared_array(num_lines)


#def parallel_function(i):
#    array[i, :] = i
#    print(array)


#class Fun:
#    def __init__(self, i):
#        self.a = i


#if __name__ == '__main__':
#    pool = pathos.pools.ProcessPool(processes=4)
        
#    args = range(num_lines[0])

#    pool.imap(parallel_function, args)
#    pool.close()
#    pool.join()
    
#    print(array)



import numpy as np
import ctypes
import array
import multiprocessing as mp
import random
from multiprocessing import get_context


class J:

    def __init__(self, a):
        self.a = a


def init_shared(ncell):
    shared_array_base = mp.Array(ctypes.py_object, ncell,lock=False)
    return(shared_array_base)


def tonumpyarray(shared_array):
    nparray= np.ctypeslib.as_array(shared_array)
    return nparray


def init_parameters(**kwargs):
    params = dict()

    for key, value in kwargs.items():
        params[key] = value
    return params


def init_worker(shared_array_,parameters_):
    global shared_array
    global shared_parr
    global dim

    shared_array = tonumpyarray(shared_array_)
    shared_parr = tonumpyarray(parameters_['shared_parr'])

    dim = parameters_['dimensions']


def worker_fun(data):
    i = data[0]
    j = data[1]

    arr = tonumpyarray(shared_array)
    parr = tonumpyarray(shared_parr)

    arr.shape = dim

    random.seed(i)
    rint = random.randint(1,10)

    parr[i] = j

    arr[i,...] = arr[i,...] * rint

    print(parr)
    print(arr)


if __name__ == '__main__':
    nrows = 10
    ncols = 3

    shared_array = init_shared(nrows*ncols)
    shared_parr = init_shared(nrows)

    params = init_parameters(shared_parr=shared_parr, dimensions=(nrows,ncols))

    arr = tonumpyarray(shared_array)
    parr = tonumpyarray(params['shared_parr'])

    arr.shape = (nrows,ncols)

    arr[...] = np.random.randint(1, nrows, size=(nrows,ncols))

    def get_data():
        for i in range(5):
            yield i, J(i)

    pool = mp.Pool(processes=8,initializer = init_worker, initargs=(shared_array,params))
    pool.map(worker_fun, get_data())
    pool.close()
    pool.join()

    print('*' * 40)
    print(arr)
    print('-' * 30)
    print(parr)
    for p in parr:
        print(p)
