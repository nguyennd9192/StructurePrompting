
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

def gp_kernel(c, l, n):
    kernel = ConstantKernel(constant_value=c)*RBF(length_scale=l) + WhiteKernel(noise_level=n)
    return kernel
