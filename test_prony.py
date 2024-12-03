import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt

cpt.set_logging("INFO")

gf = cpt.Delhommeau()

def ref_function(kh):
    @np.vectorize
    def ff(x):
        return gf.fortran_core.old_prony_decomposition.ff(x, kh*np.tanh(kh), kh) + 2
    return ff

def reconstruction(a, lamda):
    def ff(x):
        return sum(a[i] * np.exp(lamda[i] * x) for i in range(len(a)))
    return ff

x = np.linspace(0.0, 10.0, 1000)
kh = 1.0
wavelength = 2*np.pi*1.0/kh
print(wavelength)
prony_1 = gf.find_best_exponential_decomposition(kh, method="fortran")
prony_2 = gf.find_best_exponential_decomposition(kh, method="python")

print(np.linalg.norm(ref_function(kh)(x) - reconstruction(*prony_1)(x))/np.linalg.norm(ref_function(kh)(x)))
print(np.linalg.norm(ref_function(kh)(x) - reconstruction(*prony_2)(x))/np.linalg.norm(ref_function(kh)(x)))

plt.plot(x, ref_function(kh)(x), label="ref", linestyle="--")
plt.plot(x, reconstruction(*prony_1)(x), label="fortran")
plt.plot(x, reconstruction(*prony_2)(x), label="python")
plt.legend()
plt.show()
