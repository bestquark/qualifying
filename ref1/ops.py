from qiskit.quantum_info import SparsePauliOp
import numpy as np

def omegak(lam, k, n):
    return np.sqrt((lam - np.cos(2*np.pi*k/n))**2 + np.sin(2*np.pi*k/n)**2)


def bangle(lam, n, k):
    num = np.sin(2*np.pi*k/n)
    den = lam - np.cos(2*np.pi*k/n)
    return np.arctan2(num, den)

def ci_jw(i, n):
    if i == 0:
        raise ValueError("i must be greater than 0")
    if i == 1:
        ham_list = []
        ham_list.append(("X", [0], 0.5))
        ham_list.append(("Y", [0], -0.5j))
        return SparsePauliOp.from_sparse_list(ham_list, num_qubits=n)
    if i>n:
        return ci_jw(i-n, n)
    ham_list = []
    ham_list.append(("Z"*(i-1) + "X", range(i), 0.5))
    ham_list.append(("Z"*(i-1) + "Y", range(i), -0.5j))
    return SparsePauliOp.from_sparse_list(ham_list, num_qubits=n)

def ci_jw_dagger(i, n):
    if i == 0:
        raise ValueError("i must be greater than 0")
    if i == 1:
        ham_list = []
        ham_list.append(("X", [0], 0.5))
        ham_list.append(("Y", [0], 0.5j))
        return SparsePauliOp.from_sparse_list(ham_list, num_qubits=n)
    if i>n:
        return ci_jw_dagger(i-n, n)
    ham_list = []
    ham_list.append(("Z"*(i-1) + "X", range(i), 0.5))
    ham_list.append(("Z"*(i-1) + "Y", range(i), 0.5j))
    return SparsePauliOp.from_sparse_list(ham_list, num_qubits=n)

def bk_ft(k, n):
    # k in -n/2+1, ..., n/2
    op = 0
    for i in range(1, n+1):
        op += ci_jw(i, n)*np.exp(2j*np.pi*k*(i-1)/n)
    return op/np.sqrt(n)

def bk_ft_dagger(k, n):
    # k in -n/2+1, ..., n/2
    op = 0
    for i in range(1, n+1):
        op += ci_jw_dagger(i, n)*np.exp(-2j*np.pi*k*(i-1)/n)
    return op/np.sqrt(n)

def ak_bog(lam, k, n):
    bang = bangle(lam, n, k)
    return np.cos(bang/2)*bk_ft(k, n) - 1j * np.sin(bang/2)*bk_ft_dagger(-k, n)

def ak_bog_dagger(lam, k, n):
    bang = bangle(lam, n, k)
    return np.cos(bang/2)*bk_ft_dagger(k, n) + 1j * np.sin(bang/2)*bk_ft(-k, n)

def jw_hamiltonian(lam, n):
    first_term = 0
    for i in range(1, n+1):
        first_term += -(ci_jw_dagger(i,n) - ci_jw(i,n)) @ (ci_jw(i+1,n) + ci_jw_dagger(i+1,n))
    # first_term *= n*0.5

    second_term = sum([2*ci_jw_dagger(i, n)@ci_jw(i, n) - SparsePauliOp("I"*n) for i in range(1, n+1)])*lam

    return first_term + second_term

def ftjw_hamiltonian(lam, n):
    # not sure if this is ok
    ham = 0
    for k in range(-n//2+1, n//2+1):
        if k>0 and k!=0 and k<n//2:
            # ham += 2*(lam - np.cos(2*np.pi*k/n))*bk_ft_dagger(k, n)@bk_ft(k, n) + 1j*np.sin(2*np.pi*k/n)*(bk_ft_dagger(-k,n)@bk_ft_dagger(k,n) - bk_ft(-k,n)@bk_ft(k,n))
            ham +=  ( (lam - np.cos(2*np.pi*k/n))*(bk_ft_dagger(k,n)@bk_ft(k,n) + bk_ft_dagger(-k, n)@bk_ft(-k,n)) 
                    + np.sin(2*np.pi*k/n)*(bk_ft_dagger(k,n)@bk_ft_dagger(-k,n) + bk_ft(-k,n)@bk_ft(k,n)) 
                    )
        elif k == 0:
            ham += 2*(lam - 1)*bk_ft_dagger(0, n)@bk_ft(0, n)
        elif k == n//2:
            ham += 2*(lam + 1)*bk_ft_dagger(n//2, n)@bk_ft(n//2, n)
        
    return ham - SparsePauliOp("I"*n)*n*lam

def bogftjw_hamiltonian(lam, n):
    ham = 0
    for k in range(-n//2+1, n//2+1):
        ham += 2*omegak(lam, k, n)* ak_bog_dagger(lam, k, n)@ak_bog(lam, k, n) 
    return ham