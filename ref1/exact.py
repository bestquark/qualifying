import numpy as np
from qiskit.circuit import QuantumCircuit

from . import fqft
from . import bogo


def omegak(lam, k, n):
    return np.sqrt((lam - np.cos(2*np.pi*k/n))**2 + np.sin(2*np.pi*k/n)**2)


def diag_evol(n, t, lam, order=None):
    circ = QuantumCircuit(n)
    if order is not None:
        # q2ind = {i: order[i] if order[i] <= n//2 else n//2 - order[i] for i in range(n)}
        # q2ind = {i: order[i] if order[i] <= n//2 else order[i] - n for i in range(n)}
        # print("evol order:", q2ind)
        # q2ind = {i: order[i] - n//2 + 1 for i in range(n)}
        # q2ind = {i: i-n//2+1 for i in range(n)}

        ## eg. case 1
        # q2ind = {i: order[i] - n//2 + 1 for i in range(n)}

        ## eg. case 3
        q2ind = {i: order[i] if order[i] <= n//2 else order[i] - n for i in range(n)}
        print("evol order:", q2ind)
    # else:
    #     fftswaps = fqft.recursive_fft_swaps(n)
    #     q2ind = {i: fftswaps[i] for i in range(n)}
    
    for qubit in range(n):
        # print("evolving qubit ", qubit, "(ind:", q2ind[qubit], ") with t=", t)
        circ.rz(2*omegak(lam, q2ind[qubit], n)*t, qubit) 

    return circ

def exact_evolution(hamiltonian, time=1, reps=20, lam=0.001, basis='ising'):
    n = hamiltonian.num_qubits
    assert np.log2(n).is_integer(), "Number of qubits must be a power of 2"

    # order = fqft.recursive_fft_swaps(n)
    ft = fqft.fqft_cirq(n)
    # ft = fqft.fqft(n)
    # order = order[n//2:] + order[:n//2]
    # order = [i + 1 for i in order]
    order = list(range(n))
    # m = n//2
    # order = [k - m if k<m else k - m + 1 for k in list(range(n))][::-1]
    # order = [k - n if k > n//2 else k for k in list(range(n))]
    print("order:", order)
    bg = bogo.bogo(n, lam, order)
    de = diag_evol(n, time, lam, order)

    qc = QuantumCircuit(n)
    qc = qc.compose(ft)
    qc = qc.compose(bg)
    
    qc_dag = qc.inverse()

    evol = QuantumCircuit(n)

    if basis == 'ising':
        evol = evol.compose(qc)
        evol.barrier()
    elif basis == 'diagonal':
        pass
    else:
        raise ValueError("Invalid basis")

    # not sure why this helps with getting the expected sigma_z correct when lam < 1
    # if lam < 1:
    #     evol.x(n//2)

    evol.barrier()
    evol = evol.compose(de)
    evol = evol.compose(qc_dag)

    return evol
