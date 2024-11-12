from __future__ import annotations
import numpy as np

from qiskit.circuit import QuantumCircuit

def recursive_steps(lst, depth, max_depth):
    """
    Recursively separate the list into evens and odds and print each step.

    :param lst: Current list to be separated
    :param depth: Current depth level
    :param max_depth: Maximum depth to recurse
    """
    if depth == max_depth:
        return lst

    # Separate the list into even and odd indices
    evens = lst[::2]
    odds = lst[1::2]

    new_even = recursive_steps(evens, depth + 1, max_depth)
    new_odd = recursive_steps(odds, depth + 1, max_depth)

    return new_even + new_odd


def recursive_fft_swaps(n):
    """
    Perform recursive FFT-like swaps on a list of size 2^n.
    At each step, separate even indices to the first half and odd indices to the second half.
    Recursively do this for n-1 steps and print each iteration.

    :param n: The exponent such that the list size is 2^n
    """
    lst = list(range(n))
    indlist = recursive_steps(lst, 0, n - 1)

    return indlist


def fnk(circuit, q0, q1, k, n, dagger=False):
    q0, q1 = q1, q0
    sign = -1 if dagger else 1
    circuit.p(sign*2*np.pi*k/n, q0)
    circuit.cx(q0, q1)
    circuit.ch(q1, q0)
    circuit.cx(q0, q1)
    circuit.cz(q0, q1)
    

def fermionic_swap(circuit, qubit1, qubit2):
    # circuit.cx(qubit1,qubit2)
    # circuit.h(qubit1)
    # circuit.h(qubit2)
    # circuit.cx(qubit1,qubit2)
    # circuit.h(qubit1)
    # circuit.h(qubit2)
    # circuit.cx(qubit1,qubit2)

    # # this also works
    circuit.swap(qubit1, qubit2)
    circuit.cz(qubit1,qubit2)


def bit_reversed_indices(n):
    """Compute the bit-reversed index for each qubit."""
    bits = int(np.log2(n))
    reversed_indices = []
    for i in range(n):
        binary = format(i, f'0{bits}b')
        reversed_binary = binary[::-1]
        j = int(reversed_binary, 2)
        reversed_indices.append(j)
    return reversed_indices


def fqft_dagger(n, first_iter=True):
    circuit = QuantumCircuit(n)
    assert np.log2(n).is_integer(), "n must be a power of 2."

    if n == 2:
        # circuit.h(0)
        fnk(circuit, 0, 1, 0, 2)
    
    if n > 2:
        for i in range(n//2):
            print("doing fnk", i, n//2 + i, i, n)
            fnk(circuit, i, n//2 + i, i, n)
        
        subcirc = fqft_dagger(n//2, first_iter=False)

        circuit.compose(subcirc, qubits=list(range(n//2)), inplace=True)
        circuit.compose(subcirc, qubits=list(range(n//2, n)), inplace=True)

    # bit reversal
    if first_iter:
        reversed_indices = bit_reversed_indices(n)
        for i in range(n):
            j = reversed_indices[i]
            if i < j:
                fermionic_swap(circuit, i, j)

        # # apply phase gate
        # for i in range(n):
        #     circuit.p(2*np.pi*i/n, i)
            # circuit.h(i)

    return circuit

def fswap_cascade(circuit, n):
    m = n//2 - 1
    # make a pyramid of swaps
    offset = 1
    while m>0:
        print(m)
        for i in range(m):
            fermionic_swap(circuit, 2*i + offset, 2*i + offset + 1)
        m -= 1
        offset += 1

def inverse_fswap_cascade(circuit, n):
    m = 1
    offset = n//2 - 1
    while m < n//2:
        for i in range(m):
            fermionic_swap(circuit, 2*i + offset, 2*i + offset + 1)
        m += 1
        offset -= 1

def fqft_latorre(n, first_iter=True):
    circuit = QuantumCircuit(n)
    assert np.log2(n).is_integer(), "n must be a power of 2."

    if n == 2:
        fnk(circuit, 0, 1, 0, 2)

    if n > 2:
        inverse_fswap_cascade(circuit, n)
        for i in range(n//2):
            fnk(circuit, 2*i, 2*i + 1, i, n)
        fswap_cascade(circuit, n)
        
        subcirc = fqft_latorre(n//2, first_iter=False)

        circuit.compose(subcirc, qubits=list(range(n//2)), inplace=True)
        circuit.compose(subcirc, qubits=list(range(n//2, n)), inplace=True)
    
    if first_iter:
        for i in range(n//4):
            fermionic_swap(circuit, 4*i+1, 4*i+2)

    return circuit
        
def fqft(n, first_iter=True):
    circuit = fqft_dagger(n, first_iter=first_iter)
    return circuit.inverse()



def fqft_cirq(n):
    import cirq
    import openfermion
    qubits = [cirq.GridQubit(0, i) for i in range(n)]
    circuit = cirq.Circuit()
    circuit.append(openfermion.ffft(qubits))

    qasm_output = cirq.QasmOutput((circuit.all_operations()), qubits)
    qasm_circuit = QuantumCircuit().from_qasm_str(str(qasm_output))
    
    return qasm_circuit
