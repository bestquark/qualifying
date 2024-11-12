import numpy as np

from qiskit.circuit.quantumcircuit import QuantumCircuit


def bgate(circuit, angle, q0, q1):
    """bogoliubov gate"""
    circuit.x(q1)
    circuit.cx(q1, q0)
    circuit.crx(angle, q0, q1)
    circuit.cx(q1, q0)
    circuit.x(q1)

def bgate_dagger(circuit, angle, q0, q1): 
    """bogoliubov gate"""
    circuit.x(q1)
    circuit.cx(q1, q0)
    circuit.crx(-angle, q0, q1)
    circuit.cx(q1, q0)
    circuit.x(q1)

# wrap from 0,pi to -pi/2,pi/2
def wrap_angle(angle):
    if angle > np.pi/2:
        return angle - np.pi/2
    elif angle < -np.pi/2:
        return angle + np.pi/2
    else:
        return angle

# def bangle(lam, n, k):
#     """angle for bogoilubov gate"""
#     num = lam - np.cos(2*np.pi*k/n)
#     den = np.sqrt((lam - np.cos(2*np.pi*k/n))**2 + np.sin(2*np.pi*k/n)**2)
#     # if lam < 1:
#     #     print("wangle: ", np.arccos(num/den), wrap_angle(np.arccos(num/den)))
#     #     return wrap_angle(np.arccos(num/den))
#     # print("angle: ", np.arccos(num/den))
#     return np.arccos(num/den) 

# def bangle(lam, n, k):
#     num = np.sin(2*np.pi*k/n)
#     den = lam - np.cos(2*np.pi*k/n)
#     print("num, den:", num, den.round(4))
#     if lam < 1:
#         print("wangle: ", np.arctan2(num, den), wrap_angle(np.arctan2(num, den)))
#         return wrap_angle(np.arctan2(num,den))
#     return np.arctan2(num,den)

def bangle(lam, n, k):
    num = np.sin(2*np.pi*k/n)
    den = lam - np.cos(2*np.pi*k/n)
    return np.arctan2(num, den)

def bogo(n, lam, order=None):
    assert np.log2(n).is_integer(), "n must be a power of 2"
    
    circuit = QuantumCircuit(n)

    m = n//2

    if order is not None:
        # eg1 and eg 3 work for time evol when initial state is |0>^n and |1>^n
        
        ## eg 1: symmetric ##
        # order = [k - m + 1 for k in order]
        #####

        ## eg 3: shift by n
        order = [k - n if k > n//2 else k for k in order]
        #####
        
        # for k>m, let k = k-m
        # order = [k if k <= m else m-k for k in order]
        # order = [k if k<=m else k-n for k in order]
        # order = [0,-1,1,2]
        # order = order[m:] + order[:m]  
        # order = order[::-1] 
        # order = [k - m if k<m else k - m + 1 for k in order]
        # print("new order:", order)

        # apply bgate between modes k and -k
        for k in range(m):
            
            ## eg 1: symmetric ##
            k = k if k == 0 else k
            q1 = order.index(k)
            q2 = order.index(m) if k == 0 else order.index(-k)
            #####


            ## eg 2: sandwich ##
            # q1 = order.index(m-k-1)
            # q2 = order.index(m+k)
            #####

            # if k%2 == 0:
            #     q1, q2 = q2, q1
            # k = 0 if k == m else k
            # print("lam", lam, "n", n, "k", k)

            
            
            # momentum = m-k if k != 0 else 0
            # angle = bangle(lam, n, momentum)
            # q1 = order.index(k - m + 1)
            # q2 = order.index(m - k)

            angle = bangle(lam, n, k)
            print(q1, q2, "angle:", angle)
            bgate(circuit, angle, q1, q2)

    return circuit

