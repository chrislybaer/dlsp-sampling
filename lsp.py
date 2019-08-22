import numpy as np
import networkx as nx
from matplotlib import pyplot as plt

def cartesian_product(*arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    shape = [len(a) for a in arrays]
    for i, a in enumerate(np.ix_(*arrays)):
        arr[..., i] = a
    return arr

def pairwise_closure(generators, meet):
    generators = list(generators)
    intersections = set()
    for idx1 in range(len(generators)-1):
        for idx2 in range(idx1+1, len(generators)):
            intersections.add(meet(generators[idx1], generators[idx2]))
    return intersections

def close_under_meet(generators, meet):
    closed = set(generators)
    change = True
    while change:
        len_old = len(closed)
        closed |= pairwise_closure(closed, meet)
        len_new = len(closed)
        if len_old == len_new:
            change = False
    return closed

def compute_partial_order(lattice, meet):
    P = np.zeros((len(lattice), len(lattice)))
    for idx1, l1 in enumerate(lattice):
        for idx2, l2 in enumerate(lattice):
            P[idx1, idx2] = (meet(l1, l2) == l2)
    return P


def fdsft3(signal):
    N = len(signal)
    h = 1
    transform = signal.copy()
    while h < N:
        for i in range(0, N, 2*h):
            for j in range(i, i+h):
                x = transform[j]
                y = transform[j + h]
                transform[j] = x
                transform[j + h] = x - y
        h *= 2
    return transform

def popcount(arr):
    N = arr.shape[0]
    out = np.asarray([pypopcount(A) for A in arr])
    return out

def pypopcount(n):
    """ this is actually faster """
    return bin(n).count('1')

def fidsft3(transform):
    return fdsft3(transform)

def invert_po(P):
    x_leq_y = P.T
    M = np.eye(len(P))
    for row in range(len(P)):
        for col in range(row+1, len(P)):
            M[row, col] = -np.sum(M[row][np.where((x_leq_y[:, col] * x_leq_y[row, :]))[0]])
    return M.T


class LSP:
    def __init__(self, generators, meet, closed=False, use_la_inv=True):
        self.meet = meet
        self.lattice = None
        self.F_inv = None
        self.F = None
        self.use_la_inv = use_la_inv
        self.__build__(generators, closed)

    def __build__(self, generators, closed):
        print('closing under meet...')
        if not closed:
            lattice = close_under_meet(generators, self.meet)
        else:
            lattice = generators
        print('computing partial order...')
        iF = compute_partial_order(lattice, self.meet)
        G = nx.DiGraph()
        G.add_nodes_from(list(range(len(iF))))
        G.add_edges_from([(a, b) for a, b in zip(*np.where(iF - np.eye(len(iF))))])
        index = np.asarray(list(nx.topological_sort(G)))[::-1]
        self.F_inv = iF[index][:, index]
        print('inverting partial order matrix...')
        if self.use_la_inv:
            self.F = np.linalg.inv(self.F_inv)
        else:
            self.F = invert_po(self.F_inv)
        self.lattice = np.asarray(list(lattice))[index]

    def ft(self, s):
        return np.dot(self.F, s)

    def ift(self, s):
        return np.dot(self.F_inv, s)

    def fr(self, h):
        return np.dot(self.F_inv.T, h)

    def conv(self, h, s):
        return self.ift(self.fr(h) * self.ft(s))

    def sample(self, s, freq_idxs):
        """F_inv is lower triangular, therefore to sample frequencies, we just take frequencies,
           maybe make argument s a function in the future"""
        print('solving ...')
        plt.imshow(self.F_inv[freq_idxs][:, freq_idxs])
        plt.show()
        fourier_coefficients = np.linalg.solve(self.F_inv[freq_idxs][:, freq_idxs], s[freq_idxs])
        s_reconstructed = np.dot(self.F_inv[:, freq_idxs], fourier_coefficients)
        return s_reconstructed

    def sample_data_driven(self, s, freq_idxs):
        sample_pos = np.where(s != 0)[0]
        print('using %d samples and %d frequencies'%(len(sample_pos), len(freq_idxs)))
        fourier_coefficients, _, _, _ = np.linalg.lstsq(self.F_inv[sample_pos][:, freq_idxs], s[sample_pos], rcond=None)
        s_reconstructed = np.dot(self.F_inv[:, freq_idxs], fourier_coefficients)
        return s_reconstructed

    def reconstruction_plot(self, s, s_hat, n_steps=20):
        indices = np.argsort(-np.abs(s_hat))
        n_freqs = len(indices)//n_steps
        errs = []
        steps = []
        for step in range(n_steps):
            freqs = np.sort(indices[: np.maximum(step * n_freqs, 1)])
            recon = self.sample(s, freqs)
            err = np.linalg.norm(s - recon)
            errs.append(err)
            steps.append(len(freqs))

        return errs, steps


class LatticeSignalProcessing:
    def __init__(self, P):
        """
        :param P: partial order matrix of size |L| x |L|, P_ij = 1 iff j <= i
        """
        self.F_inv = P
        self.F = np.linalg.inv(P)
        self.FR = P.T
        self.n = self.F.shape[0]

    def ft(self, s):
        return np.dot(self.F, s)

    def ift(self, s):
        return np.dot(self.F_inv, s)

    def fr(self, h):
        return np.dot(self.FR, h)

    def conv(self, h, s):
        return self.ift(self.fr(h) * self.ft(s))

    def sample(self, s, freq_idxs):
        """F_inv is lower triangular, therefore to sample frequencies, we just take frequencies,
           maybe make argument s a function in the future"""
        print('solving ...')
        plt.imshow(self.F_inv[freq_idxs][:, freq_idxs])
        plt.show()
        fourier_coefficients = np.linalg.solve(self.F_inv[freq_idxs][:, freq_idxs], s[freq_idxs])
        s_reconstructed = np.dot(self.F_inv[:, freq_idxs], fourier_coefficients)
        return s_reconstructed

