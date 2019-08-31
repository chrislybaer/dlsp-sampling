# dlsp-sampling
A minimal implementation accompaning the paper "Sampling Signals on Meet/Join Lattices" to reproduce the results.

@INPROCEEDINGS{Wend1911:Sampling,
AUTHOR={Chris Wendler and Markus {P{\"u}schel}},
TITLE="Sampling Signals on {Meet/Join} Lattices",
BOOKTITLE="2019 IEEE Global Conference on Signal and Information Processing
(GlobalSIP) (GlobalSIP 2019)",
ADDRESS="Ottawa, Canada",
DAYS=11,
MONTH=nov,
YEAR=2019,
KEYWORDS="lattice signal; lattice sampling; lattice Fourier transform; lattice shift;
meet; join; graph signal processing; algebraic signal processing",
ABSTRACT="We present a novel sampling theorem, and prototypical applications, for
Fourier-sparse lattice signals, i.e., data indexed by a finite semilattice.
A semilattice is a partially ordered set endowed with a meet (or join)
operation that returns the greatest lower bound (smallest upper bound) of
two elements. Semilattices can be viewed as a special class of directed
graphs with triangular adjacency matrix, which thus cannot be diagonalized.
Our work does not build on prior graph signal processing (GSP) frameworks
but on the recently introduced discrete-lattice signal processing (DLSP),
which uses the meet as shift operator to derive convolution and Fourier
transform. DLSP is fundamentally different from GSP in that it requires
several generating shifts which capture the partial-order- rather than the
adjacency-structure and a diagonalizing Fourier transform is always
guaranteed by algebraic lattice theory. We apply and demonstrate the
utility of our novel sampling scheme in three real world settings from
computational biology, document representation, and auction design."
}
