The standalone python program main_generateChimera_generalized.py generates
the edge list of nonplanar, effectively 2D lattice graphs, based on complete
bichromatic subgraphs, in DIMACS format.


The standalone python program main_sitePercolation.py (main_bondPercolation.py)
contains function defintions for a union find based algorithm for site (bond)
percolation.  Albeit the implementation is slightly different, it is based on
the idea of the fast percolation algorithm detailed in

"Fast Monte Carlo algorithm for site and bond percolation"
Newman, M.E.J. and Ziff, R.M.
PRE 64 (2001) 016706

Exemplary usage:

python main_sitePercolation.py chimera_N32768.dimacs.gz 1

Therein, the bonfile needs to follow the dimacs edge-list format.
