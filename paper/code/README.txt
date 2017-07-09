This package implements in Matlab the following indices: 'prel', 'prec',
'rnia', 'ce', 'l&w', 'stm', 'wjac', 'wdic', 'fabia', 'u', 'e', 'ay',
'erel', 'erec', 'csi', and 'ebc'. Their references and definitions can
be found in the paper "Similarity Measures for Comparing Biclusterings".


Usage example:

ba(1).rows = [1 3];
ba(1).cols = [1 2 5];
bb(1).rows = [1 3];
bb(1).cols = [1 2];
bb(2).rows = [4 5];
bb(2).cols = [1 4 5];
nr=5;
nc=6;
external_biclustering_indices('prel', ba, bb, nr, nc)

Author: Danilo Horta
