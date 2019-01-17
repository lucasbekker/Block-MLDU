function [n, i, v] = sparse_stack(n, i, v, is, vs)
% does not check whether i(end) < is, could call sparse_invariants for this
i = [i; is];
v = [v; vs];

