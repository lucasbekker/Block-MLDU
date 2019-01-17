function [n, m, css, i, lu] = coo2css(n, m, i, j, v)
s=(i>=j);
[n, m, ccs, ic, vc] = coo2ccs(n, m, i(s), j(s), v(s));
s=(i<=j);
% The next THIS OPERATION LEADS TO coo-UNSORTED INPUT, PERHAPS SORT HERE IF 'SORTED' SHOULD BE coo-INVARIANT
[m, n, crs, jc, wc] = coo2ccs(m, n, j(s), i(s), v(s));
css = [];
i = [];
lu = zeros(0,2);
mx = max(n,m);
c = min(length(ccs),length(crs))-1;
for k=1:c
 %[mx, ind, l, u] = sparse_stretch(mx, ic(ccs(k):ccs(k+1)-1), vc(ccs(k):ccs(k+1)-1), mx, jc(crs(k):crs(k+1)-1), wc(crs(k):crs(k+1)-1));
 [mx, ind, dlu] = sparse_join(mx, ic(ccs(k):ccs(k+1)-1), vc(ccs(k):ccs(k+1)-1), mx, jc(crs(k):crs(k+1)-1), wc(crs(k):crs(k+1)-1));
 css = [css; length(ind)];
 i = [i; ind];
 %lu = [lu; [l, u]];
 lu = [lu; dlu];
end
% append the cols/rows only stored in matrix 1
css = [css; diff(ccs(c+1:end))];
i  = [i;   ic(ccs(c+1):end)]; % if i = [] -> i = 1 x 0 (octave 3.8.2)
lu = [lu; [vc(ccs(c+1):end), zeros(ccs(end)-ccs(c+1),1)]];
% append the cols/rows only stored in matrix 2
css = [css; diff(crs(c+1:end))];
i  = [i; jc(crs(c+1):end)]; % if i = 1 x 0 -> i = 2 x 0 (octave 3.8.2)
lu = [lu; [zeros(crs(end)-crs(c+1),1), wc(crs(c+1):end)]];
css = cumsum([1;css]);

% if length(i) == 0, i should be 0 x 0
if length(i) == 0; i = []; end
