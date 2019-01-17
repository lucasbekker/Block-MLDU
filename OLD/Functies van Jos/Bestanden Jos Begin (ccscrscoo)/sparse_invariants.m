function check = sparse_invariants(n, i, v)
epsilon = 2^-36;

% invariant: max(i) <= n: must also be defined also for i = []
% invariant: size(i, 1) == size(v, 1) (n(i) == n(v), i.e., same amount of rows, v can have multiple columns)
check = (n >= 0);
if (n < 0), fprintf('Err(sparse_invariants): %d = n < 0.\n', n); end;
if length(i) > 0

 %mx = max(i);
 %check = (mx <= n);
 %if (mx > n), fprintf('Err(sparse_invariants): %d = max(i) > n = %d.\n', mx, n); end;

 % invariant: entries v_k with i_k>n should be zero
 inconsistentnb = sum(i>n);
 inconsistentvalue = sum(abs(v(i>n)));
 check = (inconsistentvalue < epsilon);
 if (inconsistentnb > 0 && inconsistentvalue < epsilon)
  fprintf('Wrn(sparse_invariants): %d indices > %d related to zero values (is ok with invariants).\n', inconsistentnb, n);
 end;
 if (inconsistentnb > 0 && inconsistentvalue >= epsilon)
  fprintf('Err(sparse_invariants): %d indices > %d related to non-zero values v, of which sum abs(v(i>n)) is: %g\n', inconsistentnb, n, inconsistentvalue);
  fprintf('Violation is at indices:\n'); disp(i');
 end;

 df = abs(size(i, 1) - size(v, 1));
 check = check && (df == 0);
 if (df > 0), fprintf('Err(sparse_invariants): %d = size(i, 1) <> size(v, 1) = %d.\n', size(i, 1), size(v, 1)); end;
 %ascent = (sum(diff(i) > 0) == size(i,1) - 1)
 descent = (diff(i) <= 0);
 if sum(descent) > 0, k = find(descent); k = k(end); end;
 check = check && (sum(descent) == 0);
 if (sum(descent) > 0), fprintf('Err(sparse_invariants): %d = i_%d >= i_%d = %d.\n', i(k), k, k+1, i(k+1)); end;
end;
%if ~check, abort; end;
