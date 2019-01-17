function cxs_T_inv_cxs_B_exe()

% examples
clear, clc
L = [1, 2, 0; 0, 1, 1; 0,0,-1]'
B = speye(3);
cxs_T_inv_cxs_B_example(L, B, '1');

L = tril(magic(6))
B = indexvalue(6,10)
cxs_T_inv_cxs_B_example(L, B, '2');

L = tril(magic(6))
B = indexvalue(6,10)'
cxs_T_inv_cxs_B_example(L, B, '3');

% % A = diag(ones(6,1),0) - diag(ones(5,1),1)
%
%   1  -1  -0  -0  -0  -0
%  -0   1  -1  -0  -0  -0
%  -0  -0   1  -1  -0  -0
%  -0  -0  -0   1  -1  -0
%  -0  -0  -0  -0   1  -1
%  -0  -0  -0  -0  -0   1
%
% % is the inverse of
%
%   1   1   1   1   1   1
%  -0   1   1   1   1   1
%  -0  -0   1   1   1   1
%  -0  -0  -0   1   1   1
%  -0  -0  -0  -0   1   1
%  -0  -0  -0  -0  -0   1
L = tril(ones(6,6))
B = speye(6); full_B = full(B)
cxs_T_inv_cxs_B_example(L, B, '4');

function cxs_T_inv_cxs_B_example(L, B, nb)
epsilon = 2^-42;
% L MUST BE lower triangular
if norm(L - tril(L), 'fro') > epsilon
 fprintf('Err(cxs_T_inv_cxs_B_example): Invalid L, is not lower triangular.\n');
 abort
end
if size(L,2) == size(B,1)
 fprintf('Multiply inv(L)*B where L MUST BE LOWER TRIANGULAR L = \n'); disp(full(L))
 fprintf('Multiply inv(L)*B where B = \n'); disp(full(B))
 [k, l, crsL, iL, vL] = mat2ccs(L);
 [m, n, crsBT, iBT, vBT] = mat2ccs(B');
 [m, n, crsLinvBT, iLinvBT, vLinvBT] = cxs_T_inv_cxs_B(k, l, crsL, iL, vL, m, n, crsBT, iBT, vBT, [], 1);
 L_inv_B = ccs2mat(m, n, crsLinvBT, iLinvBT, vLinvBT)'; full_L_inv_B_css = full(L_inv_B)
 full_L_inv_B_exact = L\B; full(full_L_inv_B_exact)
 fprintf('Amount of stored entries in L_inv_B: %d\n', length(iLinvBT));
 if norm(L_inv_B - full_L_inv_B_exact) > epsilon, fprintf('Msg(cxs_T_inv_cxs_B_exe): %sA. FAILED\n',nb); else fprintf('Msg(cxs_T_inv_cxs_B_exe): %sA. PASSED\n',nb); end;
end
if size(L,1) == size(B,2)
 % now we need upper triangular U, use U = L' 
 U = L';
 fprintf('Multiply B*inv(U) where U MUST BE UPPER TRIANGULAR U = \n'); disp(full(U))
 fprintf('Multiply B*inv(U) where B   = \n'); disp(full(B))
 [l, k, ccsUT, iUT, vUT] = mat2ccs(U');
 [n, m, ccsB, iB, vB]    = mat2ccs(B);
 [n, m, crsBUinv, iBUinv, vBUinv] = cxs_T_inv_cxs_B(l, k, ccsUT, iUT, vUT, n, m, ccsB, iB, vB, [], 2);
 B_U_inv = ccs2mat(n, m, crsBUinv, iBUinv, vBUinv); full_B_U_inv_css = full(B_U_inv)
 full_B_U_inv_exact = B*inv(U); full(full_B_U_inv_exact)
 if norm(B_U_inv - full_B_U_inv_exact) > epsilon, fprintf('Msg(cxs_T_inv_cxs_B_exe): %sB. FAILED\n',nb); else fprintf('Msg(cxs_T_inv_cxs_B_exe): %sB. PASSED\n',nb); end;
end
