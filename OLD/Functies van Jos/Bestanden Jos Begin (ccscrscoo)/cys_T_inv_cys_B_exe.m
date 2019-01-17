function cys_T_inv_cys_B_exe()

% examples
clear, clc
U = [1, 2, 0; 0, 1, 1; 0,0,-1]
B = speye(3);
cys_T_inv_cys_B_example(U, B, '1');

U = triu(magic(6))
B = indexvalue(6,10)
cys_T_inv_cys_B_example(U, B, '2');

U = triu(ones(6,6))
B = speye(6); full_B = full(B)
cys_T_inv_cys_B_example(U, B, '3');

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
U = triu(ones(6,6))
B = speye(6); full_B = full(B)
cys_T_inv_cys_B_example(U, B, '4');

% reverse roles of U and B
U = tril(magic(6))'
B = indexvalue(10,6)
cys_T_inv_cys_B_example(U, B, '5');

function cys_T_inv_cys_B_example(U, B, nb)
epsilon = 2^-42;
% U MUST BE upper triangular
if norm(U - triu(U), 'fro') > epsilon
 fprintf('Err(cxs_T_inv_cxs_B_example): Invalid U, is not upper triangular.\n');
 abort
end
if size(U,2) == size(B,1)
 fprintf('Multiply inv(U)*B where U = \n'); disp(full(U))
 fprintf('Multiply inv(U)*B where B = \n'); disp(full(B))
 [l, k, crsUT, iUT, vUT] = mat2ccs(U');
 [m, n, crsBT, iBT, vBT] = mat2ccs(B');
 [m, n, crsUinvBT, iUinvBT, vUinvBT] = cys_T_inv_cys_B(l, k, crsUT, iUT, vUT, m, n, crsBT, iBT, vBT);
 U_inv_B = ccs2mat(m, n, crsUinvBT, iUinvBT, vUinvBT)'; full_U_inv_B_css = full(U_inv_B)
 full_U_inv_B_exact = U\B; full(full_U_inv_B_exact)
 fprintf('Amount of stored entries in U_inv_B: %d\n', length(iUinvBT));
 if norm(U_inv_B - full_U_inv_B_exact) > epsilon, fprintf('Msg(cys_U_inv_cys_B_exe): %sA. FAILED\n',nb); else fprintf('Msg(cys_U_inv_cys_B_exe): %sA. PASSED\n',nb); end;
end
if size(U,1) == size(B,2)
 % now we need lower triangular L, use L = U' 
 L = U';
 fprintf('Multiply B*inv(Utp) where Utp = \n'); disp(full(L))
 fprintf('Multiply B*inv(Utp) where B   = \n'); disp(full(B))
 [k, l, crsU, iU, vU] = mat2ccs(L);
 [n, m, crsB, iB, vB] = mat2ccs(B);
 [n, m, crsBUinv, iBUinv, vBUinv] = cys_T_inv_cys_B(k, l, crsU, iU, vU, n, m, crsB, iB, vB);
 B_U_inv = ccs2mat(n, m, crsBUinv, iBUinv, vBUinv); full_B_U_inv_css = full(B_U_inv)
 full_B_U_inv_exact = B*inv(L); full(full_B_U_inv_exact)
 if norm(B_U_inv - full_B_U_inv_exact) > epsilon, fprintf('Msg(cys_U_inv_cys_B_exe): %sB. FAILED\n',nb); else fprintf('Msg(cys_U_inv_cys_B_exe): %sB. PASSED\n',nb); end;
end
