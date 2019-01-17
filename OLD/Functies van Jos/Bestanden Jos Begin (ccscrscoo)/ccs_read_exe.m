function ccs_read_exe()

% Read ccs data
% examples
% square n == m and sparse
A = stokes17();
[n, m, ccs, i, v] = mat2ccs(A);
ccs_write(n, m, ccs, i, v, 'ccs.data');
[n2, m2, ccs2, i2, v2] = ccs_read('ccs.data');
B = ccs2mat(n2, m2, ccs2, i2, v2);
max_difference = max(max(abs(A - B))), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end
max_difference = norm(A - B, 'inf'), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end
%
% non-square n < m
A = magic(8); A = A(1:5,:)
[n, m, ccs, i, v] = mat2ccs(A);
ccs_write(n, m, ccs, i, v, 'ccs.data');
[n2, m2, ccs2, i2, v2] = ccs_read('ccs.data');
B = ccs2mat(n2, m2, ccs2, i2, v2);
max_difference = max(max(abs(A - B))), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end
max_difference = norm(A - B, 'inf'), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end
%
% non-square n > m
A = magic(8); A = A(:,1:5)
[n, m, ccs, i, v] = mat2ccs(A);
ccs_write(n, m, ccs, i, v, 'ccs.data');
[n2, m2, ccs2, i2, v2] = ccs_read('ccs.data');
B = ccs2mat(n2, m2, ccs2, i2, v2);
max_difference = max(max(abs(A - B))), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end
max_difference = norm(A - B, 'inf'), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end


% Read ccsc data
% examples
A = stokes17();
[n, m, ccs, i, v] = mat2ccs(A);
ccs_write(n, m, ccs, i, v, 'ccsc.data', 1);
[n2, m2, ccs2, i2, v2] = ccs_read('ccsc.data', 1);
B = ccs2mat(n2, m2, ccs2, i2, v2);
max_difference = max(max(abs(A - B))), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end
max_difference = norm(A - B, 'inf'), if max_difference == 0, fprintf('Msg(ccs_read_exe): PASSED\n'); else fprintf('Msg(ccs_read_exe): FAILED\n'); end
