function [A] = Matrix_A(n)
% test matrix as provided by Jos Maubach

A = kron(eye(n), ...
         gallery('tridiag',n,-1,2,-1)) + ...
         kron(gallery('tridiag',n,-1,2,-1), eye(n));

end

