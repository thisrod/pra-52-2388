function [L, R] = raise(A)
%RAISE raise a matrix to a superoperator
%    RAISE(A), where A is an N*N square matrix, returns N^2*N^2 sparse
%    matrices L and R such that the following identities are satisfied
%    by every N*N matrix X:
%
%        A*X = reshape(L*X(:), N, N)  and  X*A = reshape(R*X(:), N, N).
%
%    Recall that X(:) is expanded in column-major order.

assert(ndims(A) == 2)
N = size(A,1);
assert(N == size(A,2))

L = kron(speye(N), A);
R = kron(A.', speye(N));
end