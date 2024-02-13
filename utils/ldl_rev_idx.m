function [U,D] = ldl_rev_idx(A)
%ldl_rev_idx Block ldl' factorization for Hermitian indefinite matrices
%with reversed indexing
%   [U,D] = ldl_rev_idx(A) stores a block diagonal matrix D and a permuted
%   unit upper triangular matrix (i.e. a product of unit upper triangular
%   and permutation matrices) in L so that A = U*D*U'.

reverse_idx = @(X) X(end:-1:1, end:-1:1);
[U,D] = ldl(reverse_idx(A));
U = reverse_idx(U);
D = reverse_idx(D);
end

