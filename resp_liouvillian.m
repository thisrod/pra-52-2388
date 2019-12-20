function L = resp_liouvillian(lda, g, N, C)
% superoperator from Equation 10.
% gamma = 1, this is a typo where the paper omits to scale the units

M = sqrt(N*(N+1));
[i,j,n,m] = ndgrid(0:C);

terms = {
	lda/2 i i-1 2 0; ...
	lda/2 j j-1 0 2; ...
	-lda/2 i+1 i+2 -2 0; ...
	-lda/2 j+1 j+2 0 -2; ...
	2*(N+1) i+1 j+1 -1 -1; ...
	2*N i j 1 1; ...
	-2*M i+1 j -1 1; ...
	-2*M j+1 i 1 -1; ...	% bingo
	M i+1 i+2 -2 0; ...
	M j+1 j+2 0 -2; ...
	M j j-1 0 2; ...
	M i i-1 2 0; ...
};

L = 0;
for t = terms.'
	L = L + t{1}*sqrt(t{2}.*t{3}).*(i==n+t{4} & j==m+t{5});
end

L = L - g^2/2*(i.*(i-1) + j.*(j-1)).*(i==n & j==m);
L = L + g^2*sqrt((i+1).*(i+2).*(j+1).*(j+2)).*(i==n-2 & j==m-2);
L = L - (N+1)*(i+j).*(i==n & j==m);
L = L - N*(i+j+2).*(i==n & j==m);

L = reshape(L, (C+1)^2, (C+1)^2);

L = sparse(L);

end
