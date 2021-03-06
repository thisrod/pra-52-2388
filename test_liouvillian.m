% compare Liouvillian terms in paper to those constructed from raise

C = 3;  dC = 1;
[i,j,n,m] = ndgrid(0:C);

a = spdiags(sqrt(0:C+dC)', 1, C+1+dC, C+1+dC);  % annihilation operator
[aR, Ra] = raise(a);
[AR, RA] = raise(a');

terms = {
	AR^2 i i-1 2 0; ...
	Ra^2 j j-1 0 2; ...
	aR^2 i+1 i+2 -2 0; ...
	RA^2 j+1 j+2 0 -2; ...
	aR*RA i+1 j+1 -1 -1; ...
	AR*Ra i j 1 1; ...
% 	i+1 j -1 1; ...
% 	j+1 i 1 -1; ...
% 	i+1 i+2 -2 0; ...
% 	j+1 j+2 0 -2; ...
% 	j j-1 0 2; ...
% 	i i-1 2 0; ...
};

Ls = {};
for t = terms.'
	Ls = {Ls{:} sqrt(t{2}.*t{3}).*(i==n+t{4} & j==m+t{5})};
end

odd_terms = { ...
	AR^2*aR^2 i.*(i-1) .*(i==n & j==m); ...
	Ra^2*RA^2 j.*(j-1).*(i==n & j==m); ...
	aR^2*RA^2 sqrt((i+1).*(i+2).*(j+1).*(j+2)).*(i==n-2 & j==m-2); ...
	AR*aR i.*(i==n & j==m); ...
	Ra*RA j.*(i==n & j==m); ...
	aR*AR (i+1).*(i==n & j==m); ...
	RA*Ra (j+1).*(i==n & j==m); ...
};

Ls = {Ls{:} odd_terms{:,2}};
Ls = cellfun(@(L) sparse(reshape(L, (C+1)^2, (C+1)^2)), Ls, 'UniformOutput', false);

Ms = {terms{:,1} odd_terms{:,1}};

X = false(C+1+dC);
X(1:C+1, 1:C+1) = true;
peye = speye((C+1+dC)^2);
peye = peye(X(:),:);
Ms = cellfun(@(M) reshape(peye*M*peye', (C+1)^2, (C+1)^2), Ms, 'UniformOutput', false);

for i = 1:length(Ms)
	d = norm(Ls{i}-Ms{i}, 'fro');
	if d < 1e-3, continue, end
	fprintf('Term %d: |L-M| = %.3f\n', i, d)
end
