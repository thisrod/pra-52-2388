function L = first_principles_liouvillian(lda, g, N, C)
	M = sqrt(N*(N+1));
	dC = 2;	% extra Fock states to avoid truncation errors
	a = spdiags(sqrt(0:C+dC)', 1, C+1+dC, C+1+dC);  % annihilation operator
	[aR, Ra] = raise(a);
	[AR, RA] = raise(a');

	L = lda/2*(AR^2 - RA^2 -aR^2 + Ra^2) ...
		+ g^2/2*(2*aR^2*RA^2 - AR^2*aR^2 - Ra^2*RA^2) ...
		+ (N+1)*(2*aR*RA - AR*aR - Ra*RA) ...
		+ N*(2*AR*Ra - aR*AR - RA*Ra) ...
		- M*(2*aR*Ra - aR^2 -Ra^2) ...
		-M*(2*AR*RA - AR^2 - RA^2);

	% project out the extra Fock states
	% X is a density matrix where the elements to keep are true, and the
	% ones to discard are false.
	X = false(C+1+dC);
	X(1:C+1, 1:C+1) = true;
	peye = speye((C+1+dC)^2);
	peye = peye(X(:),:);
	L = peye*L*peye';
end