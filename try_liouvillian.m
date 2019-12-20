function L = try_liouvillian(lambda, g, N, C)

times = 0;
for i = 1:C+1,  for j = 1:C+1
	k = 0:C;  l = 0:C;
	[k,l] = meshgrid(k,l);
	s = super(i-1, j-1 ,k, l, lambda, g, N, C);
	times = times+1;
	s = sparse(s);
	L(times,:) = s;
end, end

end

function s = super(n,m,k,l,lambda, g, N, C)
    M = sqrt(N*(N+1));
    s1 = lambda/2*(n.*(n - 1)).^.5.*(k == n-2).*(l == m);

    s2 = lambda/2*(m.*(m - 1)).^.5.*(k == n).*(l == m-2);

    s3 = -lambda/2*((n+1).*(n+2)).^.5.*(k == n+2).*(l == m);

    s4 = -lambda/2*((m+1).*(m+2)).^.5.*(k == n).*(l == m+2);

    s5 = -(n+m).*(k == n).*(l == m);

    s6 = 2*((n+1).*(m+1)).^.5.*(k == n+1).*(l == m+1);

    s7 = -.5*g^2*(n.*(n-1)+m.*(m-1)).*(k == n).*(l == m);

    s8 = g^2*((n+1).*(n+2).*(m+1).*(m+2)).^.5.*(k == n+2).*(l == m+2);



    %% extra terms due to squeezed reservoir

    s9 = 2*(n.*m).^.5.*(k == n-1).*(l == m-1) - ...
(n+1).*(k == n).*(l == m) - (m+1).*(k == n).*(l == m);

    s10 = 2*(n+1).^.5.*m.^.5.*(k == n+1).*(l == m-1) ...
	- ((n+1).*(n+2)).^.5.*(k == n+2).*(l == m) ...
	- (m.*(m-1)).^.5.*(k == n).*(l == m-2);

    s11 = 2*(n.*(m+1)).^.5.*(k == n-1).*(l == m+1) ...
	- (n.*(n-1)).^.5.*(k == n-2).*(l == m) ...
	-((m+1).*(m+2)).^.5.*(k == n).*(l == m+2);



    s = s1+s2+s3+s4+(N+1)*(s5+s6)+s7+s8+ N*s9 - M*s10 - conj(M)*s11;

    s = reshape(s, (C+1)^2, 1);

    s = s';

end
