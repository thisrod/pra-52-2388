% Reproduce the results by Munro & Reid

C = 200;	% Fock number cutoff
% N = 1;  g = 2.5;  lambda = 100*g^2;  T = 0.01;
N = 100;  g = 2.5;  lambda = 100*g^2;  T = 0.01;
% N = 20;  g = 0.1;  lambda = 100*g^2;  T = 1;

tt = T*(0:14)/14;

L = first_principles_liouvillian(lambda, g, N, C);

R0 = zeros(C+1);  R0(1,1) = 1;  % initial vacuum state
R0 = R0(:);

tic
[tt, R] = ode45(@(~,R) L*R, tt, R0);
fprintf('Integration time: %.1f\n', toc)

[ns,ms] = ndgrid(0:C);
ns = ns(:)';  ms = ms(:)';

figure
plot(tt, sum(R(:, ns==ms), 2)-1, '.k')
xlabel t,  ylabel 'Tr \rho - 1',  title 'Norm check'

K = 150;
dz = 40/K;  zz = (0:K)'*dz - 20;
dp = 6/K;  pp = (0:K)'*dp - 3;

% set up observables
Pz = ones(length(zz), (C+1)^2);  Pp = ones(length(pp), (C+1)^2);
tic
for n = 0:C
	zn = focksum([zeros(1,n) 1] ,zz);
	pn = (-i)^n*focksum([zeros(1,n) 1] ,pp);
	
	Pz(:,ns==n) = Pz(:,ns==n).*zn;
	Pz(:,ms==n) = Pz(:,ms==n).*conj(zn);
	Pp(:,ns==n) = Pp(:,ns==n).*pn;
	Pp(:,ms==n) = Pp(:,ms==n).*conj(pn);
end
fprintf('Hermite time: %.1g\n', toc)

warning 'How big are the imaginary parts of P(z) and P(p)?'

figure
subplot 211
waterfall(zz,tt,real(Pz*R.')')
A = gca;  A.Colormap = [0 0 0];
view(-134.3,15.6),  pbaspect([1 1 0.4])
xlabel z, ylabel t, zlabel P(z)

subplot 212
waterfall(pp,tt,real(Pp*R.')')
A = gca;  A.Colormap = [0 0 0];
view(-134.3,15.6),  pbaspect([1 1 0.4])
xlabel p, ylabel t, zlabel P(p)

figure
subplot 121
p = real(Pz*R(end,:)');
tp = dz*sum(p);
plot(zz,p)
title(sprintf('Total probability = 1 %+.0e', tp-1))
xlabel z, ylabel P(z)

subplot 122
p = real(Pp*R(end,:)');
tp = dp*sum(p);
plot(pp,p)
title(sprintf('Total probability = 1 %+.0e', tp-1))
xlabel p, ylabel P(p)



