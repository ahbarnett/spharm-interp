% demo example code for interpolation on the sphere
% Barnett 8/27/22
clear

kvec = [9;2;6];            % wavevector (col) controlling f oscillation
f = @(x) exp(1i*x*kvec);   % underlying complex func on R^3 (a plane wave)
M = 1e5;                   % how many data pts on S^2
sigma = 0.1;               % additive noise level
z = rand(M,1)*2-1; phi = 2*pi*rand(M,1);   % iid rand pts on S^2
rho = sqrt(1-z.^2);
x = [rho.*cos(phi), rho.*sin(phi), z];    % M*3 coords in R^3
ynoiseless = f(x);
y = ynoiseless + sigma*randn(M,1);
figure; subplot(1,2,1); ptarea = 1.0;
scatter3(x(:,1),x(:,2),x(:,3),ptarea*ones(M,1),real(y),'filled');
colormap(jet(256)); colorbar; axis vis3d; grid off
title('data at M=10^5 sphere points (Re part)')
P = 20;                    % spherical harmonic expansion max degree
fprintf('fitting at M=%d pts, max degree P=%d...\n',M,P), tic
[coeffs resid] = lsqsolvespharm(y, z, phi, P);     % fit the expansion
fprintf('done in %.3g s\n',toc)
subplot(1,2,2); showspharmexp(coeffs);
title('fitted expansion of max degree P=20 (Re part)')
fprintf('Eval slow exact expansion at M=%d pts...\n',M), tic
fSH = spharmeval(coeffs, z, phi);
rmsr = norm(fSH-y)/sqrt(M);
fprintf('done in %.3g s. RMS residual = %.3g (cf sigma=%.3g)\n',toc,rmsr,sigma)
rmse = norm(fSH-ynoiseless)/sqrt(M);
fprintf('RMS error vs noiseless func = %.3g (cf sigma*P/sqrt(M)=%.3g)\n',rmse,sigma*P/sqrt(M))

% generate the plot used in the README:
%print -dpng demo.png
