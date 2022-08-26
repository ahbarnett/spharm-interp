function A = spharmmat(z, phi, P)
% SPHARMMAT  matrix of spherical harmonics degree<=P at arbitrary pts on S2
%
%  A = spharmmat(z, phi, P) where z is vector of N z-coords and phi an
%   length-N vector of phi coords, defining N targets (z(j),phi(j)) for j=1...N,
%   on S2, returns A the matrix of size N by (P+1)^2
%   which evaluates all spherical harmonics up to degree P at all the targets.
%
% Notes: 1) normalization matches h3dprojloc.
% 2) col ordering should match coeff output ordering of flattencnm.m
% 3) all MATLAB for now. But probably memory-access not flop bound.
% 4) this is faster than SPHARMLATMAT when need points not on same z,
%   but probably slower if matrix needed on a grid (has repeated legendre calls)

% todo: interface to Fortran

% Barnett 8/3/15
if nargin==0, test_spharmmat; return; end
N = numel(phi); if numel(z)~=N, error('phi and z must have same length!'); end
z = z(:); phi = phi(:);     % make column vectors
A = zeros(N,(P+1)^2);       % the output
allords = (-P:P)'; phases = exp(1i*phi*allords');     % precompute, N by 2P+1
ind = 0;    % index offset as step through cols of A
for p=0:P
  ym = sqrt(2)*legendre(p,z,'norm'); % match MATLAB to h3dprojloc normalization
  ym = [ym(end:-1:2,:); ym];   % stack negative orders; ym is 2p+1 by M
  ords = (-p:p)';              % col vec of orders for this degree
  ym = ym .* repmat((-1).^ords,[1 N]);   % match phase conventions
  A(:,ind+(1:2*p+1)) = ym.' .* phases(:,P+1+ords);  % all orders at once
  ind = ind+2*p+1;
end

%%%%
function test_spharmmat      % just checks against spharmgrideval
M = 30; N = 60;  % # nodes in z and phi, can be rect
[z w] = gauss(M);
phi = 2*pi*(1:N)/N; % first var is z (fast), second is phi (slow)
P = 20;         % max degree, somewhat less than 0.5 M
[pphi zz] = meshgrid(phi,z);  % makes z change fast, phi slow
tic; A = spharmmat(zz, pphi, P); fprintf('spharmmat in %.3g s\n',toc)   % do it
c = zeros(P+1,2*P+1);
%deg = 1; ord = 1; c(deg+1,ord+P+1) = 1; % check a single (n,m), or...
for n=0:P, for m=-n:n, c(n+1,m+P+1) = randn+1i*randn; end, end
cnm = flattencnm(c);
u = A*cnm;       % evaluate sph harm expansion via matvec
u = reshape(u, [M N]);
figure; showsphgrid(u,z);
ue = spharmgrideval(cnm,z,phi);   % grid eval
fprintf('grid l2 error is %.3g\n', norm(u(:)-ue(:)))



