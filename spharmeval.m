function u = spharmeval(cvec, z, phi)
% SPHARMEVAL  evaluate spherical harmonic expansion at given targets on S2
%
% u = spharmeval(cvec, z, phi) interprets cvec as an unrolled vector of (P+1)^2
% complex spherical harmonic coefficients, and evaluates the resulting
% expansion at the
% arbitrary unit sphere points (z,phi), where z and phi are vectors or arrays
% of the same size (with number of targets elements).
% Each entry of z=cos(theta) is in [-1,1], and phi is azimuthal angle.
% The returned u is the same size as z and phi.
%
% Native MATLAB code, vectorized over targets -- not a performance code!
% (Uses O(PM) storage where M = # targets)
%
% Also see (for evaluation on grids): SPHARMGRIDEVAL

% Barnett 8/25/22, to complete the basic tools.

if nargin==0, test_spharmeval; return; end

P = sqrt(numel(cvec))-1;
if abs(P-round(P))>1e-14, error('cvec is not length (P+1)^2 for integer P'); end
P = round(P);
siz = size(z); M = numel(z);
if size(z)~=size(phi), error('z and phi must have the same size!'); end
u = zeros(1,M);
ind = 0;
for n=0:P       % degree
  Pnonnegm = legendre(n,z(:),'norm');             % (n+1) by M,  orders m=0...n
  eimp = exp(1i*((0:n)'*phi(:)'));         % phases via outer prod
  c = sqrt(2)*cvec(ind+(1:2*n+1));            % coeffs for n, note normalization
  c = c(:).' .* (-1).^(-n:n);              % row vec, fix phases (-1)^m, why?
  u = u + c(n+1:end) * (Pnonnegm.*eimp);
  u = u + c(n:-1:1) * (Pnonnegm(2:end,:).*conj(eimp(2:end,:)));  % m<0 parts
  ind = ind + 2*n+1;
end
u = reshape(u,siz);

%%%%%%%%%%%%%%
function test_spharmeval
Nz = 100; Np = 200; % # nodes in each direction
[z w] = gauss(Nz);
phi = 2*pi*(0:Np-1)/Np; % first var is z (fast), second is phi (slow)
P = 30;         % max degree
PP = (P+1)^2;
cvec = randn(PP,1)+0*1i*randn(PP,1);   % pick some coeffs
                                       % P = 1; cvec = [0 1i 0 0].';  % debug
tic; ug = spharmgrideval(cvec,z,phi); fprintf('spharmgrideval in %.3g s\n',toc)
figure; showsphgrid(ug,z);
[zz pphi] = ndgrid(z, phi);
tic; u = spharmeval(cvec,zz,pphi); fprintf('spharmeval in %.3g s\n',toc)
figure; showsphgrid(u,z);
fprintf('rel l2 err from grid eval: %.3g\n', norm(ug(:)-u(:))/norm(u(:)))
%u(1), ug(1), u(1)/ug(1)  % debug

