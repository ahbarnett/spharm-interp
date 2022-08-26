function u = spharmgrideval(cnm,z,phi)
% SPHARMGRIDEVAL  evaluate spherical harmonic expansion on tensor grid on S2.
%
% u = spharmgrideval(cnm,z,phi) evaluates u a M*N complex array given size-M
%  grid of z, and size-N grid of phi, defining a tensor-product grid on the
%  sphere, and the length-(P+1)^2 list of coefficients c_nm.
% 
%     u =    sum      sum  c_nm Y_nm(theta,phi)   where z = cos(theta)
%           n=0..P  |m|<=n
%
%  The fast (row) index of u is z, the slow (col) is phi, and normalization
%  matches projloc3d.
%
% Notes: 1) We allow M to differ from N, unlike h3dprojloc
% 2) now O(K^3) where P and N and M are of order K.
% 3) self-test also tests spharmgridevalf
%
% See also: SPHARMGRIDEVALF (Fortran/MEX version)

% Barnett 8/2/15. 10x faster via less legendre calls & exps, O(K^4) 8/3/15
% another 10x faster, O(K^3), where K~P~N~M, and avoiding fft. 8/4/15
if nargin==0, test_spharmgrideval; return; end
N = numel(phi);
M = numel(z);
if ~isvector(cnm), error('cnm must be in flat vector format not matrix'); end
P = sqrt(numel(cnm))-1;
if abs(P-round(P))>1e-14, error('cnm is not length (P+1)^2'); end
P = round(P);
fouc = zeros(M,2*P+1);   % Fourier coeffs on each latitude ring
ind = 1;           % index in cnm
for p=0:P
  ym = sqrt(2)*legendre(p,z,'norm').';    % match to h3dprojloc normalization
  ym = [ym(:,end:-1:2), ym];   % stack negative orders; ym now M by 2p+1
  ords = (-p:p);          % col vec of orders for this degree
  fouc(:,P+1+ords) = fouc(:,P+1+ords) + ym .* repmat((-1).^ords .* cnm(ind+p+ords).',[M 1]);   % phase & coeffs, add to Fou coeffs of these orders
  ind = ind+2*p+1;
end
% eval Fourier series on each ring (can be underresolved unlike for FFT):
fser = exp(1i*(-P:P)'*phi(:)');  % the DFT-like matrix eval the Fourier series
u = fouc * fser;    % matmat, O(K^3) too but v fast, handles general phi grid

%%%%%
function test_spharmgrideval    % evals then projects back, needs spharmproj
Nz = 200; Np = 400; % # nodes in each direction
[z w] = gauss(Nz);
phi = 2*pi*(0:Np-1)/Np; % first var is z (fast), second is phi (slow)
P = 180;         % max degree, somewhat less than 0.5 N, for projection test
PP = (P+1)^2;
cnm = randn(PP,1)+1i*randn(PP,1);   % pick some coeffs
tic; u = spharmgrideval(cnm,z,phi); fprintf('matlab grideval in %.3g s\n',toc)
showsphgrid(u,z);
tic; cp = spharmproj(u, z, w, P); fprintf('proj in %.3g s\n',toc)
cnmp = flattencnm(cp);
disp('grid-eval-then-proj-back error:');
norm(cnm-cnmp)   % projecting back gives same coeffs?

% also test fortran version:
tic; u = spharmgridevalf(cnm,z,phi); fprintf('MEX grideval in %.3g s\n',toc)
tic; cp = spharmproj(u, z, w, P); fprintf('proj in %.3g s\n',toc)
cnmp = flattencnm(cp);
disp('grid-eval-then-proj-back error:');
norm(cnm-cnmp)   % projecting back gives same coeffs?
