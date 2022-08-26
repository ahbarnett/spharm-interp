function [cnm resvec] = lsqsolvespharm(b, z, phi, P, o)
% LSQSOLVESPHARM  least-squares solve sph harm coeffs given arbitrary pts on S^2.
%
% cnm = lsqsolvespharm(b, z, phi, P) solves the least-squares problem for a
%  truncated spherical harmonic expansion up to degree P, matching complex data
%  b given at arbitrary points (z,phi) on the sphere. Useful for cryo-EM to
%  get Fourier data on a k-shell given set of image rings at the same k.
%
%  Let S[c](z,phi) := \sum_{n=0}^P \sum_{m=-n}^n c_{nm} Y_{nm}(z) e^{im\phi}
%
%  define the a function on S^2 that is the spherical harmonic expansion with
%  length-(P+1)^2 coefficient vector c. Let (z,phi) define a point on S^2.
%
%  Then the routine returns
%
%        \hat{c}  =     arg min     \sum_{i=1}^m | b_i - S[c](z_i,\phi_i)|^2
%                   c in C^{(P+1)^2}
%
%  ie the least-squares solution to Ac=b where A is the matrix of spherical
%  harmonics evaluated at the data points.
%
%  It uses a fast algorithm that is O(m + P^3) where m is number of data points
%  (see meth='f' below). The m prefactor is roughly 10^-5, ie 1e5 points/sec
%  to fill the sparse interpolation matrix (single-core).
%
% Inputs:
%   b = complex length-m vector of values on points (the data values).
%   z,phi = length-m vectors of z and phi locations of data points on S^2.
%   P = maximum degree of spherical harmonic projections.
%
% Outputs:
%   cnm = complex coefficients returned as column vector, ordered (0,0), (1,-1),
%         (1,0), (1,1), (2,-2), ... (P,P).  Length is (P+1)^2.
%   resvec = residual history (if CG used)
%
% cnm = lsqsolvespharm(b, z, phi, P, opts) controls various options such as:
%     opts.verb = 0,1,.. amount of text or figure feedback (0=none, etc)
%     opts.meth = 'f' conjugate gradient (CG) on normal eqns, fast sparse
%                     interpolation matrix apply. (Default.) O(m + P^3)
%                     If not overdetermined, may fail to converge.
%                 'l' plain dense least-squares solve on matrix, O(m.P^4)
%                 'n' CG on normal equations, dense apply matrix, O(m.P^2)
%                     If not overdetermined, may fail to converge.
%                 'd' disc-patch-averaging & projection (obsolete)
%   To control accuracy in meth='f', you could vary:
%     opts.ppw  = point-per-wavelength for grid projection.
%     opts.q    = interpolation order for sparse matrix: eg 5,7,...
%   To control CG residual tolerance:
%     opts.tol  = tolerance in MATLAB's PCG call
%     opts.maxit = maximum # CG iterations (default 100)
%     opts.regparam = multiple of identity to regularize normal eqns by.
% Unspecified options are given sensible default values.
%
% See also: devel/explore_lsq.m for development & more tests.
%
% Barnett 8/12/15
% 0-indexed phi grid, better rect grid for proj, 9/1/15.
% regularization testing 8/30/16
if nargin==0, test_lsqsolvespharm; return; end
if nargin<5, o = []; end
if ~isfield(o,'meth'), o.meth = 'f'; end
if ~isfield(o,'ppw'), o.ppw = 10; end             % enough for 1e-5
if ~isfield(o,'q'), o.q = 7; end
if ~isfield(o,'tol'), o.tol = 1e-6; end
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'regparam'), o.regparam = 0; end
if ~isfield(o,'maxit'), o.maxit = 100; end          % for CG

if strcmp(o.meth,'l') | strcmp(o.meth,'n') | strcmp(o.meth,'c')
  tic;
  A = spharmmat(z,phi,P);                                          % O(m.P^2)
  if o.verb, fprintf('A filled %dx%d in %.3g s\n',size(A,1),size(A,2),toc), end
end

if strcmp(o.meth,'f')                      % the main method
 N = round(o.ppw*P);   % # phi nodes
 M = round(o.ppw*P/2); % # z nodes
 [zg wg] = fejer(M);
 phig = 2*pi*(0:N-1)/N;
 tic;
 L = fillsparseL(zg,phig,z,phi,o.q);
 if o.verb, fprintf('fill L in %.3g s\n',toc), end
 fastapplyA = @(c) L*reshape(spharmgrideval(c,zg,phig),[M*N 1]); % grid->colvec
 fastapplyATA = @(c) fastapplyAT(fastapplyA(c),L,P,zg,wg);
 matvec = @(c) fastapplyATA(c) + o.regparam*c;              % add multiple of Id
 tic;
 [cnm,flag,res,iter,resvec] = pcg(matvec,fastapplyAT(b,L,P,zg,wg),o.tol,o.maxit);
 if flag>0, warning('CG failed to converge: flag=%d\n',flag), end
 if o.verb, fprintf('CG done, fast applying A^*A, in %.3g s, rel linsys resid %.3g in %d iters\n',toc,res,iter), end
 
elseif strcmp(o.meth,'l')
  cnm = A\b;

elseif strcmp(o.meth,'n')
  applyATA = @(x) A'*(A*x);
  tic; [cnm flag res iter resvec] = pcg(applyATA,A'*b,o.tol,o.maxit);
  if flag>0, warning('CG failed to converge: flag=%d\n',flag), end
  if o.verb, fprintf('CG done, dense apply A^*A, in %.3g s, rel resid %.3g in %d iters\n',toc,res,iter), end

elseif strcmp(o.meth,'d')
  disp('meth=d is OBSOLETE; will be inaccurate')
  %figure; r=sqrt(1-z.^2); plot3(r.*cos(phi),r.*sin(phi),z,'.'); axis equal vis3d
  r=sqrt(1-z.^2); X = [r.*cos(phi),r.*sin(phi),z]'; % data pts as cols in R3
  M = ceil(3.0*P); N = round(M/2);   % # z and phi nodes (can differ)
  [zg wg] = gauss(M);
  phig = 2*pi*(0:N-1)/N; % first var is z (fast), second is phi (slow)
  [pphi zz] = meshgrid(phig,zg);
  rr = sqrt(1-zz.^2); G = [rr(:).*cos(pphi(:)),rr(:).*sin(pphi(:)),zz(:)]';
  dist = min(0.05, 0.6*pi/M)  % discs won't use all the data (cover S2)
  s = setupnearpts(X,dist);
  uu = 0*zz;         % grid data
  for j=1:numel(uu)
    i = getnearpts(X,s,G(:,j),dist);
    if numel(i)>0, uu(j) = mean(b(i)); end   % otherwise NaN
  end
  cnm = flattencnm(spharmproj(reshape(uu, [M N]), zg, wg, P));

end
%%%%

function test_lsqsolvespharm
P = 30;     % max degree (enough to capture the below f to 1e-6)
m = 1e6;    % # pts on S^2
z = rand(m,1)*2-1; phi = 2*pi*rand(m,1);   % iid rand pts on S2
% a blob function (width & placement is crucial for cnm coeff decay):
f = @(z,phi) exp(-0.5*((z-.2).^2+(phi-1.3).^2)/0.15^2);
[cnme err] = spharmprojfunc(f,P,'show'); drawnow
fprintf('P=%d func proj err %.3g\n',P,err)
fprintf('eval m=%d data points... ',m)
b = f(z,phi);    % noise-free data vector
%b = b + 10*randn(size(b));  % add noise?
fprintf('\nLSQ solve (wait a few secs)... ')
%o.meth = 'n';
o.verb = 1; cnm = lsqsolvespharm(b, z, phi, P, o);
fprintf('rel l2 coeff err %.3g\n',norm(cnm-cnme)/norm(cnme))

% timing: P=100, m=1e7: 4 mins (varies 100% to 800% CPU), ~20 GB RAM
