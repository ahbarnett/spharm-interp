function L = fillsparseL(zg,phig,z,phi,q,o)
% L = fillsparseL(zg,phig,z,phi,q) returns sparse interpolation matrix from
%  tensor product grid on S^2 to arbitrary points on S^2 given
%  by the length-m vectors z and phi. zg is a length-M list of z=cos(theta)
%  nodes that must be a length-M type-1 z-increasing Fejer grid
%  (z_n=cos(theta_n) where theta_n=pi*(M-n-.5)/M), while phig is a
%  length-N list of azimuths that must be uniform, but may be 0- or 1-indexed.
%  q sets the # of Lagrange points used for interpolation. Effort O(q^2 m).
%  The ordering of the MN columns of L is such that zg is "fast" and phig "slow".
%
%  Note: if zg is not close to a z-increasing type-1 Fejer grid, error is raised.
%
% No arguments does self-test.
%
% Barnett 8/19/15.
% 8/20/15 fixed small odd order m sqrt(1-z^2) singularity prob via Chebyshev.
% 8/21/15 enforced Fejer grid, interp done in -theta, ghost nodes through poles.

% input checking
if nargin==0, test_fillsparseL; return; end
m = numel(z); z = z(:);       % # output pts, make col vec
if numel(phi)~=m, error('number of phi and z values must match'); end
if mod(q,2)==0, q=q+1; warning('even q is being bumped up by 1 to be odd'); end
M = numel(zg);
if max(abs(zg(:)'-cos(pi*(M-.5:-1:.5)/M)))>1e-14
  error('zg not close to increasing type-1 Fejer grid: cannot continue')
end

% common 1d interp setup
hwid = round((q-1)/2);   % half-width in indices
x = -hwid:hwid;          % rescaled interp grid centered around 0
w = baryweights(x);

% set up phi 1d interp from regular grid (see toy1dfseries/fillsparseL) ...
phi = phi(:)-phig(1);   % make as if phig 0-indexed, and col vec
N = numel(phig);
dphi = 2*pi/N; invdphi = 1/dphi;          % phi reg grid spacing
iphi = round(phi*invdphi); tphi = phi*invdphi - iphi; % tphi = frac part
iphi = iphi-(hwid-1);     % iphi = start indices (phi nodes) for each target
Sphi = baryprojs(x,w,tphi);             % interp weights as m*q array

% set up t=-theta 1d interp from regular grid... (still 1e-14 acos acc nr poles)
tg = -acos(zg); t = -acos(z);   % use t=-theta since incr from -pi to 0 as z incr
dt = pi/M; invdt = 1/dt;        % t reg grid spacing
it = ceil((t+pi)*invdt);        % center indices (t nodes) for each target
it(find(it<1))=1; it(find(it>M))=M; % falling off ends
tt = (t+pi)*invdt - it + 1/2;  % tt = frac part relative to nodes (ie in x space)
St = baryprojs(x,w,tt);             % interp weights as m*q array

iL = zeros(q*q,m); jL = iL; valL = iL;  % alloc lists sparse mat, as rect arrays
c = 1;  % counter for target pts
for i=1:M    % ...... loop over t subintervals, finding targets in each
  ii = find(it==i); n = numel(ii);   % target pt indices, O(m) - todo: speed up?
  jt = i+(-hwid:hwid);               % local t node indices
  flipphi = jt<1 | jt>M;             % ghost pts flip phi by pi passing thru pole
  bot = jt<1; jt(bot) = 1-jt(bot);   % fell off bottom: ghost pts, reflect j back
  top = jt>M; jt(top) = 2*M+1-jt(top);  % "     top
  % vectorized creation of col indices, which are  z_index + M*(phi_index-1)...
  js = bsxfun(@plus, iphi(ii)', (0:q-1)');   % blk of phi inds, size q*n
  js = reshape(kron(js,ones(q,1)) + kron(ones(q,n),(N/2)*flipphi'),[q q n]); %yuk
  js = M*mod(js-1,N);         % wrap phi inds and scale by M so can add to z inds
  indblks = reshape(repmat(jt(:),[1 q n]) + js, [q*q n]);  % sum t & phi inds
  jL(:,c+(0:n-1)) = indblks;             % set all col inds together
  for k=1:numel(ii), iik = ii(k);        % loop over all targets in subinterval
    iL(:,c) = iik;                       % all row inds are same: iik
    valblk = St(iik,:).' * Sphi(iik,:);      % outer prod q*q w/ t fast, phi slow
    valL(:,c) = valblk(:);
    c = c+1;
  end
end         % ......
if c~=m+1, warning('problem: looks like not all target points included!'); end

L = sparse(iL(:),jL(:),valL(:),m,N*M);   % build it

%%%%
function test_fillsparseL   % test errors for bandlimited func of various orders
verb = 1;
%f = @(z,phi) exp(-0.5*((z-.2).^2+(phi-1.3).^2)/0.15^2);  % known func on S^2
% (which can be rep by P=30 expansion to err 1e-6)
pn = 30; fprintf('overall P=%d\n',pn)
ppw = 10;
N = round(ppw*pn); %150;  % >2P where P is bandlimit of f
M = round(N/2);
[zg wg] = fejer(M);    % only allowed z quadr scheme
phig =  2*pi*(1:N)/N;   % 1-indexed, but 0- also works
[pphi zz] = meshgrid(phig,zg);
m = 1e5;  % # pts on S^2
z = rand(m,1)*2-1; phi = 2*pi*rand(m,1);   % rand pts on S2
q = 7;
%profile clear; profile on   % figure out why fill L is slow
tic; L = fillsparseL(zg,phig,z,phi,q);
fprintf('q=%d sparse fill %d targets, grid M=%d N=%d, in %.3g s\n',q,m,M,N,toc)
%profile off; profile viewer

for pm = [0:5 pn]  % test orders: small odd m has worst end-singularities!
  f = @(z,phi) legendrenm(pn,pm,z) .* exp(1i*pm*phi); % a sph harm
  fg = f(zz,pphi);  fg = fg(:);      % func on grid, as col vec
  tic; fi = L*fg; tap = toc;         % apply the interpolation
  fe = f(z,phi);    % exact at targets
  err = fi(:)-fe(:);
  fprintf('Pnm order m=%d: rel sup interp err %.3g, time %.3g s\n',pm,max(abs(err))/max(abs(fe(:))),tap)
  if verb, figure; semilogy(z,abs(err),'+'); xlabel('z'); ylabel('err');
    axis([-1 1 1e-16 1]); hold on; semilogy(zg,1+0*zg,'r.'); % interp nodes in z
    title(sprintf('m=%d: fastapplyA using L, target pt errs vs z',pm))
  end
  if verb>1, figure; plot(z,[fe(:) fi(:)],'+'); xlabel('z'); ylabel('f');
    hold on; plot(zg,0*zg,'r.');
    title(sprintf('m=%d: func & its interpolant evals',pm))
  end
end
