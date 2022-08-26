function showsphgrid(u,z,sc)
% SHOWSPHGRID  color plot of a function sampled on tensor grid on sphere
%
% showsphgrid(u,z,sc)
%  u = M*N array of real or complex (in which case Re is taken) values.
%  z = size-M vector of z = cos(theta) coordinates corresponding to rows of u.
%  sc (optional) = color scale (value of u that maxes out the color range).
%                  Color range is always symmetric about zero.
%
% With no arguments, a self-test is done.
%
% Notes: 1) phi grid is now 0-indexed, as is projloc3d.
% 2) the only amusing content here is making the sphere plot seal up closed.

% Barnett 8/1/15
% 0-indexed phi grid assumed for u, 9/1/15
if nargin==0, test_showsphgrid; return; end
[M N] = size(u); if numel(z)~=M, error('z length must match size(u,1)!'); end
z = z(:);  % ensure col

% glue up the top and bottom holes with average data:
if z(end)>z(1), z = [-1;z;1]; else z = [1;z;-1]; end
u = [mean(u(2,:))*ones(1,N); u; mean(u(end-1,:))*ones(1,N)];
u = [u u(:,1)];   % append in phi by one column to close up at phi=0 (0-indexed)

rho = sqrt(1-z.^2);      % col
rr = repmat(rho,[1 N+1]);
phi = 2*pi*(0:N)/N;   % row (includes origin twice for closure)
xx = repmat(cos(phi),[M+2 1]) .* rr;
yy = repmat(sin(phi),[M+2 1]) .* rr;
zz = repmat(z,[1 N+1]);
ushow = real(u);   % what to show
%figure(1); 
surf(xx,yy,zz,'cdata', ushow); shading interp
axis vis3d; 
%xlabel('x'); ylabel('y'); zlabel('z'); 
grid off
if nargin<3, sc = max(ushow(:)); end
if sc>0, caxis(sc*[-1 1]); end
colorbar; view(135,30); light; %zoom(1.5)  % Alex's defaults

%%%%%
function test_showsphgrid
M = 20; N = 30;
z = gauss(M);
ph =  2*pi*(0:N-1)/N;
deg = 4; ord = 3; ym = legendre(deg,z); y = ym(ord+1,:);
u = bsxfun(@times, y(:), exp(1i*ord*ph(:)'));
max(real(u(:)))
figure; showsphgrid(u,z)
%showsphgrid(u,z,1); title('sc=1')
