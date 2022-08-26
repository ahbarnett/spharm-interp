% test projection of spherical harmonics on S2, interface to Fortran
% Barnett 7/31/15-8/1/15
% changed to 0-indexing & rect case 9/1/15

clear; verb = 1;  % 0= just text, 1 = u fig, 2 = more figs
n = 15; m = 30;  % # nodes in z, and in phi
[z w] = gauss(n);
ph = 2*pi*(0:m-1)/m;    % first var is z (fast), second is phi (slow)

P = 10;         % max degree. somewhat less than 0.75 n for good accuracy
u = ones(n,m);  % should project just to c(0,0), but with what coeff? 4pi? 1!
cnm = spharmproj(u, z, w, P);
cnm(1,P+1) = cnm(1,P+1)-1; fprintf('proj err for u=1: %.3g\n',norm(cnm(:)))
if verb>1, figure; imagesc(-P:P,0:P,abs(cnm)); colorbar; end

% guesses normalization & phase of Matlab's legendre rel to FMM:
deg = 5; ord = -4;   % can replace w/ anything up to P
ym = sqrt(2)*legendre(deg,z,'norm'); y = (-1)^ord*ym(abs(ord)+1,:);
u = bsxfun(@times, y(:), exp(1i*ord*ph(:)'));  % should excite just c(deg,ord)
if verb, figure; showsphgrid(u,z); end
cnm = spharmproj(u, z, w, P);
fprintf('1-|coeff|^2 err %.3g\n',abs(1-abs(cnm(1+deg,P+1+ord))^2))
fprintf('phase angle of coeff = %.15g pi\n',angle(cnm(1+deg,P+1+ord))/pi)
cnm(1+deg,P+1+ord) = 0;
fprintf('proj err for rest of coeffs for (%d,%d): %.3g\n',deg,ord,norm(cnm(:)))
if verb>1, figure; imagesc(-P:P,0:P,abs(cnm)); colorbar; hold on; plot(ord,deg,'*'); end
