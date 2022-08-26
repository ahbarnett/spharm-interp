function u = legendrenm(n,m,t)
% u = legendrenm(n,m,t) returns single normalized assoc
% Legendre poly P_n^m(t), just for testing purposes.
% m is in range -n,...,n.
if abs(m)>n, error('need |m|<n'); end
u = legendre(n,t,'norm');
u = reshape(u(1+abs(m),:),size(t));
