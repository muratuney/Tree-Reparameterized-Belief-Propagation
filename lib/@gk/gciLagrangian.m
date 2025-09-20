function [L, varargout] = gciLagrangian( gi, x )

% Murat Uney May 2021
Np1 = length( x );
N = Np1 - 1;

w = x(1:N);
lambda = x(Np1);

gf = emd( gi, w );
J = 0;
for cnt=1:N
    J = J + w(cnt)*kld( gf, gi(cnt) );
end

L = J + lambda*( 1 - sum(w) );
if nargout>=2
    varargout{1} = gf;
end


