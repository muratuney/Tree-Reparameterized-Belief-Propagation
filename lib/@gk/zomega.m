function z = zomega( gi, varargin )


if length( varargin ) == 2
    % Backward compatibility
    z = zomegaold( gi, varargin{1}, varargin{2} );
elseif length( varargin ) == 1
    
    w = varargin{1};
    [J, gf] = gcicost( gi, w );
    
    z = exp(-1*(  J - (1 - sum( w ) )*entropy( gf )    ) );
    
elseif length( varargin ) == 0
    N=length( gi );
    w = ones(N,1)/N;
    [J, gf] = gcicost( gi(:), w );
    
    z = exp(-1*(  J - (1 - sum( w ) )*entropy( gf )    ) );
    
else
    error('Incorrect input');
end
end
function z = zomegaold( pi, pj, oms )
% function z = zomega( pi, pj, oms )
% Returns the following integral 
% z = \int pi(x)^(1-w)*pj(x)^w dx
% for an array of w values in oms
% 
% See also ZOMEGA, EMD

% Murat Uney 20/08/2018

m1 = pi.m;
C1 = pi.C;

m2 = pj.m;
C2 = pj.C;

d = length(m1(:));

S1 = inv(C1);
S2 = inv(C2);

for cnt=1:length( oms )
    om = oms( cnt );
    
    fact1 = ( ( (2*pi)^(d/2) ) *det(C1)^((1-om)/2) *det(C2)^(om/2) )^-1;
    
    A = (1-om)*S1 + om*S2;
    B = (1-om)*m1'*S1 + om*m2'*S2 + (1-om)*m1'*S1' + om*m2'*S2';
    B = 0.5*B';
    
    fact2 = sqrt( (2*pi)^d/det(A) )*...
        exp( -0.5*(  (1-om)*m1'*S1*m1 + om*m2'*S2*m2  ) ...
        + 0.5*B'*inv(A)*B  );
    
    z(cnt) = fact1*fact2;
end

z = reshape( z, size(oms) );
end