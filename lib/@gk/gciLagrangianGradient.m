function [g] = gciLagrangianGradient( finp, x )

% Murat Uney May 2021
Np1 = length( x );
N = Np1 - 1;

w = x(1:N);
lambda = x(Np1);
fw = emd( finp, w );

g = zeros(  Np1, 1 );
g(Np1) = (1 - sum(w) );

for cnti=1:N
    fi = finp(cnti);
    
    kldval = kld( fw, fi );
    
% % Below is Eq. (19)
%     g(cnti) = kldval + ( 1 - sum(w) )*...
%         (  crossentropy( fw, fi ) - capitalG( fw, fi )  ) ...
%         - lambda;

% Below is Eq. (45) in the Fusion'21 paper
g(cnti) = kldval + ( 1 - sum(w) )*...
    (  crossentropy( fw, fi )*entropy(fw) - capitalG( fw, fi )  ) ...
    - lambda;
end

if nargout>=2
    varargout{1} = fw;
end
end
function g = capitalG( f1, f2 )
N = 10000;

samples = f1.gensamples( N );
ef1 = f1.evaluatelog( samples );
ef2 = f2.evaluatelog( samples );
g = (1/N)*sum( ef1.*ef2 );

end


