function [g, varargout ] = delgciObjective( finp, w )

% Murat Uney August 2021
N = length( w );
fw = emd( finp, w );

g = zeros(  N, 1 );

for cnti=1:N
    fi = finp(cnti);
    
    kldval = kld( fw, fi );
    
% Below are the terms from Eq. (45) in the Fusion'21 paper
g(cnti) = kldval + ( 1 - sum(w) )*...
    (  crossentropy( fw, fi )*entropy(fw) - capitalG( fw, fi )  );
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


