function [ kvals , varargout] = kldfromemd( gi, w )


% Murat Uney May 2021

N = length( gi );
gf = emd( gi(:), w(:) );

kvals = zeros(N,1);
for cnt=1:N
    kvals(cnt) = kld( gf, gi(cnt) );
end
if nargout>=2
    varargout{1} = gf;
end


