function [J, varargout] = gcicost( gi, w )
% FUNCTION GCICOST is the weighted Kullback Leibler Divergence sum 
% to maximise for Generalised Covariance
% Intersection

% Murat Uney May 2021

N = length( gi );
gf = emd( gi, w );

J = 0;
for cnt=1:N
    J = J + w(cnt)*kld( gf, gi(cnt) );
end

if nargout>=2
    varargout{1} = gf;
end


