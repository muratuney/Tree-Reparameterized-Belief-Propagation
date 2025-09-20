function [J, varargout] = gciObjective( finp, w )
% FUNCTION GCIOBJECTIVE is the weighted Kullback Leibler Divergence sum 
% to maximise for Generalised Covariance
% Intersection

% Murat Uney May 2021

N = length( finp );
fw = emd( finp, w );

J = 0;
for cnt=1:N
    J = J + w(cnt)*kld( fw, finp(cnt) );
end

if nargout>=2
    varargout{1} = fw;
end


