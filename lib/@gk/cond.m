function [C_pGq, a, B] = cond( p, varargin )
% returns the parameters (covariance matrix and the mean vector)
% of the required conditional distribution given the those of the joint.
%
%   [ C, a, B ] = p.cond() returns the parameters; covariance C, vector a 
% and transform B where the mean is mu = a + B x(2:N) for 
% p(x1 | x2 ... xN )= N(x1 ; C, mu ).
%
%   [ C, a, B ] = p.cond( p_ind ) returns the 
%   matrix of p( x_p ) where x_p = [x_j]_{j \in p_ind} considering the joint 
%   N(x;  p.C, p.m  ).
%

% Murat Uney

N = length(p.m);

if nargin >= 2
    p_ind = unique( varargin{1}(:) , 'stable');
    if ~isnumeric(p_ind)
        error('The second argument should be an index array of type numeric');
    end
    if ~isempty( find( p_ind<1 | p_ind>N ) ) | length(p_ind)>=N
        error('The index array contains entries that exceeds 1,...,N');
    end
    lenQ = length( p_ind );
else
    p_ind = [2:N]';
    lenQ = N-1;
end

q_ind = setdiff([1:N], p_ind );
[C_pGq, ~, a, B] = gausscond( p.C, q_ind, p.m );
