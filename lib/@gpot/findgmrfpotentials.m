function [ np, ep ] = findgmrfpotentials( m, C, E, varargin  )
% findgmrfpotentials returns the edge potentials of a Gaussian density over
% a selected graph with edges E
%
% The variables associated with the vertices are assumed to be ordered. The
% default dimensionality is 1.
%
%  [ np, ep ] = findgmrfpotentials( m, C, E ) returns the node potentials
%  in np and edge potentials in ep for the density with mean vector m and
%  covariance C over the graph with edges E.
%
% [ np, ep ] = findgmrfpotentials( m, C, E, d ) uses the scalar d as the
% node dimensionality.
%
% [ np, ep ] = findgmrfpotentials( m, C, E, D ) uses the array D to read
% the node dimensionality from each entry.

% These potentials should be used in accordance with Eq.s 2.42-2.45 in Erik
% Sudderth's master thesis.
%
% Murat Uney 01.2014 Initial implementation
% Murat Uney 03.2024 Re-implementation as part of @gpot
% Upgraded to evenly split the potentials between node and edge potentials,
% so 2.42-2.45 do not directly apply - they apply to the combination of the
% edge and node potentials over the same pair of nodes.


d = 1;


d = ones( numel(m), 1 );
if nargin>2
    inp = varargin{1}(:);
    if numel( inp ) == 1
        d =  ones( floor( numel(m)/inp(1) ), 1 );
    else
        d = inp;
    end
end

% Find the inverse covariance
R = chol(C);
Rinv = R^(-1);
Lambda = Rinv*Rinv';

nu = Lambda*m;

[ np, ep ] = gpot.findinfogmrfpotentials( nu, Lambda, E, d  ); % This function is implemented in a separate file
      

