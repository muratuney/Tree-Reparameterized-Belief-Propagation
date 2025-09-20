function [npe, epe] = subgraphpotentials( np, ep, E )
% function subgraphpotentials finds the node and edge potentials of a 
% pairwise Markov Random Field given its potentials.
%
% [npe, epe] = subgraphpotentials( np, ep, E ) returns the node and edge
% potentials over E, given all node and edge potentials of the MRF in gpot
% array np and ep

% August 2024 Murat Uney
V = unique( E(:));
Ls = zeros( numel(V) );

np = np(:);
ep = ep(:);

npe = gpot([]);
epe = gpot([]);

% Assign the diagonal entries of Ls
for cnt=1:length( V )
    for cnt2 = 1:length(np)
        if np(cnt2).ids(1) == V(cnt)
            npe(end+1) = np(cnt2);
        end
    end
end

for cnt=1:size( E,1 )
    for cnt2=1:length( ep )
        if ep(cnt2).ids(1)== E(cnt,1) && ep(cnt2).ids(2)==E(cnt,2) 
            epe(end+1) = ep(cnt2);
        end
    end
end