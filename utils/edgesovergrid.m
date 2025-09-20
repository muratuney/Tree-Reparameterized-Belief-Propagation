function [E] = edgesovergrid( N )
% edgesovergrid returns a loopy grid of N nodes.n of variables are 1.
%
% [E]=edgesovergrid(N) returns the edge set E of a grid graph.
% E is an undirected edge set, i.e. if
% (i,j) is in E, then so is (j,i).
%

% Murat Uney 03.2024

V = [1:N]';
% Number of chains:
C = floor(sqrt(N));
M = floor(N/C); % number of nodes in a chain
V = [1:N];

CV = {}; % Chains
for ccnt=1:C
    CV{ccnt} = V([(ccnt-1)*M+1:ccnt*M]);
end
if C*M<N
    CV{end+1}=V([C*M+1:N]);
end

E = [];
for ccnt=1:numel(CV)
    chain = CV{ccnt};
    for vcnt=1:length( chain )-1
        e = [chain(vcnt), chain(vcnt+1)];
        re = [chain(vcnt+1), chain(vcnt)];
        E = [E;e;re];
        
    end
end

% Add now the connector chains
for ccnt=1:numel(CV)-1
    % On the path from the connector node
    chain1 = CV{ccnt};
    chain2 = CV{ccnt+1};
    for vcnt=1:min( length( chain1 ), length( chain2 )  )

        e =  [chain1(vcnt), chain2(vcnt)];
        re = [chain2(vcnt), chain1(vcnt)];
        E = [E;e;re];
        
    end
end



