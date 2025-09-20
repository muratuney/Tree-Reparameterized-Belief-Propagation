function D = findoutdegree( E, V )
% FINDOUTDEGREE returns the in-degree of the vertices in an undirected graph
% D = findoutdegree(E,V) returns the out-degree (i.e. the number of descendants)
% of vertices V in the undirected graph with the edge set E.
%

% Murat Uney

D = zeros(size(V));
for i=1:length(V)
    D(i) = length( find(E(:,1)==V(i)) );
end
    