function D = findindegree( E, V )
% FINDINDEGREE returns the in-degree of the vertices in an undirected graph
% D = findindegree(E,V) returns the in-degree (i.e. the number of ascendants)
% of vertices V in the undirected graph with the edge set E.
%

% Murat Uney

D = zeros(size(V));
for i=1:length(V)
    D(i) = length( find(E(:,2)==V(i)) );
end
    