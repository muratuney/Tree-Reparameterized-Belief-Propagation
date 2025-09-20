function [ch, varargout] = chi(G,inp)
% ch = chi( E, i )
% returns the child nodes of the node labelled i with respect to the Edge
% array E (2xN). 
% [ch, e ] = chi( E, i)
% returns the list of the Edges from i to the children in e.
% [ch, e ] = chi( G, i) is also a valid function call where G is a cell
% array representing the graph G = {V,E} with E being the list of vertices.
% [ch, e ] = chi( E, N) returns the intersection of all children of all 
% nodes in N with respect to E.

% Updated 04.2024 to return the intersection of all children of a set of
% nodes
% Updated 17.05.2017 unique now is called with 'stable' to have the same
% behaviour
% Updated 20.01.2014 Murat Uney
% First version 2009 Murat Uney

ch = [];
enums = [];
for cnt=1:length(inp)
    i = inp(cnt);
    if iscell( G )
        enums_ = find(G{2}(:,1)== i);
        ch_ = unique( G{2}( enums_,2),'stable');
    else
        enums_ = find(G(:,1)== i);
        ch_ = unique( G( enums_,2), 'stable');
    end
    ch = [ch;ch_];
    enums = [enums;enums_];
end
ch = unique(ch, 'stable' );
enums = unique(enums, 'stable' );

if nargout>1
    varargout{1} = enums;
end

    


