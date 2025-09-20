function [p, varargout] = pa( G, inp )
% p = pa( E, i )
% returns the parent nodes of the node labelled i with respect to the Edge
% array E (2xN). 
% [p, e ] = pa( E, i)
% returns the list of the edges towards i from the parents in e.
% [p, e ] = pa( G, i) is also a valid function call where G is a cell
% array representing the graph G = {V,E} with E being the list of vertices.
% [p, e ] = pa( E, N) returns the intersection of all parents of nodes in
% N with respect to E.

% Updated 04.2024 to return the intersection of all parents of a set of
% nodes
% Updated 20.01.2014 Murat Uney
% First version 2009 Murat Uney

p = [];
enums = [];

for cnt=1:length(inp)
    i = inp(cnt);
    if iscell( G )
        enums_ = find( G{2}(:,2)== i);
        p_ = unique( G{2}( enums_, 1),'stable');
    else
        enums_ =  find( G(:,2)== i);
        p_ = unique( G( enums_, 1),'stable');
    end
    p = [p;p_];
    enums = [enums;enums_];
end
p = unique(p, 'stable' );
enums = unique(enums, 'stable' );


if nargout>1
    varargout{1} = enums;
end
    
