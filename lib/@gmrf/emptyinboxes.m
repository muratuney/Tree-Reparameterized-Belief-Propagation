function varargout = emptyinboxes( this, varargin )


global DEBUG_GMRF DEBUG_VERBOSE

if DEBUG_GMRF
    disp('inside @gmrf.emptyinboxes');
end

this.N = size( this.V, 1 );
% Traverse through the node objects
for i=1:this.N 
    nodeobj = this.nodes(i);
    nodeobj.resetioboxes;
    this.nodes(i) = nodeobj;
end

if nargout == 0
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end
else
    varargout{1} = this;
end
