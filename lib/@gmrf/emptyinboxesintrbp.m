function varargout = emptyinboxesintrbp( this, E )


global DEBUG_GMRF DEBUG_VERBOSE

if exist('DEBUG_VERBOSE')
    if DEBUG_VERBOSE
        disp('Inside @gmrf.emptyinboxesintrbp');
    end
end


this.N = size( this.V, 1 );
% Traverse through the node objects
for i=1:this.N
    nodeobj = this.nodes(i);

    parents = pa(E, nodeobj.id );
    deg = length( parents );

    % Find pointers to parents in E in the inbox
    pptrs = [];
    for k=1:deg
        pptrs(k) = find( parents(k)== nodeobj.parents );
    end
   
    nodeobj.previnbox(pptrs) = gpot( size(pptrs) );
    nodeobj.previnboxlog(pptrs) = 0;
    nodeobj.inbox(pptrs) = gpot( size(pptrs) );
    nodeobj.inboxlog(pptrs) = 0;

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
