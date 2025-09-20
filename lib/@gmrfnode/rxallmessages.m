function r = rxallmessages(this, varargin )

r = 0;
if nargin>=2
    % This condition is for the tree-reweighted bp
    E = varargin{1};
    parents = pa( E, this.id );
    indeg = length( parents );

    % Find pointers to parents in E in the inbox
    pptrs = [];
    for k=1:indeg
        pptrs(k) = find( parents(k)== this.parents );
    end
    numreceivedmessages = sum(this.inboxlog(pptrs) );
    if indeg ==  numreceivedmessages
        r=1;
    end
else   
    % The below if for loopy-bp and bp with a tailored message schedule
    if isempty( find( this.inboxlog == 0) )
        r=1;
    end
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
