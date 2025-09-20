function varargout = resetioboxes(this)

this.previnbox = gpot( [] );
this.previnboxlog = zeros(1, length(this.parents) );
this.inbox = gpot( [] );
this.inboxlog = zeros(1, length(this.parents) );


if nargout == 0
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end
else
    varargout{1} = this;
end
