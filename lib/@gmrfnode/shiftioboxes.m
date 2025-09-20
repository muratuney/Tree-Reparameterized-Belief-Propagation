function varargout = shiftioboxes(this)

this.previnbox = this.inbox;
this.previnboxlog = this.inboxlog;
this.inbox = this.inbox([]);
this.inboxlog = this.inboxlog*0;



if nargout == 0
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end
else
    varargout{1} = this;
end
