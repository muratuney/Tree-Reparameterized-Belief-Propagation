function  varargout = updatestate(this )


% Use the messages in the inbox if not empty
% degree of the node
deg = length( this.parents );
if isa( this.inbox, 'gpot' ) && sum(this.inboxlog) == deg
    % Update the state
    % The below implement equations 2.34 and 2.35 
    ds = this.nodepot.getdims;


    Lambda_sum = zeros(ds, ds);
    nu_sum = zeros( ds, 1);
    if ~isempty(this.y)
        Mmat = this.obsmat'*this.Lambdan*this.obsmat; % Measurement matrix
        mterm = this.obsmat'*this.Lambdan*this.y; % Measurement term
    else
        Mmat = zeros( ds, ds ) + this.nodepot.Lambda;
        mterm = zeros( ds, 1 ) + this.nodepot.nu;
    end

    
   
    for i=1:deg
        Lambda_sum = Lambda_sum + this.inbox(i).Lambda;
        nu_sum = nu_sum + this.inbox(i).nu;
    end
    Lambda_state = Lambda_sum + Mmat; % 2.35
    nu_state = nu_sum + mterm; % 2.34

    infostate = this.infostate;
    infostate.Lambda = Lambda_state;
    infostate.nu = nu_state;
    this.infostate = infostate;

    C_up = inv(Lambda_state);
    m_up = C_up*nu_state;
    
    try
        this.state = cpdf( gk( C_up, m_up ) );    
    catch
        warning(' Ill-conditioned inverse covariance! Could not find the pdf version of the belief state at %d', this.id)
    end 
else
    % no message to update
    disp( sprintf('Node %d does not have sufficient messages in the inbox to update its state', this.id ));
end
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

 