function  varargout = updatestateintrbp(this, E )


% Use the messages in the inbox if not empty
% degree of the node
parents = pa(E, this.id );
indeg = length( parents );

% Find pointers to parents in E in the inbolx
pptrs = [];
for k=1:indeg
    pptrs(k) = find( parents(k)== this.parents ); 
end

inboxlogonE = this.inboxlog( pptrs );

if isa( this.inbox, 'gpot' ) && sum( inboxlogonE ) == indeg
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

    % Use all incoming messages to update the state
    

    % First, add the messages on the tree
    for i=1:indeg
        Lambda_sum = Lambda_sum + this.inbox( pptrs(i) ).Lambda;
        nu_sum = nu_sum + this.inbox( pptrs(i) ).nu;
    end
    
    % % Second, add the others not on E 
    % numallmessages = sum( this.inboxlog );
    % ind2others = setdiff([1:numallmessages], pptrs );
    % 
    % for i=1:numel(ind2others)
    %     if ~isempty(this.inbox(ind2others(i)).Lambda)
    %         Lambda_sum = Lambda_sum + this.inbox(ind2others(i)).Lambda;
    %         nu_sum = nu_sum + this.inbox(ind2others(i)).nu;
    %     end
    % end

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
this.inbox(pptrs) = gpot( size(pptrs) );
this.inboxlog(pptrs) = 0;


if nargout == 0
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end
else
    varargout{1} = this;
end

 