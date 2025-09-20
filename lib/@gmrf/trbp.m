function varargout = trbp( this, messagetrees, varargin )
% trbp is the tree reparameterised BP method for @gmrf

global DEBUG_VERBOSE
global DEBUG_GMRF
global DEBUG_GMRF_CARRAY

if exist( 'DEBUG_GMRF' ) && exist('DEBUG_VERBOSE')
    if DEBUG_GMRF && DEBUG_VERBOSE
        disp(sprintf('Tree reparameterised BP - messages are sent in accordance with the leaf-root-leaf schedule of the specified spanning trees'));
    end
end

finum = this.iternum+1; % Start index of the iterations.
linum = this.itermax; % Last iteration numer

numberOfTrees = numel( messagetrees );
for cnt=1:numberOfTrees
    mschedule = leafrootleafmschedule( messagetrees{cnt} );
    messageschedules{cnt} = [mschedule];
end


for i=finum:linum
    % ith iteration
    if exist('DEBUG_VERBOSE')
        if DEBUG_VERBOSE
            disp(sprintf('iteration %d',i));
        end
    end
    
    treeptr = mod( (i-finum), numberOfTrees ) + 1;  % Pick the tree pointer
    Etilde = messagetrees{treeptr};                % Pick the current tree
    myschedule = messageschedules{treeptr}; % The messaging pattern is the entire set of edges

    this.emptyinboxesintrbp(Etilde); % Empty the messages in the message boxes to run message passing on this tree
    
    msStepNum = length(myschedule);
    for msstepcnt = 1:msStepNum
        mypat = myschedule{msstepcnt}; % Take the messaging pattern

        for k=1:size( mypat, 1 ) % Loop over edges in the pattern
            e = mypat(k,:);
            if exist('DEBUG_VERBOSE')
                if DEBUG_VERBOSE
                    disp(sprintf('Messaging from %d to %d', e(1), e(2) ));
                end
            end

            tnode = this.nodes( e(1) );
            rnode = this.nodes( e(2) );

            message2pass = tnode.compmessagemostrecentintrbp(Etilde, e(2) ); % Compute the outgoing message

            if isempty( message2pass )
                warning(sprintf('Message from %d to %d empty in iteration number %d', e(1), e(2), i ));
            end

            rnode = rnode.recmessage( e(1), message2pass ); % Save the received message

            this.nodes( e(1) ) = tnode;
            this.nodes( e(2) ) = rnode;
        end
        this.numberofmessages = this.numberofmessages + size( mypat, 1 );
        % Update states
        for k=1:length( this.V )
            unode = this.nodes( k );
            % Update the node states if all messages all received
            if unode.rxallmessages( Etilde )
                if DEBUG_VERBOSE && DEBUG_GMRF
                    disp(sprintf('Updating the state of node %d', k ));
                end
                unode.updatestateintrbp( Etilde );
             end
            this.nodes( k ) = unode;
        end
    end % End of messaging pattern over the tree
    if exist('DEBUG_GMRF') && exist('DEBUG_GMRF_CARRAY')
        if DEBUG_GMRF
            DEBUG_GMRF_CARRAY{end+1} = this;
        end
    end
    this.iternum = i;

    % Update edge potentials
    nexttreeptr = mod( (i+1-finum), numberOfTrees ) + 1;
    %Enext = messagetrees{nexttreeptr};
    %Eupdate = findcommonedges(Etilde,Enext);
    this = updatepotentialsontree(this, Etilde);

    
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end

end % End of the iteration loop


if nargout == 0
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end
else
    varargout{1} = this;
end
end
function this = updatepotentialsontree( this, E )
% This method is doing the tree reparameterisation by updating both edge
% and node potentials on the tree

global DEBUG_VERBOSE
global DEBUG_GMRF
global DEBUG_GMRF_CARRAY

if exist( 'DEBUG_GMRF' ) && exist('DEBUG_VERBOSE')
    if DEBUG_GMRF && DEBUG_VERBOSE
        disp(sprintf('Inside @gmrf.updatepotentialsontree - updating both edge and node potentials on the tree'));
    end
end

for i=1:size(E,1)
    e = E(i,:); % (i,j)
    sourcenode = e(1); % sourcenode = i
    sinknode = e(2);  % sinknode = j

    eptr = findedge( this.E, e ); % edgepointer

    sourceptr =  find( this.V == sourcenode );
    sinkptr = find( this.V == sinknode );

    % Update the edgepotential from the source to the sink

    % i) Fetch the last message to the source node on this edge
    inptr = find( this.nodes( sourceptr ).parents == sinknode );
    messageji = this.nodes( sourceptr ).previnbox(inptr);
   

    % ii) Fetch the last message from the sink node on this edge
    outptr = find( this.nodes( sinkptr ).parents == sourceptr );
    messageij = this.nodes( sinkptr ).previnbox(outptr);
     
  
    % iii) Update the edge potential \psi(x_i,x_j) \leftarrow \psi(x_i,x_j)/(m_ij m_ji)
    epot = this.edgepots( eptr );

    dimi = epot.dims(1);
    dimj = epot.dims(2);

    Lambdaii = epot.Lambda(1:dimi, 1:dimi);
    epot.Lambda(1:dimi,1:dimi) = ( Lambdaii )*(1-this.treeweight) +  (Lambdaii - messageji.Lambda )*(this.treeweight) ;

    Lambdajj = epot.Lambda(dimi+1:dimi+dimj, dimi+1:dimi+dimj);
    epot.Lambda(dimi+1:dimi+dimj,dimi+1:dimi+dimj) = ( Lambdajj)*(1-this.treeweight) + (Lambdajj - messageij.Lambda )*(this.treeweight) ;

    epot.nu = ( epot.nu)*(1-this.treeweight)+ (epot.nu - [messageji.nu;messageij.nu])*(this.treeweight);
    % epot.nu = epot.nu*0;

    this.edgepots( eptr ) = epot;

    
    eptr2 = find( this.nodes( sourceptr ).children == sinknode ); 
    this.nodes( sourceptr ).edgepotentials(eptr2) = epot;
end

V = unique(E(:));
for cnt=1:length( V )
    ptr = find( this.V == V(cnt) );

    nodepot = this.nodes( ptr ).nodepot;
    infostate = this.nodes( ptr ).infostate;

    nodepot.Lambda  = nodepot.Lambda*(1-this.treeweight) + infostate.Lambda*this.treeweight;
    nodepot.nu = nodepot.nu*(1-this.treeweight) +  infostate.nu*this.treeweight;
    % nodepot.nu = nodepot.nu*0;

    % this.nodepots(ptr) = nodepot; % This array is vacant
    this.nodes(ptr).nodepot = nodepot;
end




end