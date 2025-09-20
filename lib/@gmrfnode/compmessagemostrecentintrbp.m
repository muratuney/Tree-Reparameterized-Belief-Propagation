function m = compmessagemostrecentinterbp(this, E, l )
% This method computes the message from node "this" to l over the tree E by
% using the most recent messages from the neighbours - this is needed in
% leaf-root-leaf messaging.

children = chi( E, this.id );
outdeg = length( children );
% Find pointers to children in E in the edge potentials
cptrs = [];
for k=1:outdeg
    cptrs(k) = find( children(k)== this.children );
end

% Find the local child node index number for node l
s = cptrs( find( children == l ) );
if isempty(s)
    warning(sprintf('Node %d is not a child of node %d \n returning empty message',l, this.id));
    m = [];
    return;
end

dt = this.nodepot.getdims; % dimensionality of x_t

% The node potential contribution
if ~isempty(this.y)
    Mmat = this.obsmat'*this.Lambdan*this.obsmat; % Measurement matrix
    mterm = this.obsmat'*this.Lambdan*this.y; % Measurement term
else
    Mmat = zeros( dt, dt )+ this.nodepot.Lambda;
    mterm = zeros( dt, 1 )+ this.nodepot.nu;
end

% Combine the incoming messages on the tree E
parents = pa(E, this.id );
indeg = length( parents );
% Find pointers to parents in E in the inbox
pptrs = [];
for k=1:indeg
    pptrs(k) = find( parents(k) == this.parents ); 
end

Lambda_sum = zeros(dt, dt);
nu_sum = zeros( dt, 1);
% Pick the most recent messages in the inbox
sumind = setdiff( [1:indeg], find( parents == l ) );

% Add the message contributions from the current inbox
inthecurrent = find( this.inboxlog(pptrs(sumind)) == 1);
for i=1:length( inthecurrent )
    Lambda_sum = Lambda_sum + this.inbox( pptrs(sumind(inthecurrent(i))) ).Lambda;
    nu_sum = nu_sum + this.inbox( pptrs(sumind(inthecurrent(i))) ).nu;
end

% Add the message contributions from the previous inbox
intheprevious = find( this.previnboxlog(pptrs(sumind)) == 1);
diffprevious = setdiff( intheprevious, inthecurrent);
for i=1:length(diffprevious)
    Lambda_sum = Lambda_sum + this.previnbox(  pptrs( sumind( diffprevious(i) ) ) ).Lambda;
    nu_sum = nu_sum + this.previnbox(  pptrs( sumind( diffprevious(i) ) ) ).nu;
end

ds = this.edgepotentials(s).dims(2);

% Below are 2.39 and 2.40 in ES's MSc Thesis
% except, the edgepotential here is in the [x_t, x_s] order with respect to
% 2.15, i.e. here the edgepotential parameters are
% Lambda_ts = [ J_t(s), J_t,s;
%               J_s,t , J_s(t)];
%
%
%
J_ts = this.edgepotentials(s).Lambda(1:dt, dt+1:end);
J_st = this.edgepotentials(s).Lambda(dt+1:end,1:dt);

J_tofs = this.edgepotentials(s).Lambda(1:dt, 1:dt );
J_soft = this.edgepotentials(s).Lambda(dt+1:end,dt+1:end);

Lambda1 = J_tofs + Mmat + Lambda_sum; % The paranthesis term in 2.40
invLambda1 = inv(Lambda1);
Lambda_message = J_soft - J_st*invLambda1*J_ts;
%if Lambda_message<0
%    warning(sprintf('Negative Lambda_message evaluating at %g', Lambda_message ));
%    Lambda_message = -Lambda_message;
%end
nu_message = this.edgepotentials(s).nu(dt+1:end) - J_st*invLambda1*( mterm + nu_sum + this.edgepotentials(s).nu(1:dt) );

m = gpot;

m.Lambda = Lambda_message;
m.nu = nu_message;
m.dims = ds;
m.numvariables = 1;
