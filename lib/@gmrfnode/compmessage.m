function m = compmessage(this, l )

% Find the local child node index number for node l
s = find( this.children == l );
if isempty(s)
    warning(sprintf('Node %d is not a child of node %d \n returning empty message',s, this.id));
    m = [];
    return;
end

dt = this.nodepot.getdims; % dimensionality of x_t

% Combine the incoming messages
indeg = length( this.parents); % degree of the node
if isa( this.previnbox, 'gpot' ) && sum( this.previnboxlog ) == indeg
    
    sumind = setdiff( [1:indeg], find( this.parents == l ) );
       
    Lambda_sum = zeros(dt, dt);
    nu_sum = zeros( dt, 1);
    for i=1:length( sumind )
        Lambda_sum = Lambda_sum + this.previnbox( sumind(i) ).Lambda;
        nu_sum = nu_sum + this.previnbox( sumind(i) ).nu;
    end


    if ~isempty(this.y)
        Mmat = this.obsmat'*this.Lambdan*this.obsmat; % Measurement matrix
        mterm = this.obsmat'*this.Lambdan*this.y; % Measurement term
    else
        Mmat = zeros( dt, dt )+ this.nodepot.Lambda;
        mterm = zeros( dt, 1 )+ this.nodepot.nu;
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
   % if Lambda_message<0
   %     warning(sprintf('Negative Lambda_message evaluating at %g', Lambda_message ))
   % end
    nu_message = this.edgepotentials(s).nu(dt+1:end) - J_st*invLambda1*( mterm + nu_sum + this.edgepotentials(s).nu(1:dt) );

    

else
    % First messages
    if ~isempty(this.y)
        Mmat = this.obsmat'*this.Lambdan*this.obsmat; % Measurement matrix
        mterm = this.obsmat'*this.Lambdan*this.y; % Measurement term
    else
        Mmat = zeros( dt, dt )+ this.nodepot.Lambda;
        mterm = zeros( dt, 1 ) + this.nodepot.nu;
    end

    ds = this.edgepotentials(s).dims(2);
    % [Lambda_message, nu_message ] = gaussmarginf( this.edgepotentials(s).Lambda, [dt+1:dt+ds], this.edgepotentials(s).nu );
    J_ts = this.edgepotentials(s).Lambda(1:dt, dt+1:end);
    J_st = this.edgepotentials(s).Lambda(dt+1:end,1:dt);

    J_tofs = this.edgepotentials(s).Lambda(1:dt, 1:dt );
    J_soft = this.edgepotentials(s).Lambda(dt+1:end,dt+1:end);

    Lambda1 = J_tofs + Mmat; % 
    invLambda1 = inv(Lambda1);
    Lambda_message = J_soft - J_st*invLambda1*J_ts;
   % if Lambda_message<0
   %     warning(sprintf('Negative Lambda_message evaluating at %g', Lambda_message ))
   % end
    nu_message = this.edgepotentials(s).nu(dt+1:end) - J_st*invLambda1*( mterm + this.edgepotentials(s).nu(1:dt) );
end

m = gpot;


m.Lambda = Lambda_message;
m.nu = nu_message;
m.dims = ds;
m.numvariables = 1;
