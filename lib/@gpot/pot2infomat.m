function [Lambda, nu, varargout] = pot2infomat( np, ep )


% Find V and E
V = [];
d = [];
for cnt=1:length(np)
    V(end+1,1) = np(cnt).ids(1);
    d(end+1,1) = np(cnt).dims(1);   
end
% Variable indices are below:
stind = [1; cumsum( d )+1 ];
stind = stind(1:end-1); % start indices
endind = cumsum( d );

N = endind(end);


nu = zeros( N, 1 );
Lambda = zeros( N, N );

E = [];
for cnt=1:length(ep)
    E(cnt+1,:) = ep(cnt).ids(1:2);
end
E = sue( E ); % Remove one if both (i,j) and (j,i) is an element of E

% Find pointers in the edge potentials corresponding to edges in E
eptrs = [];
for cnt=1:size(E,1)
    e = E(cnt, : );
    for cnt2=1:numel(ep)
        if ep(cnt2).ids(1) == e(1) && ep(cnt2).ids(2) == e(2)
            eptrs(end+1) = cnt2;
            break;
        end
    end
end 

% Now, find nu by adding both the node and edge potential contributions
for cnt=1:length(V)
    stptr = stind(cnt);
    endptr = endind(cnt);
    nu( stptr:endptr ) = nu( stptr:endptr ) +  np( cnt ).nu;
end

for cnt=1:size(E,1)
    e = E(cnt,:);

    epobj = ep( eptrs(cnt) );

    ind1 = find( V == e(1) );
    if ~isempty( ind1 )
        stptr = stind(ind1);
        endptr = endind(ind1);
        nu_ = epobj.nu( 1:epobj.dims(1) );
        nu( stptr:endptr  ) = nu( stptr:endptr  ) + nu_;
    end
    ind2 = find( V == e(2) );
    if ~isempty( ind2 )
        stptr = stind(ind2);
        endptr = endind(ind2);
        nu_ = epobj.nu( epobj.dims(1)+1:epobj.dims(1)+epobj.dims(2) );
        nu( stptr:endptr  ) = nu( stptr:endptr  ) + nu_;
    end
end

% Now, find Lambda by adding both the node and edge potential contributions
for cnt=1:numel( V )
    stptr = stind(cnt);
    endptr = endind(cnt);
    Lambda( stptr:endptr, stptr:endptr ) = ...
        Lambda( stptr:endptr, stptr:endptr ) + np( cnt ).Lambda;
end

for cnt=1:size(E,1)
    e = E(cnt,:);
    % Find the edge potential for e
    ind1 = find( V == e(1) );
    ind2 = find( V == e(2) );
   
    stptr1 = stind(ind1);
    endptr1 = endind(ind1);

    stptr2 = stind(ind2);
    endptr2 = endind(ind2);

    EdgeLambda = ep( eptrs(cnt) ).Lambda;
    
    Lambda(stptr1:endptr1, stptr1:endptr1 ) = Lambda(stptr1:endptr1, stptr1:endptr1 ) + EdgeLambda(1:d(ind1),1:d(ind1));
    Lambda(stptr1:endptr1, stptr2:endptr2 ) = Lambda(stptr1:endptr1, stptr2:endptr2 ) + EdgeLambda(1:d(ind1), d(ind1)+1:d(ind1)+d(ind2));
    Lambda(stptr2:endptr2, stptr1:endptr1 ) = Lambda(stptr2:endptr2, stptr1:endptr1 ) + EdgeLambda(d(ind1)+1:d(ind1)+d(ind2),1:d(ind1) );
    Lambda(stptr2:endptr2, stptr2:endptr2 ) = Lambda(stptr2:endptr2, stptr2:endptr2 ) + EdgeLambda(d(ind1)+1:d(ind1)+d(ind2), d(ind1)+1:d(ind1)+d(ind2) );
end

if nargout>=3
    varargout{1} = E;
end
if nargout>=4
    varargout{2} = V;
end


