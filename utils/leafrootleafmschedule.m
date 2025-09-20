function M = leafrootleafmschedule(E)
% leafrootleafmschedule returns a leaf-root-leaf message schedule for an
% undirected tree. 
%
% M = leafrootleafmschedule(E) returns the message schedule in a 1xt cell 
% array M. Each entry is an array of directed edges indicating message directions.
% E is an undirected edge set, i.e. if (i,j) is in E, so is (j,i).

% Murat Uney

% Find the vertices V
V = sort(unique(E(:),'legacy'));
% Find their out-degrees
D = findoutdegree( E, V );
allDs = sort( unique( D(:),'legacy'));

M = {};
mpat = []; % MessagePattern

Ecopy = E;
Vsource = V;
Dsource = D;
% Leaf-root messaging
prevsinknodes = [];
prevsourcenodes = [];
cnt = 1;
% Find the leaves
dinds = find( Dsource == 1 );
% Pick the leaves
leafnodes = Vsource(dinds); % Pick the leaves
% Start from the leaves and traverse upwards in the tree
sourcenodes = leafnodes;
sinknodes = chi( E, sourcenodes );
while(~isempty(setdiff( sinknodes, prevsourcenodes )) )
    
    mpat = [];
    newsinknodes = [];
    for scnt = 1:length(sourcenodes)
        % Find the sink nodes
        inds = find( E( : , 1) == sourcenodes(scnt) ); % all edges to sink nodes indices
        [newsn, indns] = setdiff( E(inds,2), prevsourcenodes ); % The sink nodes should not be previous source nodes
       
        newEdges = Ecopy(inds(indns), : ) ; % New edges in the message pattern
        mpat = [mpat;newEdges];       

        newsinknodes = [newsinknodes; newEdges( : , 2)];
    end
    uniquesinknodes = unique( newsinknodes ,'legacy');
    prevsourcenodes = [prevsourcenodes;sourcenodes];
    prevsinknodes = [prevsinknodes; uniquesinknodes];


    M{cnt} = mpat;
    cnt = cnt+1;

    % Prune the sink nodes that have not received all messages but one
    sourcenodes = prune( uniquesinknodes, E, M  );
    sinknodes = chi( E, sourcenodes );
end

% Now, move back in the tree
% for root-leaf messaging
lastsinknodes = unique( M{end}(:,2),'legacy');
prevsinknodes = [];
prevsourcenodes = [];

cnt = 1;

sourcenodes = lastsinknodes;
sinknodes = pa( E, sourcenodes );
Md = {};
cnt = length(M) +1;
mpat = []; % MessagePattern
while(~isempty(setdiff( sinknodes, prevsourcenodes )) )
    
    mpat = [];
    newsinknodes = [];
    for scnt = 1:length(sourcenodes)
        % Find the sink nodes
        inds = find( E( : , 1) == sourcenodes(scnt) ); % all edges to sink nodes indices
        [newsn, indns] = setdiff( E(inds,2), prevsourcenodes ); % The sink nodes should not be previous source nodes
       
        newEdges = Ecopy(inds(indns), : ) ; % New edges in the message pattern
        mpat = [mpat;newEdges];       

        newsinknodes = [newsinknodes; newEdges( : , 2)];
    end
    uniquesinknodes = unique( newsinknodes ,'legacy');
    prevsourcenodes = [prevsourcenodes;sourcenodes];
    prevsinknodes = [prevsinknodes; uniquesinknodes];

    M{cnt} = mpat;
    cnt = cnt+1;

    % Prune the sink nodes that have not received all messages
    sourcenodes = prune2( setdiff( uniquesinknodes, leafnodes ), E, M );
    sinknodes = pa( E, sourcenodes );
end


% M = [M,Md];


end

function out = prune( inp, E, M )
E2 = cell2mat(M(:));
out = [];
for cnt=1:length(inp)
    if findindegree( E, inp(cnt) ) == findindegree( E2, inp(cnt) )+1
        out = [out; inp(cnt)];
    end
end
end
function out = prune2( inp, E, M )
E2 = cell2mat(M(:));
out = [];
for cnt=1:length(inp)
    if findindegree( E, inp(cnt) ) == findindegree( E2, inp(cnt) )
        out = [out; inp(cnt)];
    end
end

end