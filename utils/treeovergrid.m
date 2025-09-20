function [E, varargout] = treeovergrid( N, varargin )
% treeovergrid returns a tree of N nodes.
% An undirected edge set is returned. The tree is obtained by connecting chains
% over a rectangular grid of N nodes.
%
% [E]=treeovergrid(N) returns a randomly selected edge set E. The tree E represents
% is undirected, i.e. if (i,j) is in E, then so is (j,i).
%
% [E]=gausstreeovergrid(N, 'row' ) returns E with row-wise chains.
% 
% [E]=gausstreeovergrid(N,  'col' ) returns E with coloumn-wise chains.
%

% Murat Uney 03.2024


isRow = 0;
isCol = 0;

nvarargin = length(varargin);
if nvarargin>=1
    if sum( strcmp(varargin{1},{'coloumn','col','column'}) )
        isRow = 0;
        isCol = 1;
    elseif sum( strcmp(varargin{1},{'row','rows'}) )
        isRow = 1;
        isCol = 0;
    else
        error('Unknown inputs. Type help gausstreeovergrid');
    end
end


V = [1:N]';
% Number of chains:
C = floor(sqrt(N));
M = floor(N/C); % number of nodes in a chain
if isRow
    inds = [1:N]';
elseif isCol
    inds = [reshape( [1:C*M], C, M )]';
    inds = [inds(:);[C*M+1:N] ];
else
    inds = randperm(N);
end
V = V(inds);

CV = {}; % Chains
ConnectorVs = []; % Connector vertices
for ccnt=1:C
    CV{ccnt} = V([(ccnt-1)*M+1:ccnt*M]);
    if isRow || isCol
        connectorInd = round( (ccnt-1)*M+M/2 );
    else
        connectorInd = (ccnt-1)*M+randi(M);
    end
    ConnectorVs(ccnt) = V( connectorInd );
end

if C*M<N
    CV{end+1}=V([C*M+1:N]);
    if isRow || isCol
        connectorInd = floor(C*M+(N-CM)/2);
    else
        connectorInd = C*M + randi(N-CM);
    end
    ConnectorVs(end+1) = V( connectorInd );
end

E = [];

for ccnt=1:numel(CV)
    % On the path from the connector node
    chain = CV{ccnt};
    connode = ConnectorVs(ccnt);
    conind = find( chain==connode );

    for vcnt=1:conind-1
        e = [chain(vcnt), chain(vcnt+1)];
        re = [chain(vcnt+1), chain(vcnt)];
        E = [E;e;re];
        
    
    end

    for vcnt=conind:length( chain )-1
        e = [chain(vcnt+1), chain(vcnt)];
        re = [chain(vcnt), chain(vcnt+1) ];
        E = [E;re;e];
        
      
    end
end

% Add now the connector chain
for vcnt = 1:length( ConnectorVs )-1
    e = [ConnectorVs(vcnt),ConnectorVs(vcnt+1)];
    re = [ConnectorVs(vcnt+1), ConnectorVs(vcnt)];
    E = [E;e;re];
end

if nargout>=2
    varargout{1}= inds;
end
