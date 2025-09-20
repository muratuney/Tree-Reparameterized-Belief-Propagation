function [Es, varargout] = rowtreesovergrid( N )
% rowtreesovergrid returns row-wise trees of N nodes with different connector
% edges across rows. An undirected edge set is returned. The trees are
% obtained by connecting chains over a rectangular grid of N nodes.
%
% [E]=rowtreesovergrid(N) returns a randomly selected edge set E. The tree E represents
% is undirected, i.e. if (i,j) is in E, then so is (j,i).
%

% Murat Uney 06.2024



V = [1:N]';
% Number of chains:
C = floor(sqrt(N));
M = floor(N/C); % number of nodes in a chain
inds = [1:N]';
V = V(inds);

for cnt=1:M
    Es = rowtreesovergrid( N, )


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
end
