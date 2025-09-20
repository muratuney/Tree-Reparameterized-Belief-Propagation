function [Lambda, E, varargout] = gausstreeovergrid( N, varargin )
% gausstreeovergrid returns a tree Gaussian Markov random field of N nodes.
% An undirected MRF is produced. The tree is obtained by connecting chains
% over a rectangular grid of N nodes. The default cannonical correlation
% coeffcient rho = 0.5 and the standard deviation of variables are 1.
%
% [L, E]=gausstreeovergrid(N) returns the information matrix (or, the
% inverse covariance) L along with the edge set E. Note that the off-diagonal
% non-zero entries of L are encoded in E. E is an undirected edge set, i.e. if
% (i,j) is in E, then so is (j,i).
%
% [L, E]=gausstreeovergrid(N, 'rho', r ) 
% [L, E]=gausstreeovergrid(N,  r ) uses the canonical correlation coefficient r
% for each edge E. 
%
% [L, E]=gausstreeovergrid(N, 'variance', v ) 
% [L, E]=gausstreeovergrid(N,  v ) uses the variance value v for each variable.

% Murat Uney 03.2024


sigmasq = 1;
rho = 0.5;
isRow = 1;

nvarargin = length(varargin);
argnum = 1;
if nvarargin>=1
    if isa( varargin{1}, 'char' )
        while argnum<=nvarargin
            if isa( lower(varargin{argnum}), 'char')
                switch lower(varargin{argnum})
                    case {'rho' }
                        if argnum + 1 <= nvarargin
                            rho = varargin{argnum+1};
                            argnum = argnum + 1;
                        end
                    case {'variance','sigmasq','var'}
                        if argnum + 1 <= nvarargin
                            sigmasq = varargin{argnum+1};
                            argnum = argnum + 1;
                        end
                    case {'coloumn','col','column'}
                        isRow = 0;
                        
                    otherwise
                        error('Wrong input string');
                end
            end
            argnum = argnum + 1;
        end
    elseif isnumeric( varargin{1} )
        if nvarargin>=1
            rho = varargin{1};
        end
        if nvarargin>=2
            sigmasq = varargin{2};
        end
        if nvarargin>=3
            if sum( strcmp(varargin{3},{'coloumn','col','column'}) )
                isRow = 0;
            end
        end
        
    else
        error('Unknown inputs. Type help gausstreeovergrid');
    end
end
sigma = sqrt( sigmasq );


V = [1:N]';
% Number of chains:
C = floor(sqrt(N));
M = floor(N/C); % number of nodes in a chain
if isRow
    V = [1:N];
else
    Vt = [reshape( V(1:C*M), C, M )]';
    V = [Vt(:);V(C*M+1:end)];
end

CV = {}; % Chains
ConnectorVs = []; % Connector vertices
for ccnt=1:C
    CV{ccnt} = V([(ccnt-1)*M+1:ccnt*M]);
    ConnectorVs(ccnt) = V(round( (ccnt-1)*M+M/2 ) );
end
if C*M<N
    CV{end+1}=V([C*M+1:N]);
    ConnectorVs(end+1) = V( floor(C*M+(N-CM)/2) );
end

E = [];
Lambda = zeros(N,N);
for ccnt=1:numel(CV)
    % On the path from the connector node
    chain = CV{ccnt};
    connode = ConnectorVs(ccnt);
    conind = find( chain==connode );

    for vcnt=1:conind-1
        e = [chain(vcnt), chain(vcnt+1)];
        re = [chain(vcnt+1), chain(vcnt)];
        E = [E;e;re];
        
        % Below is addition of the info mat for e(1) given e(2)
        Lambda( e(1), e(1) ) = Lambda( e(1), e(1) ) + 1/( sigmasq*(1-rho^2) );
        Lambda( e(1), e(2) ) = Lambda( e(1), e(2) ) - rho/( sigma*sigma*(1-rho^2) );
        Lambda( e(2), e(1) ) = Lambda( e(2), e(1) ) - rho/( sigma*sigma*(1-rho^2) );
        Lambda( e(2), e(2) ) = Lambda( e(2), e(2) ) + rho^2/( sigmasq*(1-rho^2) );     
    end

    for vcnt=conind:length( chain )-1
        e = [chain(vcnt+1), chain(vcnt)];
        re = [chain(vcnt), chain(vcnt+1) ];
        E = [E;re;e];
        
        % Below is addition of the info mat for e(1) given e(2)
        Lambda( e(1), e(1) ) = Lambda( e(1), e(1) ) + 1/( sigmasq*(1-rho^2) );
        Lambda( e(1), e(2) ) = Lambda( e(1), e(2) ) - rho/( sigma*sigma*(1-rho^2) );
        Lambda( e(2), e(1) ) = Lambda( e(2), e(1) ) - rho/( sigma*sigma*(1-rho^2) );
        Lambda( e(2), e(2) ) = Lambda( e(2), e(2) ) + rho^2/( sigmasq*(1-rho^2) );
    end
end

% Add now the connector chain
for vcnt = 1:length( ConnectorVs )-1
    e = [ConnectorVs(vcnt),ConnectorVs(vcnt+1)];
    re = [ConnectorVs(vcnt+1), ConnectorVs(vcnt)];
    E = [E;e;re];

    % Below is addition of the info mat for e(1) given e(2)
    Lambda( e(1), e(1) ) = Lambda( e(1), e(1) ) + 1/( sigmasq*(1-rho^2) );
    Lambda( e(1), e(2) ) = Lambda( e(1), e(2) ) - rho/( sigma*sigma*(1-rho^2) );
    Lambda( e(2), e(1) ) = Lambda( e(2), e(1) ) - rho/( sigma*sigma*(1-rho^2) );
    Lambda( e(2), e(2) ) = Lambda( e(2), e(2) ) + rho^2/( sigmasq*(1-rho^2) );
end
lastnode = ConnectorVs(end);
Lambda(lastnode,lastnode)= Lambda(lastnode,lastnode) + 1/sigmasq;




Cx = Lambda^(-1);

if nargout>=3
    varargout{1} = Cx;
end

