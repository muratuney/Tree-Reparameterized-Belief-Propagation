function [Lambda, E, varargout] = attractrepulsegmrf( N, varargin )
% attractrepulsegmrf returns a loopy grid Gaussian Markov random field (GMRF)
% of N nodes. An undirected MRF is produced. The default cannonical correlation
% coeffcient rho = 0.5 and the standard deviation of variables are 1. 
% For rho>0, the model is attractive whereas for rho<0, the model is repulsive.
%
% [L, E]=attractrepulsegmrf(N) returns the information matrix (or, the
% inverse covariance) L along with the edge set E. Note that the off-diagonal
% non-zero entries of L are encoded in E. E is an undirected edge set, i.e. if
% (i,j) is in E, then so is (j,i).
%
% [L, E]=attractrepulsegmrf(N, 'rho', r ) 
% [L, E]=attractrepulsegmrf(N,  r ) uses the canonical correlation coefficient r
% for each edge E. 
%
% [L, E]=attractrepulsegmrf(N, 'variance', v ) 
% [L, E]=attractrepulsegmrf(N,  v ) uses the variance value v for each variable.

% Murat Uney 09.2025


sigmasq = 1;
rho = 0.5;

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
        
    else
        error('Unknown inputs. Type help gausstreeovergrid');
    end
end
sigma = sqrt( sigmasq );


V = [1:N]';
% Number of chains:
C = floor(sqrt(N));
M = floor(N/C); % number of nodes in a chain
V = [1:N];

CV = {}; % Chains
for ccnt=1:C
    CV{ccnt} = V([(ccnt-1)*M+1:ccnt*M]);
end
if C*M<N
    CV{end+1}=V([C*M+1:N]);
end

E = [];
Lambda = zeros(N,N);
for ccnt=1:numel(CV)
    chain = CV{ccnt};
    for vcnt=1:length( chain )-1
        e = [chain(vcnt), chain(vcnt+1)];
        re = [chain(vcnt+1), chain(vcnt)];
        E = [E;e;re];
        
           
    end
end

% Add now the connector chains
for ccnt=1:numel(CV)-1
    % On the path from the connector node
    chain1 = CV{ccnt};
    chain2 = CV{ccnt+1};
    for vcnt=1:min( length( chain1 ), length( chain2 )  )
        e =  [chain1(vcnt), chain2(vcnt)];
        re = [chain2(vcnt), chain1(vcnt)];
        E = [E;e;re];
    end
end

% Lambda = I - R

% 1) Lambda = I
for cnt=1:N
    Lambda(cnt,cnt) = sigmasq;
end

% 2) Lambda = Lambda - R
numE = size(E,1);
for cnt=1:numE
    Lambda( E(cnt,1), E(cnt,2) ) = - rho*sqrt( Lambda( E(cnt,1), E(cnt,1) )*Lambda( E(cnt,2), E(cnt,2) ) );
end



Cx = Lambda^(-1);

if nargout>=3
    varargout{1} = Cx;
end

