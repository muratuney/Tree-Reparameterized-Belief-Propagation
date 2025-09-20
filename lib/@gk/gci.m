function [gf, varargout ] = gci( gi, varargin )
% GCI is the Generalised Covariance Intersection for Gaussian
% distributions. The input is an array of Gaussian densities and the output
% is a fused density that is the barycentre with respect to the 
% Kullback-Leibler divergence. 
%
% gf = gci([g1, g2, ..., gN ] ) returns the fused density as a Gaussian Kernel
% @gk object gf when input N Gaussian densities g1, g2, ..., gN that are also
% @gk objects.
%
% [gf, w] = gci([g1, g2, g3, g4 ] ) returns also the optimal weights in an 
% arrray w that characterises the KLD barycentre.
%
% gf = gci([g1, g2, g3, g4 ], '' )
%

% Murat Uney May 2021


isInitialPointIn = 0; % Is the initial point passed in the input arguments?
isDiagnostics = 0;   % Is diagnostics message passed
diagMessage = 'off'; % (default is diagnosics off)
isDisplay = 0; % Is display message passed
dispMessage = 'off'; % (default is display off)

isOptions = 0; % Is further options passed
optStr = ''; % (default is empty)


isInteriorPoint = 1; % interior-point algorithm by default
tolerance = 1.0e-3; % Iteration stopping tolerance
itermax = 2000;

nvarargin = length(varargin);
argnum = 1;
while argnum<=nvarargin
    if isa( varargin{argnum} , 'char')
        switch lower(varargin{argnum})
           
            case {'tolerance'}
                if argnum + 1 <= nvarargin
                    tolerance = varargin{argnum+1};
                    argnum = argnum + 1;
                end
            case {'display'}
                if argnum + 1 <= nvarargin
                    dispMessage = varargin{argnum+1};
                    argnum = argnum + 1;
                end
            case {'algorithm'}
                if argnum + 1 <= nvarargin
                    if strcmp( lower(varargin{argnum+1}), 'interior-point')
                        isInteriorPoint = 1;
                    elseif strcmp( lower(varargin{argnum+1}), 'trust-region-reflective' )
                        isInteriorPoint = 0;
                    else
                        error('Unknown algorithm %s selected',lower(varargin{argnum+1}) );
                    end
                    argnum = argnum + 1;
                end
                
           case {'diagnostic','diagnostics'}
                if argnum + 1 <= nvarargin
                    diagMessage = varargin{argnum+1};
                    argnum = argnum + 1;
                end
            case {'options'}
                if argnum + 1 <= nvarargin
                    optStr = varargin{argnum+1};
                    argnum = argnum + 1;
                end
            case {'initial-point'}
                if argnum + 1 <= nvarargin
                    Theta0 = varargin{argnum+1};
                    if ~isempty(Theta0)
                        isInitialPointIn = 1;
                    end
                    argnum = argnum + 1;
                end
            otherwise
                error('Wrong input string');
        end
    end
    argnum = argnum + 1;   
end

gi = cpdf( gi(:) ); % Ensure that the Gaussian kernels 
N = length( gi );

if isInitialPointIn==0
    % Initial point
    w0 = ones( N, 1)/N;
    lambda0 = avkld2emd( gi, w0 );
    Theta0 = [w0; lambda0 ];
end

gradfun = @gciLagrangianGradient; % log-likelihood gradient
funHessian = @gciLagrangianHessian; % log-likelihood Hessian
gcicostfun = @gciLagrangian;

Thetaprev = Theta0 + tolerance*ones( N+1, 1 );
Theta = Theta0;

itercnt= 0;
Thetas = Theta;

fevals(1) = gciLagrangian( gi, Theta  );
residual = inf;
while( itercnt < itermax && residual > tolerance )
%while( ( norm( Theta - Thetaprev  ) > tolerance ) && itercnt < itermax && residual > tolerance )
    g =   gciLagrangianGradient( gi, Theta   );
    H =   gciLagrangianHessian( gi, Theta );
    
    Thetaprev = Theta;
    Theta = Theta - 0.05*inv( H )*g;
    %Theta = Theta + 0.001*g;
    w = Theta(1:N);
    lambda = Theta(end);
    
    if ~isempty( find(w<0) )
        disp(sprintf('Step %d, negative weight, projecting onto U', itercnt ) );
        w( w< 0 ) = eps;
        w = w/sum(w);
        lambda = avkld2emd( gi, w );
    end
    if ~isempty( find( w> 1) )
        disp(sprintf('Step %d, weight>1, projecting onto U', itercnt ) );
        w( w>1 ) = 1-eps;
        w = w/sum(w);
        lambda = avkld2emd( gi, w );
    end
    
    Theta(1:N) = w;
    Theta(end) = lambda;
    
    itercnt = itercnt+1;  
    Thetas(:,itercnt+1) = Theta;
    fevals(itercnt+1) = gciLagrangian( gi, Theta  );
    
    kldvals = kldfromemd(gi,w);
    residual = sum( abs( kldvals - lambda) );
    
end

if residual > tolerance 
   [mval, ind] = max( fevals );
   Theta = Thetas(:, ind );
   %Thetas = Thetas(:, 1:ind );
   %fevals = fevals(1:ind);
end

w = Theta(1:N);
gf = emd( gi, w );
lambda = Theta(end);

if nargout>=2
    varargout{1} = w;
end
if nargout>=3
    varargout{2} = lambda;
end
if nargout>=4
    varargout{3} = Thetas;
end
if nargout>=5
    varargout{4} = fevals;
end




