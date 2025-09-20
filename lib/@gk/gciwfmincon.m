function [gf, varargout ] = gciwfmincon( gi, varargin )
% GCIWFMINCON is the Generalised Covariance Intersection for Gaussian
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

% Murat Uney August 2021


isInitialPointIn = 0; % Is the initial point passed in the input arguments?
isDiagnostics = 0;   % Is diagnostics message passed
diagMessage = 'off'; % (default is diagnosics off)
isDisplay = 0; % Is display message passed
dispMessage = 'off'; % (default is display off)

isOptions = 0; % Is further options passed
optStr = ''; % (default is empty)


isInteriorPoint = 1; % interior-point algorithm by default
tolerance = 1.0e-5; % Iteration stopping tolerance
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
                    w0 = varargin{argnum+1};
                    if ~isempty(w0)
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
end

gradfun = @delgcicost; % log-likelihood gradient
funHessian = @hessiangcicost; % log-likelihood Hessian
gcifun = @gcicostfun;

w = w0;
fevals(1) = gcicostfun( gi, w0  );

optionsStr = [];
if isInteriorPoint
    optionsStr = ['''fmincon'',''Algorithm'',''interior-point'''];
else
    optionsStr = ['''fmincon'',''Algorithm'',''trust-region-reflective'''];
end

optionsStr = [optionsStr, ','];
optionsStr = [optionsStr, '''Diagnostics'',''',diagMessage,''''];

optionsStr = [optionsStr, ','];
optionsStr = [optionsStr, '''Display'',''',dispMessage,''''];

if isOptions
    optionsStr = [optionsStr, ',', optStr ];
end

eval( ['options = optimoptions(',optionsStr,');']);

% Below the Hessian of the negative log likelihood is passed to the options
% object as the Hessian of the cost to be minimised.
options = optimoptions( options, 'OutputFcn',@outfun, ...
    'GradObj','on','Hessian','user-supplied', 'HessFcn', @(w, lambda)hessiangcicost( gi, w ),...
    'StepTolerance', tolerance );

% Below are the linear transforms for the contraints
LB = [];
UB = [];
A = -eye(N);
b = zeros(N,1);
Aeq = ones(1,N);
beq = 1;

% Set up shared variables with outfun
history.x = [];
history.fval = [];
searchdir = [];

[wstar, gcicostval,EXITFLAG,OUTPUT,lambda, grad, hessian]  = fmincon( @(w)gcicostfun( gi, w ), w0,A,b,Aeq,beq,LB,UB, [], options );
    

w = wstar(1:N);
gf = emd( gi, w );

if nargout>=2
    varargout{1} = w;
end
if nargout>=3
    varargout{2} = lambda.eqlin;
end
if nargout>=4
    varargout{3} = history.x;
end
if nargout>=5
    varargout{4} = -history.fval';
end
function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x, x(:)];
         % Concatenate current search direction with 
         % searchdir.
           searchdir = [searchdir;... 
                        -optimValues.gradient'];
         
         otherwise
     end
end
end




