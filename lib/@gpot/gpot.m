classdef gpot % Gaussian potential class gpot
    properties
        numvariables % number of variables
        ids          % identification number of the variables
        dims         % dimensionalities of varibles
        nu           % The vector parameter of the information form of the Gaussian potential
        Lambda       % The "inverse covariance" parameter of the Gaussian potential
        % labels       % labels of variables (reserved)
      end
    properties (setAccess =protected)
        
    end
    properties ( setAccess = private )
    
    end
    
    methods
        function e = gpot(varargin)
            e.numvariables = 0;
            e.dims = []';
            e.nu = []';
            e.Lambda = [];
            
            if nargin==1
                if isempty( varargin{1} )
                    e = e([]);
                elseif isnumeric( varargin{1} )
                    e = repmat( e, varargin{1} );
                elseif isa( varargin{1}, 'gpot')
                    e = varargin{1};
                else
                    error('Unknown input!');
                end
            end
            
        end
        function d = getdims(these)
            d = [];
            for cnt=1:numel(these)
                dims = these(cnt).dims;
                d = [d; dims(:)];
            end
        end
        function those = gpot2pdf(these)
            those = gk([]);
            for i=1:length(these)
                Lambda = these(i).Lambda;
                nu = these(i).nu;

                R = chol(Lambda);
                Rinv = R^(-1);
                C = Rinv*Rinv';
                m = C*nu; 
                those(i) = cpdf(gk( C, m ));
            end
        end
    end
    methods (Static)
        [ np, ep, E ] = findinfogmrfpotentials( nu, Lambda, E, varargin  ) % This function is implemented in a separate file
        [ np, ep, E ] = findgmrfpotentials( m, C, E, varargin  ) % This function is implemented in a separate file
        [npe, epe] = subgraphpotentials( np, ep, E ) % This function is implemented in a separate file
        [ Lambda, nu, varargout] = pot2infomat( np, ep ) % This function is implemented in a separate file
    end
end
