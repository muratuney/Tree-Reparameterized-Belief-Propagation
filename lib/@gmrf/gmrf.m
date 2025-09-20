classdef gmrf 
    properties
        cfg
        nodes % Variable node objects
        nodepots % Node potential objects
        edgepots % Edge potential objects
        V % list of node ids
        dims 
        E % List of edges        
        N % number of nodes
        M % number of edges  
        mschedule %messaging schedule
        uschedule % update schedule
        iternum
        itermax
        numberofmessages 
        treeweight   % Tree weight in tree reparameterised updates
        
    end
    properties (setAccess =protected)
        
    end
    properties ( setAccess = private )
    
    end
    
    methods
        function g = gmrf(varargin)
            
            if nargin == 1
                if isempty( varargin{1} )
                    g = g([]);
                elseif isnumeric( varargin{1} )
                    g = repmat( g, varargin{1} );
                elseif isa( varargin{1}, 'gmrf' )
                    g = varargin{1};
                elseif isa( varargin{1}, 'gmrfcfg' )
                    % initialize with this config
                    mygmrfcfg = varargin{1};
                    g.cfg = mygmrfcfg;
                    g = init( g );
                else
                    error('Unknown variable input');
                end
            elseif nargin == 0
                mygmrfcfg = gmrfcfg;
                g.cfg = mygmrfcfg;
                g = init(g);
            else
                error('Unknown variable input');
            end
            
        end
    end
end
