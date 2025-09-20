classdef gmrfnode
    properties
        cfg
        id
        parents
        children
        nodepot       % Node potential
        state         % belief state
        infostate     % Information form of the belief state
        initstate     % pdf version of the node potential
        dim
        Lambdan        % This is inverse covariance of the noise distribution
        Cn             % This is the noise covariance
        obsmat         % this is the measurement transform H in y = Hx+n
        y              % measurement at node
        ydim           % dimensionality
        edgepotentials
        iternum
        inbox
        previnbox
        outbox
        inboxlog
        previnboxlog
        outboxlog
    end
    properties (setAccess =protected)
        
    end
    properties ( setAccess = private )
    
    end
    
    methods
        function o = gmrfnode(varargin)
          if nargin>=1
                if isa( varargin{1}, 'gmrfnodecfg' )
                    % initialize with this config
                    o.init( varargin{1} );
                elseif isempty(  varargin{1} )
                    o = gmrfnode;
                    o = o([]);
                else
                    error('Unknown variable input');
                end
          else
                o.cfg = gmrfnodecfg;
                o.init;
           end 
        end
        
    end
end
