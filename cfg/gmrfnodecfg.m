classdef gmrfnodecfg
    properties
        id = 1;
        parents = [];
        children = [];
        neighbours = [];
        iternum = 0;
        nodepot =  gpot;         % initial node state (density)      
        edgepotentials = []; % edge potentials (over the edges towards the children)
        Lambdan        % This is inverse covariance of the noise distribution 
        Cn             % This is the noise covariance
        obsmat         % this is the measurement transform H in y = Hx+n
        y              % measurement at node
        ydim           % dimensionality
    end
    methods
        function o = gmrfnodecfg(varargin)
            if nargin>=1
                if isa( varargin{1}, 'gmrfnodecfg' )
                    % initialize with this config
                    o = varargin{1};
                elseif isempty(  varargin{1} )

                    o = o([]);
                else
                    error('Unknown variable input');
                end
            end
        end
    end
end
