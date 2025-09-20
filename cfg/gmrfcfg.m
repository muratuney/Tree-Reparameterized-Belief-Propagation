classdef gmrfcfg 
    properties
        V = [];
        E = [];
        mschedule
        nodes = [];
        edgepots = [];
        itermax = 10;
        treeweight = 1.0; % 0.25;
    end
    methods
        function g = gmrfcfg
            g.V = [];
            g.E = [];
            g.mschedule =  {} ; % Messaging schedule
            g.nodes = gmrfnodecfg;
            g.edgepots = gpot([]);
        end
    end
end
