function varargout = init( this, varargin )


global DEBUG_GMRF DEBUG_VERBOSE

if nargin>=2
    cfg = varargin{1};
    if ~isa( cfg, 'gmrfcfg' )
        error('Unknown configuration object!');
    end
    this.cfg = cfg;
else
    if ~isa( this.cfg, 'gmrfcfg' )
        error('Unknown configuration object provided!');
    end
    cfg = this.cfg;
end



if exist('DEBUG_VERBOSE')
    if DEBUG_VERBOSE
        disp('inside @gmrf.init');
    end
end

% Here perform the assignments
this.V = unique( cfg.V(:), 'legacy' );

% Verify that the edge set satisfies undirectedness,
% i.e., for every (i,j) in E, (j,i) is in E:
[iu, ei ] = isundirected( cfg.E );
if iu == 0
    % Edge set is not undirected
    warning(sprintf('Edge set is not undirected, i.e., the following (i,j), do not have their reverse (j,i) in E:%d,', ei));
end

this.E = cfg.E;
this.edgepots = cfg.edgepots;
this.itermax = cfg.itermax;
this.iternum = 0;
this.treeweight = cfg.treeweight;
this.numberofmessages = 0;
if isempty( cfg.mschedule)
    this.mschedule = scheduler('pattern', {this.E} );
else
    this.mschedule = scheduler('pattern', cfg.mschedule );  
end


this.nodes = gmrfnode([]);
this.N = size( this.V, 1 );
this.M = size( this.E, 1 );
% Create the node objects
for i=1:this.N 
    mynodecfg = this.cfg.nodes(i); % Get node configuration object
    mynodecfg.id = this.V(i);
    
    % Find the children
    [chi_, einds ] = chi( this.E, mynodecfg.id );
    mynodecfg.children = chi_;
    mynodecfg.edgepotentials = this.cfg.edgepots(einds );
    
    pa_ = pa( this.E, mynodecfg.id );
    mynodecfg.parents = pa( this.E, mynodecfg.id );

    nei_ = nei( this.E, mynodecfg.id );
    mynodecfg.neighbours = nei_;
    
    mynodecfg.iternum = this.iternum;
    
    nodeobj = gmrfnode( mynodecfg );
    this.nodes(i) = nodeobj;
    this.dims(i) = nodeobj.getstatedim;
end

if nargout == 0
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end
else
    varargout{1} = this;
end
