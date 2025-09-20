function varargout = init( this, varargin )


if nargin>=2
    cfg = varargin{1};
    if ~isa( cfg, 'gmrfnodecfg' )
        error('Unknown configuration object!');
    end
    this.cfg = cfg;
else
    if ~isa( this.cfg, 'gmrfnodecfg' )
        error('Unknown configuration object provided!');
    end
    cfg = this.cfg;
end

% Here perform the assignments
this.id = cfg.id;
this.parents = cfg.parents;
this.children = cfg.children;
this.nodepot = cfg.nodepot;

Cmat = inv( this.nodepot.Lambda );
mvect = Cmat*this.nodepot.nu;
this.infostate = this.nodepot;
this.initstate = cpdf(gk( Cmat, mvect ));

this.dim = this.nodepot.getdims;
if ~isempty( cfg.Lambdan )
    this.Lambdan = cfg.Lambdan;
    this.Cn = inv( this.Lambdan );
elseif ~isempty(cfg.Cn)
    this.Cn = cfg.Cn;
    this.Lambdan = inv( this.Cn );
end
this.obsmat = cfg.obsmat;
this.y = cfg.y;
this.ydim = numel(this.y);
this.iternum = cfg.iternum;

this.edgepotentials = cfg.edgepotentials;

this.inboxlog = zeros(1, length(this.parents) );

this.inbox = gpot( [] );
this.previnbox = gpot( [] );

this.previnboxlog = this.inboxlog;

if nargout == 0
    if ~isempty( inputname(1) )
        assignin('caller',inputname(1),this);
    else
        error('Could not overwrite the instance; make sure that the argument is not in an array!');
    end
else
    varargout{1} = this;
end
