function h = newfigure(varargin)
% NEWFIGURE opens a new figure, configures it and returns a handle
%
% h = newfigure returns a handle h to a figure with white background
%
% h = newfigure([n,m]) returns a figure of size suitable for an nxm square
% shape axis
%
% See also newaxes


n = 1;
m = 1;
if nargin>=1
    n = varargin{1}(1);
    if length( varargin{1} )> 1
        m = varargin{1}(2);
    end
end
if nargin>=2
    m = varargin{2};
end
    
h = figure;
set( h, 'Color', [1 1 1], 'Units','normalized' );
hs = get(h, 'Position' );

hs(3) = min( hs(3)*m, 1 );
hs(4) = min( hs(4)*n, 1 );

hs(1) = max( min(1 - hs(3), hs(1) ), 0 );
hs(2) = max( min(1 - hs(4), hs(2) ), 0 );

set( h, 'Position', hs );


