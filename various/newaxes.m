function h = newaxes(varargin)
% NEWAXES opens a new axes, configures it and returns a handle
%
% h = newaxes returns a handle h to an axes with 15 point fonts
%
% h = newaxes( fontSize ) returns an axes with fontSize fonts
%
% h = newaxes( [n m l], fontsize ) returns the lth axes in an n-by-m
% tiling of the figure
%
% See also newfigure

fontSize = 15;
n = 1;
m = 1;
l = 1;
if nargin==1
    if length( varargin{1} ) == 1
        fontSize = varargin{1}(1);
    else
        n = varargin{1}(1);
        m = varargin{1}(2);
        l = varargin{1}(3:end);
    end
end
if nargin==2
    if length( varargin{2} ) == 1
        fontSize = varargin{2}(1);
    end

    n = varargin{1}(1);
    m = varargin{1}(2);
    l = varargin{1}(3:end);
end
if nargin==3
    n = varargin{1}(1);
    m = varargin{2}(1);
    l = varargin{3}(1:end);
end

    
h = subplot(n,m,l);
set( h, 'FontSize', fontSize );
hold on
grid on


