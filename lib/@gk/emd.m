function gf = emd( gi, varargin )
%
%
% 


if length(varargin) == 2 
    if isa( varargin{1}, 'gk' )
        if isnumeric( varargin{2} )
            % Backward compatibility
            gf = emdold( gi, varargin{1}, varargin{2} );
            return;
        else
            error('Function call with three arguments should a weight coefficient as the third input');
        end      
    else
        error('Function call with three arguments should have a @gk object as the second input.');
    end
elseif length(varargin) == 1
    if isnumeric( varargin{1} )
        gf = cpdf( prod( gi.^varargin{1} ) );
    else
        error('Function call with two artguments should the EMD weight array as the second input');
    end
elseif length(varargin) == 0
      N = length( gi );
       gf = cpdf( prod( gi(:).^(ones(N,1)/N) ) );   
else
    error('Incorrect inputs!');
end



end
function d = emdold( a, b, w)
% d = emd( a, b, w)
% returns the exponential mixture density of a and b with weight w, i.e.
% d = (1/Z)*a^(1-w)*b^(w)
% where d is a new Gaussian. Z can be found by zomega( a, b, w )
%
% See also ZOMEGA, CISG1

d = cpdf( a^(1-w)*b^(w) );
end