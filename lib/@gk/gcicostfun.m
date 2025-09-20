function [mJ, Del, H] = gcicostfun( fi, w )

mJ = -1*gciObjective( fi, w );

if nargout>=2
    Del = -1*delgciObjective( fi, w  );
end


if nargout>=3
    H = -1*hessiangciObjective( fi, w  );
end