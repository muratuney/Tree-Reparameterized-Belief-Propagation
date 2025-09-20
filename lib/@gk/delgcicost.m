function [g, varargout ] = delgcicost( finp, w )

% Murat Uney August 2021
[g, fw ] = delgciObjective( finp, w )

g = -1*g;
if nargout>=2
    varargout{1} = fw;
end
end