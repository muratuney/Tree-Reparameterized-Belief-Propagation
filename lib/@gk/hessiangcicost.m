function [H, varargout] = hessiangcicost( finp, w )
% Murat Uney August 2021

[H, fw] = hessiangciObjective( finp, w );
H = -1*H;

if nargout>=2
    varargout{1} = fw;
end
end