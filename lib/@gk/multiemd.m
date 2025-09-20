function gf = multiemd( gi, w)
% Multi-EMD
% d = emd( a, b, w)
% returns the exponential mixture density of a and b with weight w, i.e.
% d = (1/Z)*a^(1-w)*b^(w)
% where d is a new Gaussian. Z can be found by zomega( a, b, w )
%
% See also EMD, ZOMEGA, CISG1

d = cpdf( a^(1-w)*b^(w) );