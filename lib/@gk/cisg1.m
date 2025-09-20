function  [ popt , varargout] = cisg1( pi, pj, varargin )
% function CISG1 is the stochastic gradient implementation of the optimal
% fusion Newton iterations for Gaussians. Equivalently cisg1 is a
% Covariance Intersection (CI) algorithm with the weight selection optimised
% with respect to the equal-KLD criteria.
%
% [ popt ] = cisg1( pi, pj ) outputs the EMD of pi and pj (popt)
% such that popt = 1/Z * pi^(1-w*)pj^w*
% and D(popt||pi) = D(popt||pj)
%
% [ popt, ws ] = cisg1( pi, pj ) also returns the optimal weight parameter
% w* along with the optimal density popt.
%
% [ popt, ws, r ] = cisg1( pi, pj ) returns the residual error
% r = abs(  D(popt||pi) - D(popt||pj) ) for all iterations
%
%  [ popt, ws ] = cisg1( pi, pj, eps ) input eps as the precision of the
% omega when optimising. Technically, this is the termination threshold for
% the Newton iterations which is by default set to 10e-2
%
% [ popt, ws ] = cisg1( pi, pj, eps, n ) input n as the number of samples used
% in the Monte Carlo estimates involved in the Newton iterations. This value
% controls the precision of the omega when optimising.
%
% For more about the Newton CI, Uney, Houssineau, Delande, Julier, Clark,
% "Fusion of finite set distributions: Pointwise consistency and global
% cardinality", IEEE TAES 2018, particularly Algorithm 1.
%
% For more about CI, see Julier and Uhlmann "A non-divergent estimation
% algorithm in the presence of unknown correlations," Americal Control
% Conference 1997.
%
% For more about the weight selection criteria D(popt||pi) = D(popt||pj),
% see Hurley, ""
%
% See also Julier, Uhlmann, "General decentralized data fusion with covariance
% intersection," in Handbook of Multisensor Data Fusion, Chapter 14. CRC Press 2009.
%
% See also EMDCARD, EMDBERN, EMDPOIS, CI, SCI

% Murat Uney 20/08/2018


itermax = 1000; % Maximum number of iterations

Epsilon = 1e-2;
if nargin>=3
    Epsilon = varargin{1}(1);
end


NumSamples = 10000;
if nargin>=4
    NumSamples = varargin{2}(1);
end


Omega_prev = Inf;
Omega_crrnt = 0.5;
if nargin>=5
    if varargin{3}(1)>= 0 && varargin{3}(1)<=1
        Omega_crrnt = varargin{2}(1);
    end
end



k = 1;
Omegas(k) = Omega_crrnt;
rerr = [];
while( abs(Omega_crrnt - Omega_prev) >= Epsilon )
   
   Omega_prev = Omega_crrnt; % Save the current Omega as the previous
  
   % Now, update the current Omega
   
   % 0) Find p_omega
   popt = emd( pi, pj, Omega_crrnt); 
   
   
   % i) Find the Monte Carlo estimate of zdash/z
   
   % Generate samples from p_omega 
   X = popt.gensamples( NumSamples );
   
   % Evaluate log pj(X)/pi(X)
   logratio = pj.evaluatelog(X) - pi.evaluatelog(X);
   
   % MC estimates:
   
   % q1 = sum( logratio )/NumSamples; % Quotient one is zdash/z
   q1 = kld( popt, pi ) - kld( popt, pj );
   
   q2 = sum( logratio.^2 )/NumSamples; % Quotient one is zdashdash/z
   
   rerr(k) = abs( q1 );
   
%    % Update if lnd is not all zero
   if abs( ( q2 - q1^2  ) ) > eps
       Omega_crrnt = Omega_crrnt - q1/( q2 - q1^2  );
       
   end
   % Comment out below to use zdasg = q1*zomega when computing the update
   %  zom =  zomega( pi.m, pi.C, pi.S, pj.m, pj.C, pj.S, Omega_crrnt )
   %  Omega_crrnt = Omega_crrnt - (q1*zom^2)/( (q2*zom^2) - (q1*zom)^2  );
  
   % Update the index
   k = k+1;
  
   % Save the current Omega
   Omegas(k) = Omega_crrnt;
   
   if k-1 >= itermax
       warning(sprintf('Precision of %g not reached, quitting after %d iterations',Epsilon, k-1 ));
       break;
   end
   
end
popt = emd( pi, pj, Omega_crrnt); 


if nargout>=2
    varargout{1} = Omega_crrnt;
end

if nargout>=3
    varargout{2} = rerr;
end

if nargout>=4
    varargout{3} = Omegas;
end



end

function z = zomega( m1, C1, S1, m2, C2, S2, om )
% Find the integral
% z = \int pi(x)^(1-w)*pj(x)^w dx
d = length(m1(:));
fact1 = ( ( (2*pi)^(d/2) ) *det(C1)^((1-om)/2) *det(C2)^(om/2) )^-1;

A = (1-om)*S1 + om*S2;
B = (1-om)*m1'*S1 + om*m2'*S2 + (1-om)*m1'*S1' + om*m2'*S2';
B = 0.5*B';

fact2 = sqrt( (2*pi)^d/det(A) )*...
    exp( -0.5*(  (1-om)*m1'*S1*m1 + om*m2'*S2*m2  ) ...
    + 0.5*B'*inv(A)*B  );

z = fact1*fact2;
end



