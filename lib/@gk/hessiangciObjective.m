function [H, varargout] = hessiangciObjective( finp, w )

% Murat Uney August 2021
N = length( w );

fw = emd( finp, w );

H = zeros(  N, N );

for cnti=1:N
    
    for cntj = 1:N
        
        fi = finp( cnti );
        fj = finp( cntj );
        
        H_fw_fi = crossentropy( fw, fi );
        H_fw_fj = crossentropy( fw, fj );
        
        F_fw_fi_fj = capitalF( fw, fi, fj );
        
        G_fw_fj = capitalG( fw, fj );
        G_fw_fi = capitalG( fw, fi );
        
        H_fw = entropy( fw );
        
        sumw = sum( w );
        
        % Below is the partial derivative of Eq.(45)
%         hval = (2- sumw )*H_fw_fj*H_fw_fi...
%             -( (1- sumw)*H_fw + (1-sumw) + 1 )*F_fw_fi_fj...
%             -( 1- (1-sumw)  )*H_fw_fj*H_fw...
%             -H_fw_fi*H_fw...
%             + 2*(1-sumw)*H_fw_fj*H_fw_fi*H_fw...
%             + (  1- (1-sumw) - (1-sumw)*H_fw_fi  )*G_fw_fj...
%             + G_fw_fi;
        
        % % Below is the partial derivative of (45) before collecting similar terms
        % % together
        hval2 = H_fw_fj*H_fw_fi - F_fw_fi_fj ...
            - H_fw_fj*H_fw + G_fw_fj...
            - H_fw_fi*H_fw + (1-sumw)*H_fw_fj*H_fw_fi*H_fw...
            - (1-sumw)*F_fw_fi_fj*H_fw...
            + (1-sumw)*H_fw_fi*H_fw_fj*H_fw...
            -(1-sumw)*H_fw_fi*G_fw_fj...
            + G_fw_fi ...
            - (1-sumw)*( G_fw_fj - H_fw*H_fw_fj - H_fw_fi*H_fw_fj + F_fw_fi_fj  );
        
        H( cntj, cnti ) = hval2;
        
    end
end
% Comment out below to ensure symmetric H
%H = ( H + H' )/2;

if nargout>=2
    varargout{1} = fw;
end
end
function g = capitalG( f1, f2 )
N = 10000;

samples = f1.gensamples( N );
ef1 = f1.evaluatelog( samples );
ef2 = f2.evaluatelog( samples );
g = (1/N)*sum( ef1.*ef2 );

end
function g = capitalF( f1, f2, f3 )
N = 10000;

samples = f1.gensamples( N );
ef2 = f2.evaluatelog( samples );
ef3 = f3.evaluatelog( samples );
g = (1/N)*sum( ef2.*ef3 );

end


