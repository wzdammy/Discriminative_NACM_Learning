function [ Wc ] = ClassfierUpdate_Sep( Xc,Fc,Wc,Maxiter )
% classifier update
% min ||Xc-WcFc||_2^2 s.t. |w_i|<=1
%  
rho=1.5;
FFt=Fc*Fc';
FInv=FFt+rho*eye(size(FFt));
iter=1;
XFt=Xc*Fc';
Lambda=zeros(size(Wc));
Zc=Lambda;
while iter<Maxiter 
    Wc=(XFt+rho*Zc+Lambda)/FInv;
    Zc=matrixNormalize(Wc-Lambda/rho,1,1);
    R=Zc-Wc;
    Lambda=Lambda+rho*(R);
    if MatInf(R)<1e-3
        break;
    else
        iter=iter+1;
    end    
end
end
function [ infnorm ] = MatInf( X )

% compute max(max(|X|))


infnorm=max(max(abs(X)));
end

%--------------------------------------------------------------------------

function [Yn] = matrixNormalize(Y,r,mode)
colnorm=sqrt(sum(Y.^2,1));
Colmat=repmat(colnorm,[size(Y,1),1]);
if mode==1
    [~,label]=find(colnorm>r);
    Yn=Y;
    Yn(:,label)=Y(:,label)./Colmat(:,label);
else
  Yn=Y./Colmat;
end
end
