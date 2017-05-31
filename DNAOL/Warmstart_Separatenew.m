function [ Wc,Ac,Zc,Fc,lambda ] = Warmstart_Separatenew( Wc,Xc,Xcin,Xrest,Ac,Zc,alpha,rho,lambda,Maxiter )
%UNTITLED3 此处显示有关此函数的摘要
% warm start: alternating
%min
%|Xc-WcFc|_F^2+alpha*|AXcrest|_F^2+rho|S(Zc,lambda)-Fc|_F^2+rho|Zc-AcXc|_F^2
%s.t. |wc|_2<=1
%   此处显示详细说明

 Ic=eye(size(Wc,2));
iter=1;
func=zeros(1,Maxiter);
while iter<=Maxiter    
    % update F
    Fc=(Wc'*Wc+rho*Ic)\(rho*Sel(Zc,lambda)+Wc'*Xc);    
    % update Z     
    [ Zc ] =UpdateZ( lambda,1, Fc,Ac*Xc );     
    % update lambda
    [ lambda ] = Updatelambda_inside( Zc,Fc,lambda );    
    % update A
    Ac=rho*Zc*Xcin;     
    % update W
   [ Wc ] = ClassfierUpdate_Sep( Xc,Fc,Wc,10 ); 
%       func(iter)=norm(Xc-Wc*Fc,'fro')^2+alpha*norm(Ac*Xrest,'fro')^2+rho*norm(Sel(Zc,lambda)-Fc,'fro')^2+rho*norm(Zc-Ac*Xc,'fro')^2;
    iter=iter+1;
end% end while



end

