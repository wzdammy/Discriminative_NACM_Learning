function [ Zout ] =UpdateZ( lambda,alpha, L1,L2 )
% Sub-Optimization with respect to Z
% Solve: min_Z ||S(Z,lambda)-L1||_2^2+alpha||Z-L2||_2^2
%
a=(lambda*L1+alpha*L2)./(lambda^2+alpha);
b=1/(1+alpha/lambda^2);
Z1=max(1,a+b);
Z2=max(-1,min(L2,1));
Z3=min(-1,a-b);
fun1 = func(Z1,lambda,alpha, L1,L2);
fun2 = func(Z2,lambda,alpha, L1,L2);
fun3 = func(Z3,lambda,alpha, L1,L2);

Zout=Z2;
minimal12=find(fun1<fun2);

Zout(minimal12)=Z1(minimal12);
fun4 = func(Zout,lambda,alpha, L1,L2);

minimal23=find(fun3<fun4);

Zout(minimal23)=Z3(minimal23);

end
function out = func(Z,lambda,alpha, L1,L2)
out= (Sel(Z,lambda)-L1).^2+alpha*(Z-L2).^2;

end

