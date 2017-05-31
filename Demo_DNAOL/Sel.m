function [ xout ] = Sel( x,lambda )
% zeromat=zeros(size(x));
% xout=zeromat;
xout=sign(x).*(lambda*max(abs(x)-1,0));

end

