function [ xout ] = Sel( x,lambda )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
% zeromat=zeros(size(x));
% xout=zeromat;
xout=sign(x).*(lambda*max(abs(x)-1,0));

end

