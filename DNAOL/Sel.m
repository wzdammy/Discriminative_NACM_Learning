function [ xout ] = Sel( x,lambda )
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% zeromat=zeros(size(x));
% xout=zeromat;
xout=sign(x).*(lambda*max(abs(x)-1,0));

end

