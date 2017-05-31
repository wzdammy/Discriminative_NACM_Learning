function [ lambda ] = Updatelambda_inside( Z,F,lambda )
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
%     update lambda classwise
% min ||Sel(Z,lambda)-F||_F^2
%     Z=A*X;
    Posi1=find(Z>1);
    ap=Z(Posi1)-1;
    bp=F(Posi1);
    Neg1=find(Z<-1);
    aq=Z(Neg1)+1;
    bq=F(Neg1);
    up=((sum(ap.*bp))+(sum(aq.*bq)));
    down=(((sum(ap.^2)))+((sum(aq.^2))));
if down~=0
    lambda=up/down;
end
% if lambda<-1
%     lambda=-1;
% end
% if lambda>0
%     lambda=1;
% end   

end

