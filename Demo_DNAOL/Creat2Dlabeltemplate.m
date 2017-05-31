function [ Template ] = Creat2Dlabeltemplate( Traindata,classnum)
% Crete 2D template
% Template2D=zeros(classnum,classnum*eachnum);
temp=0;
for i=1:classnum
    [~,eachnum]=size(Traindata{i});
    Template(i,temp+1:temp+eachnum)=1;
    temp=temp+eachnum;
  
end


end

