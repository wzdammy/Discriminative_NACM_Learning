function [ totalnum,rate ] = recog_rate( result,template )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
a=result-template;
totalnum=length(find(a==0));
rate=totalnum./length(template);

end

