function [ template ] = Createtemplatenew( Xtest)
%produce the gndtruth
for i=1:length(Xtest)
    
    if i==1
        template(1:size(Xtest{i},2))=i;
        N=size(Xtest{i},2);
    else
        template(N+1:N+size(Xtest{i},2))=i;
        N=N+size(Xtest{i},2);
    end
end


end

