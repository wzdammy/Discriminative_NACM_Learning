% Demo of Discriminative Nonlinear Analysis Cosparse Opeartor Learning on AR database
% Test: Sep_DNAOL function
% copyright: Zaidao Wen
% please cite the corresponding paper: ¡°Discriminative nonlinear analysis operator learning: When cosparse model meets image classification,¡± IEEE Trans. Image Process., vol. 26, no. 7, pp. 3449¨C3462, July 2017.
clc
clear
close all
load ARdataset.mat
%% parameters setting
n=size(Xtrain{1},1);% dimension of the input domain
Classnum=size(Xtrain,2);
trainnum=20;
subdicnum=20;
p=subdicnum*Classnum;
Maxiter=20;
alpha=1e-3;
tau=1e-4;
sigma=4;
rho=1;
%% load data and init.
% data in this demo is already randomly shuffled and separated into training and testing part
temp=0;
all=1:trainnum*Classnum;
Xall=cell2mat(Xtrain);
XXTall=Xall*Xall';
Xtestall=cell2mat(Xtest);
In=eye(n);
for i=1:Classnum
    XXt{i}=Xtrain{i}*Xtrain{i}';
    W{i}=randn(n,subdicnum);
    A{i}=sigma*randn(subdicnum,n);
    Z{i}=randn(subdicnum,trainnum);
    F{i}=randn(subdicnum,trainnum);
    Lambda1{i}=zeros(subdicnum,trainnum);
    Lambda2{i}=zeros(subdicnum,trainnum);
    Xin{i}=Xtrain{i}'/(alpha*XXTall+(rho-alpha)*XXt{i}+tau*In);
end
lambda=rand(1,Classnum);
[ templatetrain ] = Createtemplatenew( Xtrain);
[ template ] = Createtemplatenew( Xtest);
%% main loop 
for i=1:15% using fixed iterations
    for c=1:Classnum% parallel computing is also encouraged
        [ W{c} ] = ClassfierUpdate_Sep( Xtrain{c},F{c},W{c},20 );
        WWt{c}=W{c}'*W{c}+rho*eye(subdicnum);
        WTy{c}=W{c}'*Xtrain{c};
        [ A{c},F{c},Z{c},lambda(c),Lambda1{c},Lambda2{c} ] = Sep_DNAOL_ADMM( WWt{c},WTy{c},Xin{c},Xtrain{c},A{c},Z{c},Lambda1{c},Lambda2{c},lambda(c),Maxiter,rho );
%         loss(c)=norm(Xtrain{c}-W{c}*F{c},'fro')^2+alpha*norm((A{c}*Xall(:,grouprest{c})),'fro')^2+tau*norm(A{c},'fro')^2;   
    end
%     r(i)=sum(loss)
%     if i>1
%         if abs(r(i)-r(i-1))/r(i)<1e-2
%             break;
%         end
%     end
% end
% % classification illustratiion during iteration
% in general this part should be placed out of the main loop if computing
% the time cost

for c=1:Classnum
    Retest=Xtestall-W{c}*Sel(A{c}*Xtestall,lambda(c));
    Resitest(c,:)=sum(Retest.^2,1);
end
[~,indextest]=min(Resitest,[],1);
[~,rate_ARtest(i)]=recog_rate(indextest,template)
end
