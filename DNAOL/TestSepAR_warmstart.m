clc
clear
close all
load randomfaces4ar.mat
%  matlabpool open local 8
%% parameters init
n=size(featureMat,1);
In=eye(n);
Classnum=size(labelMat,1);
trainnum=20;
subdicnum=20;
p=subdicnum*Classnum;
Maxiter=20;
alpha=0.9e-3;
tau=2e-4;
rho=1;
% lambda=0.5*ones(1,Classnum);
lambda=rand*ones(1,Classnum);
%% load data
temp=0;
all=1:trainnum*Classnum;
for i=1:Classnum
  X{i}=matrixNormalize(featureMat(:,labelMat(i,:)==1),1,0);
 %  X{i}=featureMat(:,labelMat(i,:)==1);
    Position=randperm(size(X{i},2));
    Xtrain{i}=X{i}(:,Position(1:trainnum));
    Xtest{i}=X{i}(:,Position(trainnum+1:size(X{i},2)));
    group{i}=(i-1)*subdicnum+1:subdicnum*i;
    Classgroup{i}=temp+1:temp+size( Xtrain{i},2);
    Classrest{i}=all;
    Classrest{i}(Classgroup{i})=[];
    W{i}=matrixNormalize(randn(n,subdicnum),1,1);
    A{i}=1.5*randn(subdicnum,n);
    temp=temp+size(Xtrain{i},2);
    grouprest{i}=1:subdicnum*Classnum;
    grouprest{i}(group{i})=[];
    Z{i}=randn(subdicnum,trainnum);
    F{i}=Sel(Z{i},lambda(i));
    Lambda1{i}=zeros(subdicnum,trainnum);
    Lambda2{i}=zeros(subdicnum,trainnum);
end
Xall=cell2mat(Xtrain);
Xtestall=cell2mat(Xtest);
 tic

[ templatetrain ] = Createtemplatenew( Xtrain);
[ template ] = Createtemplatenew( Xtest);
% warm start

for c=1:Classnum
Xin{c}=Xtrain{c}'/(alpha*Xall(:,grouprest{c})*Xall(:,grouprest{c})'+rho*Xtrain{c}*Xtrain{c}'+tau*In);
  [  W{c},A{c},Z{c},F{c},lambda(c)] = Warmstart_Separatenew( W{c},Xtrain{c},Xin{c},Xall(:, Classrest{c}),A{c},Z{c},alpha,rho,lambda(c),Maxiter );
% [ W{c},A{c},Z{c},F{c},lambda(c) ] = Warmstart_Separate( W{c},Xtrain{c},Xin{c},Xall(:, Classrest{c}),A{c},Z{c},alpha,tau,rho,lambda(c),100);
%   loss(c)=norm(Xtrain{c}-W{c}*F{c},'fro')^2+alpha*norm(A{c}*Xall(:, Classrest{c}),'fro')^2+rho*norm(Sel(Z{c},lambda(c))-F{c},'fro')^2+rho*norm(Z{c}-A{c}*Xtrain{c},'fro')^2;
end

for i=1:50
   
for c=1:Classnum
%     W{c}=Xtrain{c}*F{c}'/(F{c}*F{c}'+tau*eye(subdicnum));
    [ W{c} ] = ClassfierUpdate_Sep( Xtrain{c},F{c},W{c},10 ); 
    WWt{c}=W{c}'*W{c}+rho*eye(subdicnum);
    WTy{c}=W{c}'*Xtrain{c};
  
     [ A{c},F{c},Z{c},lambda(c),Lambda1{c},Lambda2{c} ] = DNAOL_Separable3_ADMM( WWt{c},WTy{c},Xin{c},Xtrain{c},A{c},Z{c},Lambda1{c},Lambda2{c},lambda(c),Maxiter,rho );
%       [ lambda(c) ] = Updatelambda_independent( A{c},Xtrain{c},F{c},lambda(c) );
   
%      
    loss(c)=norm(Xtrain{c}-W{c}*F{c},'fro')^2/2+alpha/2*norm((A{c}*Xall(:,grouprest{c})),'fro')^2+tau/2*norm(A{c},'fro')^2;

    
end
    


r(i)=sum(loss)
if i>1
    if abs(r(i)-r(i-1))/r(i)<1e-3
        break;
    end
end

% % classification
for c=1:Classnum
%     Grest=1:p;
%     Grest(group{i})=[];
%     Re=Xall-W{c}*Sel(A{c}*Xall,lambda);
    Retest=Xtestall-W{c}*Sel(A{c}*Xtestall,lambda(c));
%     Resi(c,:)=sum(Re.^2,1);
    Resitest(c,:)=sum(Retest.^2,1);
end
% [~,index]=min(Resi,[],1);
[~,indextest]=min(Resitest,[],1);
% [~,rate_Gradpara_EYaleB]=recog_rate(index,templatetrain)
[~,rate_Gradpara_EYaleBtest(i)]=recog_rate(indextest,template)
end

toc
%  matlabpool close