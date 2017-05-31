% Demo of Discriminative Nonlinear Analysis Cosparse Opeartor Learning on AR database
% Test: NonSep_DNAOL function
% copyright: Zaidao Wen
% please cite the corresponding paper: ¡°Discriminative nonlinear analysis operator learning: When cosparse model meets image classification,¡± IEEE Trans. Image Process., vol. 26, no. 7, pp. 3449¨C3462, July 2017.
clc
clear
close all
load ARdataset.mat

%% parameters init
n=size(Xtrain{1},1);
In=eye(n);
Classnum=size(Xtrain,2);
trainnum=20;
p=Classnum*trainnum;% feature dimension
Maxiter=30;
tau=3e-4;
lambda=rand;
%% load data and init.
all=1:trainnum*Classnum;
Xall=cell2mat(Xtrain);
XXT=Xall*Xall';
Xtestall=cell2mat(Xtest);
W=randn(Classnum,p);
%  W=matrixNormalize(W,1,1);
A=2*randn(p,n);
Zall=randn(p,size(Xall,2));
Fall=randn(p,size(Xall,2));
Lambda1=zeros(size(Zall));
Lambda2=zeros(size(Zall));
Yall=Creat2Dlabeltemplate( Xtrain,Classnum);
rho=1;
[ templatetrain ] = Createtemplatenew( Xtrain);
[ template ] = Createtemplatenew( Xtest);

Xallin=Xall'/(tau*In+rho*XXT);


for i=1:15
[ W ] = ClassfierUpdate_Sep( Yall,Fall,W,20 );
WWtall=W'*W+rho*eye(size(W,2));
WtYall=W'*Yall;
[ A,Fall,Zall,lambda,Lambda1,Lambda2] = NonSep_DNAOL_ADMM( WWtall,WtYall,Xallin,Xall,A,Zall,Lambda1,Lambda2,lambda,Maxiter,rho );
% r(i)=norm(Yall-W*Fall,'fro')^2+norm(W,'fro')*beta+tau*norm(A,'fro');
% % r(i)=sum(loss)
% if i>1
%     if abs(r(i)-r(i-1))/r(i)<1e-2
%         break;
%     end
% end

% % classification
Resitest=W*Sel(A*Xtestall,lambda);

[~,indextest]=max((Resitest),[],1);

[~,rateARtest(i)]=recog_rate(indextest,template)
end


%  matlabpool close