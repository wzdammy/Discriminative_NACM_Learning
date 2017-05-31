function [ A,Fall,Zall,lambda,Lambda1,Lambda2] = NonSep_DNAOL_ADMM( WWtall,WtYall,Xallin,Xall,A,Zall,Lambda1,Lambda2,lambda,Maxiter,rho )

% main function of NonSep-DNAOL
% loss function 0.5*||Yc-WcFc||_F^2+tau/2||Ac||_F^2+beta/2||wc||_F^2 s.t.
% Fc=S(AcXc,lambda)
% learning separable analsis operator for class c:Ac, Fc, with Wc,lambda
% lambda and Zc updated simultaneously
% fixed
%% init

% [~,N]=size(Xc);
% [p,~]=size(Ac);
% F=zeros(p,N);
% Lambda1=zeros(p,N);
% Lambda2=Lambda1;
%  Zc=Ac*Xc;
% Zc=randn(p,N);
%main
iter=1;
% Inc=eye(size(group{1},2));
% In=eye(n);
PR1=zeros(1,Maxiter);
PR2=PR1;
PR=PR1;
DR=PR1;
DR1=PR1;
DR2=PR2;
while iter<=Maxiter
    %update Fc
    Fcosparse=Sel(Zall,lambda);
    Fall=WWtall\(rho*Fcosparse+Lambda2+WtYall);  
  
    %check feasible
    if iter==1
        AX=A*Xall;
    else
        AX=AXnew;
    end
    Lambda1tilde=AX-Lambda1/rho;
    Lambda2tilde=Fall-Lambda2/rho;
   [ Zallnew ] =UpdateZ( lambda,1, Lambda2tilde,Lambda1tilde );  
    [ lambda ] = Updatelambda_inside( Zallnew,Lambda2tilde,lambda );
        % update Ac
    A=(rho*Zallnew+Lambda1)*Xallin;
   
    % compute primal residual
    AXnew=A*Xall;
    Presi1=Zallnew-AXnew;
    Selznew=Sel(Zallnew,lambda);
    Presi2=Selznew-Fall;
    % compute dual residual
%     Dresi1=rho*(Selznew-Fcosparse);
%     Dresi2=rho*(AX-AXnew);
    PR1(iter)=norm(Presi1,inf);
    PR2(iter)=norm(Presi2,inf);
%     DR2(iter)=norm(Dresi2,inf);
%     DR1(iter)=norm(Dresi1,inf);
    PR(iter)=max(PR1(iter),PR2(iter));
%     DR(iter)=max(DR1(iter),DR2(iter));
       % Update lagrangian
    Lambda1=rho*Presi1+Lambda1;
    Lambda2=rho*Presi2+Lambda2;
%     if iter>1
%       if PR(iter)<1e-2&&DR(iter)<1e-2
%         break;
%       end
%       if abs(PR(iter)-PR(iter-1))<1e-3&&abs(DR(iter)-DR(iter-1))<1e-3
%         break;
%       end 
%         
%     end
    if iter>1
      if PR(iter)<1e-3
        break;
      end
      if abs(PR(iter)-PR(iter-1))<1e-3
        break;
      end 
        
    end
    iter=iter+1;
    Zall=Zallnew;

 
end% end while main

    Zall=Zallnew;

end

