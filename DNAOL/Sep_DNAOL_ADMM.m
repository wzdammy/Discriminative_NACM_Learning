function [ Ac,Fc,Zc,lambda,Lambda1,Lambda2] = DNAOL_Separable3_ADMM( WWtc,WTyc,Xcin,Xc,Ac,Zc,Lambda1,Lambda2,lambda,Maxiter,rho )
% min function of Sep-DNAOL
% loss function 0.5*||Yc-WcFc||_F^2+alpha/2||AcXcrest||_F^2+tau/2||Ac||_F^2 s.t.
% Fc=S(AcXc,lambda), ||wc||<=1
% learning separable analsis operator for class c:Ac, Fc, with Wc,lambda
% lambda and Zc updated simultaneously
% fixed
iter=1;
PR1=zeros(1,Maxiter);
PR2=PR1;
PR=PR1;
DR=PR1;
DR1=PR1;
DR2=PR2;
while iter<=Maxiter
    %update Fc
    Fcosparse=Sel(Zc,lambda);
    Fc=WWtc\(rho*Fcosparse+Lambda2+WTyc);   
    if iter==1
        AX=Ac*Xc;
    else
        AX=AXnew;
    end
    Lambda1tilde=AX-Lambda1/rho;
    Lambda2tilde=Fc-Lambda2/rho;   
    [ Znew ] =UpdateZ( lambda,1, Lambda2tilde,Lambda1tilde );  
    [ lambda ] = Updatelambda_inside( Znew,Lambda2tilde,lambda );   % 
        % update Ac        
    Ac=(rho*Znew+Lambda1)*Xcin;   
    % compute primal residual
    AXnew=Ac*Xc;
    Presi1=Znew-AXnew;
    Selznew=Sel(Znew,lambda);
    Presi2=Selznew-Fc;
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
    Zc=Znew;

 
end% end while main



end

