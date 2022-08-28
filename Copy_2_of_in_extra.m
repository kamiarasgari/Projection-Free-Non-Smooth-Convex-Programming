clc
clear all
close all
%%
SSigma=20;
TTT=40000;
tau = 1;

for N=[5,10]
    for M = [N,2*N]
        name = ['W=in_','N=',num2str(N),'_M=',num2str(M),'_sigma=',num2str(0),'.mat'];
        load(name)
        %[a,b,c]=svd(rand(N,M));
        %W =0.5*tau*(a*b*c')/sum(sum(b)); %INSIDE
    
        
        Sigma=SSigma;
        D = 2*tau;
        G = sqrt(M*N);
        B = sqrt(G^2+M*N*Sigma^2);
        TT= floor(linspace(200,TTT,200));
        
    %% OUR
    parfor i=1:size(TT,2)
        i
        T = TT(i);
        
        alpha = 2*B*sqrt(T)/D;
        eta = G/sqrt(T)/D;
    
        X = zeros(N,M);
        XX = zeros(N,M);
        Y = X;
        Q = zeros(N,M);
        g = zeros(N,M);  
        
        for k = 1:(T-1)
            XX = XX+X;
            Q = Q+Y-X;
            g = sign(Y-W)+Sigma*randn(N,M);
            [Uq,Sq,Vq] = svd(Q);
            X = Uq(:,1) * tau * (Vq(:,1))';
            Y = (alpha*Y+eta*X-eta*Q-g)/(alpha+eta);  
        end
        XX = XX+X;
        X_OUR = XX/T;
        E_OUR(i)=sum(sum(abs(X_OUR-W)));
    end
    %% PGD
    parfor i=1:size(TT,2)
        i
        T = TT(i);
        beta = D/sqrt(T)/B/2;
        X = zeros(N,M);
        XX = zeros(N,M);  
        for k = 1:(T-1)
            XX = XX+X;
            g = sign(X-W)+Sigma*randn(N,M);
            yy = X - beta*g;  
            [Uy,Sy,Vy]=svd(yy);
            if sum(sum(Sy))<=tau
                X = yy;
            else
                X = Uy* max(0, Sy  -  f1(spdiags(Sy),tau)  )*Vy';
            end
        end
        XX = XX+X;
        X_PGD = XX/T;
        E_PGD(i) = sum(sum(abs(X_PGD-W)));
    end
    name = ['W=in_','N=',num2str(N),'_M=',num2str(M),'_sigma=',num2str(1000*Sigma)];
    save(name,"W","E_PGD","E_OUR","TT")
    %%
        
    end
end
%%
function s = f1(sigma,tau)
    for i=1:length(sigma)
        s = (sum(sigma(1:i))-tau)/i;
        if i<length(sigma)
            if (s>=sigma(i+1))
                break
            end
        end
    end
end

