clc
clear all
close all
%%
steps = 400;
SSigma=[0,1,5,10,30];
TTT=[10000,20000,40000,50000,70000];
for n=[10,100,250,500]
    w = 1*(1-2*rand(n,1));
    X_opt = satlin(w')-satlin(-w');
for s=[1,2,3,4,5]

Sigma=SSigma(s);
D = 2*sqrt(n);
G = sqrt(n);
B = sqrt(G^2+n*Sigma^2);
%%
TT= floor(linspace(200,TTT(s),steps));
for i=1:size(TT,2)
    i
    T = TT(i);
    %%
    alpha = 2*B*sqrt(T)/D;
    eta = G/sqrt(T)/D;
    beta = D/sqrt(T)/B/2;

    x = zeros(n,T);
    y = zeros(n,1);
    Q = zeros(n,1);
    %% OURS
    x(:,1) = zeros(n,1);
    y = x(:,1);
    
    for k = 1:(T-1)
        Q = Q+y-x(:,k);
        g = grad(y,w,n,Sigma);
        
        x(:,k+1) = sign(Q);
        y = (alpha*y+eta*x(:,k+1)-eta*Q-g)/(alpha+eta);  
    end
    X_T_OUR = mean(x');
    E_OUR(i) = norm1(w'-X_T_OUR)-norm1(w'-X_opt);
    %% PGD/SGD
    x = zeros(n,T);
    for k = 1:(T-1)
        g = grad(x(:,k),w,n,Sigma);
        yy = x(:,k) - beta*g;    
        x(:,k+1) = satlin(yy')-satlin(-yy'); 
    end
    X_T_1 = mean(x');
    E_1(i) = norm1(w'-X_T_1)-norm1(w'-X_opt);
    %%
end


name = ['W=in_','dim=',num2str(n),'_sigma=',num2str(Sigma)];
save(name,"X_opt","E_1","E_OUR","X_T_1","X_T_OUR","TT")

end
end




%OUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUTOUT



%%

for n=[10,100,250,500]
    
    w =3*(1-2*rand(n,1));

    X_opt = satlin(w')-satlin(-w');
for s=[1,2,3,4,5]

Sigma=SSigma(s);
D = 2*sqrt(n);
G = sqrt(n);
B = sqrt(G^2+n*Sigma^2);
%%
TT= floor(linspace(200,TTT(s),steps));
for i=1:size(TT,2)
    i
    T = TT(i);
    %%
    alpha = 2*B*sqrt(T)/D;
    eta = G/sqrt(T)/D;
    beta = D/sqrt(T)/B/2;

    x = zeros(n,T);
    y = zeros(n,1);
    Q = zeros(n,1);
    %% OURS
    x(:,1) = zeros(n,1);
    y = x(:,1);
    
    for k = 1:(T-1)
        Q = Q+y-x(:,k);
        g = grad(y,w,n,Sigma);
        
        x(:,k+1) = sign(Q);
        y = (alpha*y+eta*x(:,k+1)-eta*Q-g)/(alpha+eta);  
    end
    X_T_OUR = mean(x');
    E_OUR(i) = norm1(w'-X_T_OUR)-norm1(w'-X_opt);
    %% PGD/SGD
    x = zeros(n,T);
    for k = 1:(T-1)
        g = grad(x(:,k),w,n,Sigma);
        yy = x(:,k) - beta*g;    
        x(:,k+1) = satlin(yy')-satlin(-yy'); 
    end
    X_T_1 = mean(x');
    E_1(i) = norm1(w'-X_T_1)-norm1(w'-X_opt);
    %%
end


name = ['W=out_','dim=',num2str(n),'_sigma=',num2str(Sigma)];
save(name,"X_opt","E_1","E_OUR","X_T_1","X_T_OUR","TT")

end
end















%%
% plot(TT,E_OUR,'b',TT,E_1,'r')
% legend('Our Algorithm','PGD/SGD')

%%
function g = grad(x,C,dim,Sigma)
    g = sign(x-C)+Sigma*randn(dim,1);
end
%%
function N = norm1(x)
    N=sum(sum(abs(x)));
end
