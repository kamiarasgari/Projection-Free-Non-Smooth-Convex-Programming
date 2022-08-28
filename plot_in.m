clc
close all
clear all
%%
Sigma = [0,0.001,0.005,0.01,0.1,1,5,20];
for N=[5,10]
    for M = [N,2*N]

fig = figure('units','inch','position',[0,0,8,11]);
for i=[1,2,3,4,5,6,7,8]
    sigma=Sigma(i);
    
    name = ['W=in_','N=',num2str(N),'_M=',num2str(M),'_sigma=',num2str(1000*sigma),'.mat'];
    load(name)
    subplot(4,2,i);
    plot(TT,E_OUR,'b',TT,E_PGD,'r');
    ylabel('${f(\bar{x})-f(x^*)}$','interpreter','latex')
    xlabel('${T}$','interpreter','latex')
    legend('Our Algorithm','PGD/SGD');
    title(['\sigma =',num2str(sigma)])
    %ylim([0,0.01])
    xlim([5000,TT(end)])
    %set(gcf,'position',[0,0,1000,5000])
        name=['W=in_','N=',num2str(N),'_M=',num2str(M),'.eps'];
        saveas(fig,name,'epsc');
end 

    end
end
