clc
close all
clear all
%%
Sigma = [0,1,10,30];
for dim=[10,100,250,500]
    fig = figure;
    for i=1:4
        sigma=Sigma(i);
        
        name=['W=in_dim=',num2str(dim),'_sigma=',num2str(sigma),'.mat'];
        load(name)
        subplot(2,2,i);
        plot(TT,E_OUR,'b',TT,E_1,'r');
        ylabel('${f(\bar{x})-f(x^*)}$','interpreter','latex')
        xlabel('${T}$','interpreter','latex')
        legend('Our Algorithm','PGD/SGD');
        title(['$\sigma =$',num2str(sigma)])
        xlim([0,TT(end)])
    
        name=['W=in_dim=',num2str(dim),'.eps'];
        saveas(fig,name,'epsc');
    end    
end
