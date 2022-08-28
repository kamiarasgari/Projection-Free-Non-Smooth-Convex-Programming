clc
close all
clear all
%%
Sigma = [0,1,5,10];
NN=[5,5,10,10];
MM=[5,10,10,20];

for k=1:4
    N=NN(k);
    M=MM(k);
        fig = figure;
        for i=[1,2,3,4]
            sigma=Sigma(i);
            
            name = ['W=out_','N=',num2str(N),'_M=',num2str(M),'_sigma=',num2str(1000*sigma),'.mat'];
            load(name)
            subplot(3,2,i);
            plot(TT,E_OUR,'b',TT,E_PGD,'r');
            ylabel('${f(\bar{x})}$','interpreter','latex')
            xlabel('${T}$','interpreter','latex')
            legend('Our Algorithm','PGD/SGD');
            title(['\sigma =',num2str(sigma)])
            %ylim(YY(k,:))
            xlim([5000,TT(end)])
            
            set(gcf,'position',[0,0,1000,2000])
            
            name=['W=out_','N=',num2str(N),'_M=',num2str(M),'.eps'];
            saveas(fig,name,'epsc');
        end 
end
