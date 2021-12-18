thresh = linspace(0,1,1000);
[FPR_SPIRAL, TPR_SPIRAL] = NOVroc_gen(fhat_SPIRAL, f_true, thresh);
AUC_SPIRAL = trapz(FPR_SPIRAL,TPR_SPIRAL)
[FPR_NEBULA, TPR_NEBULA] = NOVroc_gen(fhat_NEBULA, f_true, thresh);
AUC_NEBULA = trapz(FPR_NEBULA,TPR_NEBULA)

figure
plot(FPR_SPIRAL,TPR_SPIRAL,'r-.', 'LineWidth',3);
hold on
plot(FPR_NEBULA,TPR_NEBULA,'b--.', 'LineWidth',3);
xlabel('False Positive Rate','FontSize',16);
ylabel('True Positive Rate','FontSize',16);
% title({'ROC Curves for C1 Signal',['tau = ' num2str(tau(1)),...
%     '  tau_n = ' num2str(tau(2))],['Reconstructed from Child and Parent']},'FontSize',16)
%legend1=legend(['Novel Method Negative Binomial SPIRAL AUC=' num2str(AreavecSPNB(1))],['Novel Method Poisson SPIRAL AUC=' num2str(AreavecSP(1))]);
%legend1=legend( 'Novel Method Negative Binomial SPIRAL','Novel Method Poisson SPIRAL');
legend1=legend('SPIRAL','NEBULA');
set(legend1,...
    'Location','southeast',...
    'FontSize',14);
% set(legend2,...
%     'Location','east',...
%     'FontSize',14);
set(gca,'fontsize',14);
hold off