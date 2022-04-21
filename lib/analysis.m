% Analyze results of given file


% Sparsity - i.e. nonzeros
[y_St,x_St] = sparsity_search(fhatSPIRAL);
[y_Nt,x_Nt] = sparsity_search(fhatNEBULA);
[y_Sp,x_Sp] = sparsity_search(fhatSPIRAL_p);
[y_Np,x_Np] = sparsity_search(fhatNEBULA_p);
[y_Sc,x_Sc] = sparsity_search(fhatSPIRAL_c);
[y_Nc,x_Nc] = sparsity_search(fhatNEBULA_c);

[y_T,x_T] = sparsity_search(f_true);
[y_P,x_P] = sparsity_search(f_p);
[y_C,x_C] = sparsity_search(f_c);

figure('Position', [500 500 1500 400])

subplot(1,3,1);
loglog(x_T,y_T,'--'); hold on; 
loglog(x_Nt,y_Nt,'-', 'LineWidth',1); 
loglog(x_St,y_St,'-', 'LineWidth',1); hold off
%xlabel('Threshold','FontSize',16); 
ylabel('Nonzero count','FontSize',16);
title({'Total Signal Sparsity'},'FontSize',16)
ylim([0, 10^3]);

subplot(1,3,2);
loglog(x_P,y_P,'--'); hold on; 
loglog(x_Sp,y_Sp,'-', 'LineWidth',1);
loglog(x_Np,y_Np,'-', 'LineWidth',1); hold off;
xlabel('Threshold','FontSize',16); %ylabel('Nonzero count','FontSize',16);
title({'Parent Signal Sparsity'},'FontSize',16)
ylim([0, 10^3])

subplot(1,3,3);
loglog(x_C,y_C,'--'); hold on; 
loglog(x_Sc,y_Sc,'-', 'LineWidth',1);
loglog(x_Nc,y_Nc,'-', 'LineWidth',1); hold off;
%xlabel('Threshold','FontSize',16); ylabel('Nonzero count','FontSize',16);
title({'Child Signal Sparsity'},'FontSize',16);
legend('Truth','NEBULA', 'SPIRAL', 'FontSize',14, 'Location','southwest');
ylim([0, 10^3])