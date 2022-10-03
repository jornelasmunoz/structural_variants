% Analyze results of given file
%clc; clear all; close all;
addpath([genpath("/Users/jocelynornelasmunoz/Desktop/Research/structural_variants/"), ...
         genpath("/Users/jocelynornelas/iCloud Drive (Archive)/Desktop/UC Merced/Research/structural_variants/"),...
         genpath("/home/jornelasmunoz/structural_variants/")])

% files = dir("results/10000n_5k/2Lp_4Lc/");
% idx = 7;
% result_file = files(idx+3).name;
% %result_file = "results/10000n_5k/2Lp_4Lc/0.01tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/0.1tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/1tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/10.1tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/100tau_500gamma_RESULTS.mat";
% result_file = "results/10000n_5k/2Lp_4Lc/1000tau_500gamma_RESULTS.mat";

% load(result_file)
% disp(result_file)

%save results
% sub_folder = char(fileparts(result_file));
% fig_folder = strcat('figs/',sub_folder(6:end));
% if ~exist(fig_folder, 'dir')
%     disp('Making folder');disp(fig_folder)
%     mkdir(fig_folder)
%end
%save_path = sprintf(strcat(fig_folder,'/%stau_%sgamma_RESULTS.mat'), num2str(tau), num2str(gamma));
%save(save_path)
%%
%-------------------------------------------------------------------------
% ---------------------  Sparsity - i.e. nonzeros ------------------------
%-------------------------------------------------------------------------
% [y_St,x_St] = sparsity_search(fhatSPIRAL);
% [y_Nt,x_Nt] = sparsity_search(fhatNEBULA);
% [y_Sp,x_Sp] = sparsity_search(fhatSPIRAL_p);
% [y_Np,x_Np] = sparsity_search(fhatNEBULA_p);
% [y_Sc,x_Sc] = sparsity_search(fhatSPIRAL_c);
% [y_Nc,x_Nc] = sparsity_search(fhatNEBULA_c);
% 
% [y_T,x_T] = sparsity_search(f_true);
% [y_P,x_P] = sparsity_search(f_p);
% [y_C,x_C] = sparsity_search(f_c);
% 
% figure('Position', [500 500 1500 400])
% 
% subplot(1,3,1);
% loglog(x_T,y_T,'--'); hold on; 
% loglog(x_Nt,y_Nt,'-', 'LineWidth',1); 
% loglog(x_St,y_St,'-', 'LineWidth',1); hold off
% %xlabel('Threshold','FontSize',16); 
% ylabel('Nonzero count','FontSize',16);
% title({'Total Signal Sparsity'},'FontSize',16);
% legend('Truth','NEBULA', 'SPIRAL', 'FontSize',14, 'Location','southwest');
% ylim([0, 10^4]);
% 
% subplot(1,3,2);
% loglog(x_P,y_P,'--'); hold on; 
% loglog(x_Sp,y_Sp,'-', 'LineWidth',1);
% loglog(x_Np,y_Np,'-', 'LineWidth',1); hold off;
% xlabel('Threshold','FontSize',16); %ylabel('Nonzero count','FontSize',16);
% title({'Parent Signal Sparsity'},'FontSize',16);
% legend('Truth','NEBULA', 'SPIRAL', 'FontSize',14, 'Location','southwest');
% ylim([0, 10^4])
% 
% subplot(1,3,3);
% loglog(x_C,y_C,'--'); hold on; 
% loglog(x_Sc,y_Sc,'-', 'LineWidth',1);
% loglog(x_Nc,y_Nc,'-', 'LineWidth',1); hold off;
% %xlabel('Threshold','FontSize',16); ylabel('Nonzero count','FontSize',16);
% title({'Child Signal Sparsity'},'FontSize',16);
% legend('Truth','NEBULA', 'SPIRAL', 'FontSize',14, 'Location','southwest');
% ylim([0, 10^4])

%-------------------------------------------------------------------------
%            ROC curves
%-------------------------------------------------------------------------
% overall plot
figure('Position', [500 500 1500 400])
subplot(2,4,1);
plot(X_nt,Y_nt, '-r', 'LineWidth',2); hold on
plot(X_st, Y_st, '--b', 'LineWidth',1.5);
legend(strcat('NEBULA = ', num2str(AUC_nt)), strcat('SPIRAL  = ', num2str(AUC_st)), 'FontSize',14, 'Location', 'southeast');
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
title({'ROC Curves for total reconstruction',['\tau = ' num2str(tau),...
    ' \gamma = ' num2str(gamma)]},'FontSize',16)


%parent reconstruction - homozygous
subplot(2,4,2);
plot(X_np_2,Y_np_2, '-r', 'LineWidth',2); hold on
plot(X_sp_2, Y_sp_2, '--b', 'LineWidth',1.5);
legend(strcat('NEBULA = ', num2str(AUC_np_2)), strcat('SPIRAL  = ', num2str(AUC_sp_2)), 'FontSize',14, 'Location', 'southeast');
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
title({'ROC Curves for parent reconstruction - homozygous',['\tau = ' num2str(tau),...
    ' \gamma = ' num2str(gamma)]},'FontSize',16)

%child reconstruction - inherited homozygous
subplot(2,4,3);
plot(X_nc_2h,Y_nc_2h, '-r', 'LineWidth',2); hold on
plot(X_sc_2h, Y_sc_2h, '--b', 'LineWidth',1.5);
legend(strcat('NEBULA = ', num2str(AUC_nc_2h)), strcat('SPIRAL  = ', num2str(AUC_sc_2h)), 'FontSize',14, 'Location', 'southeast');
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
title({'ROC Curves for child reconstruction - inherited homozygous',['\tau = ' num2str(tau),...
    ' \gamma = ' num2str(gamma)]},'FontSize',16)

%child reconstruction - novel homozygous
subplot(2,4,4);
plot(X_nc_2n,Y_nc_2n, '-r', 'LineWidth',2); hold on
plot(X_sc_2n, Y_sc_2n, '--b', 'LineWidth',1.5);
legend(strcat('NEBULA = ', num2str(AUC_nc_2n)), strcat('SPIRAL  = ', num2str(AUC_sc_2n)), 'FontSize',14, 'Location', 'southeast');
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
title({'ROC Curves for child reconstruction - novel homozygous',['\tau = ' num2str(tau),...
    ' \gamma = ' num2str(gamma)]},'FontSize',16)

%parent reconstruction - heterozygous
subplot(2,4,6);
plot(X_np_1,Y_np_1, '-r', 'LineWidth',2); hold on
plot(X_sp_1, Y_sp_1, '--b', 'LineWidth',1.5);
legend(strcat('NEBULA = ', num2str(AUC_np_1)), strcat('SPIRAL  = ', num2str(AUC_sp_1)), 'FontSize',14, 'Location', 'southeast');
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
title({'ROC Curves for parent reconstruction - heterozygous',['\tau = ' num2str(tau),...
    ' \gamma = ' num2str(gamma)]},'FontSize',16)

%child reconstruction - inherited heterozygous
subplot(2,4,7);
plot(X_nc_1h,Y_nc_1h, '-r', 'LineWidth',2); hold on
plot(X_sc_1h, Y_sc_1h, '--b', 'LineWidth',1.5);
legend(strcat('NEBULA = ', num2str(AUC_nc_1h)), strcat('SPIRAL  = ', num2str(AUC_sc_1h)), 'FontSize',14, 'Location', 'southeast');
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
title({'ROC Curves for child reconstruction - inherited heterozygous',['\tau = ' num2str(tau),...
    ' \gamma = ' num2str(gamma)]},'FontSize',16)

%child reconstruction - novel heterozygous
subplot(2,4,8);
plot(X_nc_1n,Y_nc_1n, '-r', 'LineWidth',2); hold on
plot(X_sc_1n, Y_sc_1n, '--b', 'LineWidth',1.5);
legend(strcat('NEBULA = ', num2str(AUC_nc_1n)), strcat('SPIRAL  = ', num2str(AUC_sc_1n)), 'FontSize',14, 'Location', 'southeast');
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
title({'ROC Curves for child reconstruction - novel heterozygous',['\tau = ' num2str(tau),...
    ' \gamma = ' num2str(gamma)]},'FontSize',16)
%-------------------------------------------------------------------------
%------------------   Precision-Recall Curves  ---------------------------
%-------------------------------------------------------------------------
% total reconstruction
[p_N,r_N,no_skill_N,th_N] = compute_precision_recall(f_true,fhatNEBULA);
[p_S,r_S,no_skill_S,th_S] = compute_precision_recall(f_true,fhatSPIRAL);

% parent reconstruction -- homogeneous (2) and heterogeneous(1)
[p_Np_2,r_Np_2,no_skill_Np_2,th_Np_2] = compute_precision_recall(z_p,fhatNEBULA_zp);
[p_Np_1,r_Np_1,no_skill_Np_1,th_Np_1] = compute_precision_recall(y_p,fhatNEBULA_zp);
[p_Sp_2,r_Sp_2,no_skill_Sp_2,th_Sp_2] = compute_precision_recall(z_p,fhatSPIRAL_zp);
[p_Sp_1,r_Sp_1,no_skill_Sp_1,th_Sp_1] = compute_precision_recall(y_p,fhatSPIRAL_yp);

% child reconstruction -- inherited
[p_Nc_2h,r_Nc_2h,no_skill_Nc_2h,th_Nc_2h] = compute_precision_recall(z_h,fhatNEBULA_zh);
[p_Nc_1h,r_Nc_1h,no_skill_Nc_1h,th_Nc_1h] = compute_precision_recall(y_h,fhatNEBULA_yh);
[p_Sc_2h,r_Sc_2h,no_skill_Sc_2h,th_Sc_2h] = compute_precision_recall(z_h,fhatSPIRAL_zh);
[p_Sc_1h,r_Sc_1h,no_skill_Sc_1h,th_Sc_1h] = compute_precision_recall(y_h,fhatSPIRAL_yh);

% child reconstruction -- novel
[p_Nc_2n,r_Nc_2n,no_skill_Nc_2n,th_Nc_2n] = compute_precision_recall(z_n,fhatNEBULA_zn);
[p_Nc_1n,r_Nc_1n,no_skill_Nc_1n,th_Nc_1n] = compute_precision_recall(y_n,fhatNEBULA_yn);
[p_Sc_2n,r_Sc_2n,no_skill_Sc_2n,th_Sc_2n] = compute_precision_recall(z_n,fhatSPIRAL_zn);
[p_Sc_1n,r_Sc_1n,no_skill_Sc_1n,th_Sc_1n] = compute_precision_recall(y_n,fhatSPIRAL_yn);


figure('Position', [500 500 1500 400])
% total reconstruction
subplot(2,4,1);
plot([r_N;0],[p_N;1],'LineWidth',1); hold on;
plot([r_S;0],[p_S;1],'LineWidth',1);
plot([0;1], [no_skill_N;no_skill_N], '--r'); hold off
title({'Total Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');
ylabel('Precision','FontSize',16)

% parent - homozygous
subplot(2,4,2);
plot([r_Np_2;0],[p_Np_2,;1],'LineWidth',1); hold on;
plot([r_Sp_2;0],[p_Sp_2;1],'LineWidth',1);
plot([0;1], [no_skill_Np_2;no_skill_Np_2], '--r'); hold off
title({'Homozygous Parent Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');
xlabel('Recall','FontSize',16)

% parent - heterozygous
subplot(2,4,6);
plot([r_Np_1;0],[p_Np_1,;1],'LineWidth',1); hold on;
plot([r_Sp_1;0],[p_Sp_1;1],'LineWidth',1);
plot([0;1], [no_skill_Np_2;no_skill_Np_2], '--r'); hold off
title({'Heterozygous Parent Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');
xlabel('Recall','FontSize',16)

% child - homozygous
subplot(2,4,3);
plot([r_Nc_2h;0],[p_Nc_2h;1],'LineWidth',1); hold on;
plot([r_Sc_2h;0],[p_Sc_2h;1],'LineWidth',1);
plot([0;1], [no_skill_Nc_2h;no_skill_Nc_2h], '--r'); hold off
title({'Inherited Homozygous Child Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');

subplot(2,4,4);
plot([r_Nc_2n;0],[p_Nc_2n;1],'LineWidth',1); hold on;
plot([r_Sc_2n;0],[p_Sc_2n;1],'LineWidth',1);
plot([0;1], [no_skill_Nc_2n;no_skill_Nc_2n], '--r'); hold off
title({'Novel Homozygous Child Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');

% child - heterozygous
subplot(2,4,7);
plot([r_Nc_1h;0],[p_Nc_1h;1],'LineWidth',1); hold on;
plot([r_Sc_1h;0],[p_Sc_1h;1],'LineWidth',1);
plot([0;1], [no_skill_Nc_1h;no_skill_Nc_1h], '--r'); hold off
title({'Inherited Heterozygous Child Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');

subplot(2,4,8);
plot([r_Nc_1n;0],[p_Nc_1n;1],'LineWidth',1); hold on;
plot([r_Sc_1n;0],[p_Sc_1n;1],'LineWidth',1);
plot([0;1], [no_skill_Nc_1n;no_skill_Nc_1n], '--r'); hold off
title({'Novel Heterozygous Child Signal Precision-Recall'},'FontSize',16);
legend('NEBULA', 'SPIRAL','No Skill', 'FontSize',14, 'Location','southwest');

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% order [tau, gamma, AUC_nt, AUC_st, AUC_np, AUC_sp, AUC_nc, AUC_sc];
% params_AUCs = reshape(params_AUCs, 6*6,8);
% params_AUCs = params_AUCs(params_AUCs(:,1) == tau,:);
% %params_AUCs(params_AUCs ~=0);
%params_AUCs = reshape(params_AUCs, length(params_AUCs)/8,8);

% figure()
% scatter(params_AUCs(:,3),params_AUCs(:,4),25, 'k');hold on
% scatter(params_AUCs(:,5),params_AUCs(:,6),25, 'r');
% scatter(params_AUCs(:,7),params_AUCs(:,8),25, 'b');
% plot([0.5,1],[0.5,1], '--'); hold off
% xlabel('NEBULA'); ylabel('SPIRAL')
% legend('Total signal', 'Parent signal', 'Child signal', 'Location','southeast','FontSize',14);
% title('ROC AUC comparison between NEBULA and SPIRAL','FontSize',16)
% 
% figure
% h = heatmap(params_AUCs(:,3:end),'FontSize',13);
% h.YDisplayLabels = string(params_AUCs(:,2));
% h.XDisplayLabels = ["N_{total}", "S_{total}", "N_{parent}", "S_{parent}", "N_{child}", "S_{child}"];
% h.Title = strcat('ROC AUC comparisons for \tau = ', num2str(tau));