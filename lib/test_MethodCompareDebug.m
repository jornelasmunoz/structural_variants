% MethodCompare_Child.m
% This code computes and plots the roc curves of child signal for
% 3 different cases, with the following methods:
% Novel Detection (new method)
% Regular Constrained (one individual at a time)
% GASV - Minimum Support ROC (this is based on a minimum threshold)

%This code uses Poisson Based SPIRAL to compare Poisson Data and Negative
%Binomial Data

clc;
clear;
close all

% =========================================================================
% =         Comparison of methods for the dataset: Simulated Data
% =========================================================================

% Load Simulated Data
% --------------------------------------------------------
addpath([genpath('/Users/jocelynornelasmunoz/Desktop/Research/structural_variants/'), ...
         genpath('/Users/jocelynornelas/iCloud Drive (Archive)/Desktop/UC Merced/Research/structural_variants/') ])

%filename = 'data/haploid_20pctNovel_10k_100n.mat';
%filename = 'data/haploid_20pctNovel_10k_100n_reproducedAPL.mat'; %reproduced data
%filename = 'lib/old/neg_binom_nov_p4_c4_5perNov.mat'; %Andrew's 5%nov 10^6n
filename = 'lib/old/neg_binom_nov_p4_c4_20perNov.mat'; %Andrew's 20%nov
%filename = 'lib/old/dummy_data.mat';
load(filename)


% Define tauvals to be used in loop
tauvals= 1;
for i=1:length(tauvals)
    
 t= tauvals(i);
 g= 10;
 tau=[t; t*g];

% tau(1) is for f_c and f_p
% tau(2) is for f_n (novel in child). Generally tau(1) << tau(2).

% =====================================================================
% =    Set up for Novel Method reconstruction    =
% =====================================================================
if contains(filename, 'neg_binom') || contains(filename, 'dummy')
    fprintf('Using Andrews data \n')
else    
    f_c_inh_neg_binom = double(f_h);
    f_c_nov_neg_binom = double(f_n); 
    f_p_neg_binom = double(f_p);
    y_c_neg_binom = double(s_c); 
    y_p_neg_binom = double(s_p);
    A_c_neg_binom = double(A_c);
    A_p_neg_binom = double(A_p);
end


% define signals
f_neg_binom = [f_c_inh_neg_binom; f_c_nov_neg_binom; f_p_neg_binom];
y_neg_binom =[y_c_neg_binom; y_p_neg_binom];


N = length(f_neg_binom);
n = N/3;

% set up the block diagonal matrix A
A_neg_binom = sparse(2*n, 3*n);
A_neg_binom(1:n,1:n) = A_c_neg_binom;
A_neg_binom(1:n,n+1:2*n) = A_c_neg_binom;
A_neg_binom(n+1:2*n, 2*n+1:3*n) = A_p_neg_binom;


% Setup function handles for computing A and A^T:
AT  = @(x) A_neg_binom'*x;
A   = @(x) A_neg_binom*x;

% set maximum number of iterations, tol, and when to print to
% screen
maxiter = 1;
tolerance = 1e-8;
verbose = 1;


% Simple initialization:
% initialization of f to start with
% AT(y) rescaled to a least-squares fit to the mean intensity
finit = (sum(sum(y_neg_binom)).*numel(AT(y_neg_binom)))...
    ./(sum(sum(AT(y_neg_binom))) .*sum(sum((AT(ones(size(y_neg_binom)))))))...
    .*AT(y_neg_binom);

%tic
%%%% POISSON Reconstruction
% Run the algorithm:
% Demonstrating all the options for our algorithm:
% fhatSPIRAL_1p1c_nov is the name of the reconstruction
% SPIRALTAP_novel - the new method
reg_params_all = [tau(1), tau(2), tau(1)];
noisetype = 'poisson'; 
subvectors = 3;
[fhatSPIRAL_1p1c_nov_neg_binom_SP, iterationsSPIRAL, objectiveSPIRAL,...
    reconerrorSPIRAL, cputimeSPIRAL] ...
    = NEBULA_debug(y_neg_binom,A_neg_binom,tau,reg_params_all,noisetype, subvectors,...
    'maxiter',maxiter,...
    'Initialization',finit,...
    'AT',AT,...
    'miniter',1,...
    'stopcriterion',3,...
    'tolerance',tolerance,...
    'alphainit',1,...
    'alphamin', 1e-30,...
    'alphamax', 1e30,...
    'alphaaccept',1e30,...
    'logepsilon',1e-10,...
    'saveobjective',1,...
    'savereconerror',1,...
    'savecputime',1,...
    'savesolutionpath',0,...
    'truth',f_neg_binom,...
    'verbose',verbose);
% %toc
% 
% % separate reconstruction for each true signal
% fhatSPIRAL_c_inh_neg_binom_SP = fhatSPIRAL_1p1c_nov_neg_binom_SP(1:n);
% fhatSPIRAL_c_nov_neg_binom_SP = fhatSPIRAL_1p1c_nov_neg_binom_SP(n+1:2*n);
% fhatSPIRAL_p_neg_binom_SP = fhatSPIRAL_1p1c_nov_neg_binom_SP(2*n+1:3*n);
% fhatSPIRAL_c_neg_binom_SP =fhatSPIRAL_c_inh_neg_binom_SP + fhatSPIRAL_c_nov_neg_binom_SP;

% =====================================================================
% =    Set up for Novel Method reconstruction    =
% =====================================================================

% f_poisson = [f_c_inh_neg_binom; f_c_nov_neg_binom; f_p_neg_binom];
% y_poisson =[y_c_neg_binom; y_p_neg_binom];
% 
% N = length(f_neg_binom);
% n = N/3;
% 
% % set up the block diagonal matrix A
% A_neg_binom = sparse(2*n, 3*n);
% A_neg_binom(1:n,1:n) = A_c_neg_binom;
% A_neg_binom(1:n,n+1:2*n) = A_c_neg_binom;
% A_neg_binom(n+1:2*n, 2*n+1:3*n) = A_p_neg_binom;
% 
% 
% % Setup function handles for computing A and A^T:
% AT  = @(x) A_neg_binom'*x;
% A   = @(x) A_neg_binom*x;
% 
% % set maximum number of iterations, tol, and when to print to
% % screen
% maxiter = 1000;
% %maxiter = 500;
% tolerance = 1e-8;
% %tolerance = 1e-4;
% verbose = 100;
% %verbose = 30;

% Simple initialization:
% initialization of f to start with
% AT(y) rescaled to a least-squares fit to the mean intensity
% finit = (sum(sum(y_neg_binom)).*numel(AT(y_neg_binom)))...
%     ./(sum(sum(AT(y_neg_binom))) .*sum(sum((AT(ones(size(y_neg_binom)))))))...
%     .*AT(y_neg_binom);

%tic
% Run the algorithm:
% Demonstrating all the options for our algorithm:
% fhatSPIRAL_1p1c_nov is the name of the reconstruction
% SPIRALTAP_novel - the new method
% noisetype = 'negative binomial'; 
% [fhatSPIRAL_1p1c_nov_neg_binom_SPNB, iterationsSPIRAL, objectiveSPIRAL,...
%     reconerrorSPIRAL, cputimeSPIRAL] ...
%     = NEBULA_debug(y_neg_binom,A_neg_binom,tau,reg_params_all,noisetype, subvectors,... %SPIRALTAP_NegBinom_Novel(y_neg_binom,A_neg_binom,tau,...
%     'maxiter',maxiter,...
%     'Initialization',finit,...
%     'AT',AT,...
%     'miniter',1,...
%     'stopcriterion',3,...
%     'tolerance',tolerance,...
%     'alphainit',1,...
%     'alphamin', 1e-30,...
%     'alphamax', 1e30,...
%     'alphaaccept',1e30,...
%     'logepsilon',1e-10,...
%     'saveobjective',1,...
%     'savereconerror',1,...
%     'savecputime',1,...
%     'savesolutionpath',0,...
%     'truth',f_neg_binom,...
%     'verbose',verbose);
% %toc
% 
% % separate reconstruction for each true signal
% fhatSPIRAL_c_inh_neg_binom_SPNB = fhatSPIRAL_1p1c_nov_neg_binom_SPNB(1:n);
% fhatSPIRAL_c_nov_neg_binom_SPNB = fhatSPIRAL_1p1c_nov_neg_binom_SPNB(n+1:2*n);
% fhatSPIRAL_p_neg_binom_SPNB = fhatSPIRAL_1p1c_nov_neg_binom_SPNB(2*n+1:3*n);
% fhatSPIRAL_c_neg_binom_SPNB =fhatSPIRAL_c_inh_neg_binom_SPNB + fhatSPIRAL_c_nov_neg_binom_SPNB;
% 

% %-----------------------------------------------------
% % THRESHOLDING OF SIGNALS
% %-----------------------------------------------------
% 
% thresh = linspace(-.001,1.001,100); % make thresholding vector for ROC curve
% 
% 
% % ------------------------ Novel Method ROC --------------------------
% 
% % ROC Information for Negative Binomial (APL)
% % NOVroc_gen( reconstruction, true signal, threshold)
% trueSIGNAL_neg_binom = [f_c_inh_neg_binom + f_c_nov_neg_binom; f_p_neg_binom];
% 
% %Child and Parent Signal
% fhatSPIRAL_nov_cp_neg_binom_SPNB = [fhatSPIRAL_c_neg_binom_SPNB; fhatSPIRAL_p_neg_binom_SPNB];
% [FPRv_cp_neg_binom_SPNB,TPRv_cp_neg_binom_SPNB] = NOVroc_gen(fhatSPIRAL_nov_cp_neg_binom_SPNB, trueSIGNAL_neg_binom, thresh);
% Area_nov_cp_neg_binom_SPNB = trapz(FPRv_cp_neg_binom_SPNB,TPRv_cp_neg_binom_SPNB);
% 
% %Parent Signal
% [FPRv_p_neg_binom_SPNB, TPRv_p_neg_binom_SPNB] = NOVroc_gen(fhatSPIRAL_nov_cp_neg_binom_SPNB(n+1:2*n),trueSIGNAL_neg_binom(n+1:2*n),thresh);
% Area_nov_p_neg_binom_SPNB = trapz(FPRv_p_neg_binom_SPNB, TPRv_p_neg_binom_SPNB);
% 
% %Child Signal
% [FPRv_c_neg_binom_SPNB, TPRv_c_neg_binom_SPNB] = NOVroc_gen(fhatSPIRAL_nov_cp_neg_binom_SPNB(1:n),trueSIGNAL_neg_binom(1:n),thresh);
% Area_nov_c_neg_binom_SPNB = trapz(FPRv_c_neg_binom_SPNB, TPRv_c_neg_binom_SPNB);
% 
% 
% 
% %Child and Parent Signal
% fhatSPIRAL_nov_cp_neg_binom_SP = [fhatSPIRAL_c_neg_binom_SP; fhatSPIRAL_p_neg_binom_SP];
% [FPRv_cp_neg_binom_SP,TPRv_cp_neg_binom_SP] = NOVroc_gen(fhatSPIRAL_nov_cp_neg_binom_SP, trueSIGNAL_neg_binom, thresh);
% Area_nov_cp_neg_binom_SP = trapz(FPRv_cp_neg_binom_SP,TPRv_cp_neg_binom_SP);
% 
% %Parent Signal
% [FPRv_p_neg_binom_SP, TPRv_p_neg_binom_SP] = NOVroc_gen(fhatSPIRAL_nov_cp_neg_binom_SP(n+1:2*n), trueSIGNAL_neg_binom(n+1:2*n), thresh);
% Area_nov_p_neg_binom_SP = trapz(FPRv_p_neg_binom_SP,TPRv_p_neg_binom_SP);
% 
% %Child Signal
% [FPRv_c_neg_binom_SP, TPRv_c_neg_binom_SP] = NOVroc_gen(fhatSPIRAL_nov_cp_neg_binom_SP(1:n), trueSIGNAL_neg_binom(1:n), thresh);
% Area_nov_c_neg_binom_SP = trapz(FPRv_c_neg_binom_SP,TPRv_c_neg_binom_SP);
% 
% AreavecSPNB=[Area_nov_cp_neg_binom_SPNB; Area_nov_p_neg_binom_SPNB; Area_nov_c_neg_binom_SPNB];
% AreavecSP=[Area_nov_cp_neg_binom_SP; Area_nov_p_neg_binom_SP; Area_nov_c_neg_binom_SP];
% 
% % %Child and Parent Plot
% % figure
% % plot(FPRv_cp_neg_binom_SPNB,TPRv_cp_neg_binom_SPNB,'r-.', 'LineWidth',3);
% % hold on
% % plot(FPRv_cp_neg_binom_SP,TPRv_cp_neg_binom_SP,'b--.', 'LineWidth',3);
% % xlabel('False Positive Rate','FontSize',16);
% % ylabel('True Positive Rate','FontSize',16);
% % title({'ROC Curves for C1 Signal',['tau = ' num2str(tau(1)),...
% %     '  tau_n = ' num2str(tau(2))],['Reconstructed from Child and Parent']},'FontSize',16)
% % %legend1=legend(['Novel Method Negative Binomial SPIRAL AUC=' num2str(AreavecSPNB(1))],['Novel Method Poisson SPIRAL AUC=' num2str(AreavecSP(1))]);
% % %legend1=legend( 'Novel Method Negative Binomial SPIRAL','Novel Method Poisson SPIRAL');
% % legend1=legend('NEBULA','SPIRAL');
% % set(legend1,...
% %     'Location','southeast',...
% %     'FontSize',14);
% % % set(legend2,...
% % %     'Location','east',...
% % %     'FontSize',14);
% % set(gca,'fontsize',14);
% % hold off
% % 
% % 
% % 
% % 
% % 
% % figure
% % %Parent Plot
% % plot(FPRv_p_neg_binom_SPNB,TPRv_p_neg_binom_SPNB,'r-.', 'LineWidth',2);
% % hold on
% % plot(FPRv_p_neg_binom_SP,TPRv_p_neg_binom_SP,'b--.', 'LineWidth',2);
% % xlabel('False Positive Rate','FontSize',16);
% % ylabel('True Positive Rate','FontSize',16);
% % title({'ROC Curves for C1 Signal',['tau = ' num2str(tau(1)),...
% %     '  tau_n = ' num2str(tau(2))],['Reconstructed from Parent']},'FontSize',16)
% % %legend1=legend(['Novel Method Negative Binomial SPIRAL AUC=' num2str(AreavecSPNB(2))],['Novel Method Poisson SPIRAL AUC=' num2str(AreavecSP(2))]);
% % %legend1=legend( 'Novel Method Negative Binomial SPIRAL','Novel Method Poisson SPIRAL');
% % legend2=legend(['AUC=' num2str(AreavecSPNB(2))],['AUC=' num2str(AreavecSP(2))]);
% % %legend1=legend('NEBULA','SPIRAL');
% % %set(legend1,...
% % %    'Location','southeast',...
% % %    'FontSize',14);
% % set(legend2,...
% %     'Location','east',...
% %     'FontSize',14);
% % set(gca,'fontsize',14);
% % hold off
% % 
% % 
% % 
% % figure
% % %Child Plot
% % plot(FPRv_c_neg_binom_SPNB,TPRv_c_neg_binom_SPNB,'r-.', 'LineWidth',3);
% % % hold on
% % % plot(FPRv_c_neg_binom_SP,TPRv_c_neg_binom_SP,'b--.', 'LineWidth',3);
% % xlabel('False Positive Rate','FontSize',16);
% % ylabel('True Positive Rate','FontSize',16);
% % %title({'ROC Curves for C1 Signal',['tau = ' num2str(tau(1)),...
% %     %'  tau_n = ' num2str(tau(2))],['Reconstructed from Child']},'FontSize',16)
% % %legend1=legend(['Novel Method Negative Binomial SPIRAL AUC=' num2str(AreavecSPNB(3))],['Novel Method Poisson SPIRAL AUC=' num2str(AreavecSP(3))]);
% % %legend2=legend(['AUC=' num2str(AreavecSPNB(3))],['AUC=' num2str(AreavecSP(3))]);
% % %legend1=legend('NEBULA','SPIRAL');
% % % set(legend1,...
% % % 'Location','southeast',...
% % %     'FontSize',14);
% % % set(legend2,...
% % %     'Location','east',...
% % %     'FontSize',14);
% % set(gca,'fontsize',14);
% % hold off
% 
% %------------------Jocelyn edits -----------
% 
% % % overall plot
% % figure
% % [X_n,Y_n,T_n,AUC_n] = perfcurve(f_neg_binom,fhatSPIRAL_1p1c_nov_neg_binom_SPNB, 1); 
% % [X_s,Y_s,T_s,AUC_s] = perfcurve(f_neg_binom,fhatSPIRAL_1p1c_nov_neg_binom_SP, 1);
% % plot(X_n,Y_n, '-r', 'LineWidth',2); hold on
% % plot(X_s, Y_s, '--b', 'LineWidth',1.5)
% % legend(strcat('NEBULA = ', num2str(AUC_n)), strcat('SPIRAL  = ', num2str(AUC_s)), 'FontSize',12)
% % xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
% % title({'ROC Curves for total reconstruction',['\tau = ' num2str(tau(1)),...
% %     ' \gamma = ' num2str(tau(2))]},'FontSize',16)
end