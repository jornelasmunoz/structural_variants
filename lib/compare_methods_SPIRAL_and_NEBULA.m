%compare_methods_SPIRAL_and_NEBULA.m
% This script compares the simulated diploid, 1 parent-1 child data using
% SPIRAL-- Poisson-based method and
% NEBULA-- Negative Binomial-based method

% This code was adapted from Melissa Spence and Mario Banuelos's SV work

clc;
clear;
close all

% =========================================================================
% =         Preparation of data variables: Simulated Data
% =========================================================================


addpath([genpath('/Users/jocelynornelasmunoz/Desktop/Research/structural_variants/'), ...
         genpath('/Users/jocelynornelas/iCloud Drive (Archive)/Desktop/UC Merced/Research/structural_variants/') ])
% -------------------------  Load Simulated Data  -------------------------
%filename = 'data/diploid_5pctNovel_5k_100n.mat'; 
%filename = 'data/diploid_5pctNovel_500k_100000n.mat';
%filename = 'data/diploid_5pctNovel_50k_10000n.mat'; 
%filename = 'data/diploid_2pctNovel_10k_100n.mat';
filename = 'dummy_2pctNovel_2k_6n.mat';

%filename = 'data/haploid_5pctNovel_10k_100n.mat';
%filename = 'data/haploid_20pctNovel_10k_100n_reproducedAPL.mat'; %reproduced data
%filename = 'lib/old/neg_binom_nov_p4_c4_5perNov.mat'; %Andrew's 5%nov 10^6n
%filename = 'lib/old/neg_binom_nov_p4_c4_20perNov.mat'; %Andrew's 20%nov
%filename = 'lib/old/dummy_data.mat';
load(filename)

% By default, kind is diploid
% Will be changed if using Andrew's data, or manually 
kind = 'diploid';
if contains(filename, 'neg_binom') || contains(filename, 'dummy_data')
    fprintf('Using Andrews data \n')
    kind = 'haploid';
    f_p = f_p_neg_binom;
    f_h = f_c_inh_neg_binom;
    f_n = f_c_nov_neg_binom;
    f_c = f_h + f_n;
    s_p = y_p_neg_binom;
    s_c = y_c_neg_binom;
    A_c = A_c_neg_binom;
    A_p = A_p_neg_binom;
elseif contains(filename,'haploid')
    kind = 'haploid';
else
    disp('');
end

% ---------------------  Regularization parameters  -----------------------
% Define parameters regularization parameters 
tau = 1;
gamma = 10;
reg_params =[tau; tau*gamma];

% -------------------------  Define variables  ----------------------------
% Define true signal f, observed signal s
s_obs = double([s_p; s_c]);
if lower(kind) == 'haploid'
    % Andrew order
    s_obs = double([s_c; s_p]); 
    f_true = double([f_h; f_n; f_p]);
    reg_params_all =[tau; tau*gamma; tau]; %[f_h; f_n; f_p]
    A(1:n,1:n) = A_c;
    A(1:n,n+1:2*n) = A_c;
    A(n+1:2*n, 2*n+1:3*n) = A_p;

    %s_obs = double([s_p; s_c]);
    %f_true = double([f_p; f_h; f_n]);
    subvectors = 3;
    N = length(f_true);
    n = N/subvectors;
    
    % Set up diagonal block matrix A (coverage matrix)
%     A = sparse(2*n, 3*n);
%     A(1:n,1:n) = A_p;
%     A(n+1:(2*n),n+1:(2*n)) = A_c;
%     A(n+1:(2*n), (2*n)+1 :(3*n)) = A_c;
    
    %reg_params_all corresponds to the penalty for 
    % each corresponding subvector in f_true
    %reg_params_all =[tau; tau; tau*gamma]; %[f_p; f_h; f_n]
    B = A;
    

elseif lower(kind) == 'diploid'
    f_true = double([z_p; z_h; z_n; y_p; y_h; y_n]);
    subvectors = 6;
    N = length(f_true);
    n = N/subvectors; 
    
    % Set up diagonal block matrix A (coverage matrix)
    A = sparse(2*n, 6*n);
    A(1:n,1:n) = A_zp;
    A(1:n, 3*n+1:4*n) = A_yp;
    A(n+1:2*n, n+1:2*n) = A_zc;
    A(n+1:2*n, 2*n+1:3*n) = A_zc;
    A(n+1:2*n, 4*n+1:5*n) = A_yc;
    A(n+1:2*n, 5*n+1:6*n) = A_yc;
    B = A;

    %reg_params_all corresponds to the penalty for [z_p; z_h; z_n; y_p; y_h; y_n]
    % in that order
    reg_params_all =[tau; tau; tau*gamma; tau; tau; tau*gamma]; 
    
%     % ALTERNATE ORDER
%     f_true = double([z_h; z_n; z_p; y_h; y_n; y_p]);
%     subvectors = 6;
%     N = length(f_true);
%     n = N/subvectors; 
%     
%     % Set up diagonal block matrix A (coverage matrix)
%     A = sparse(2*n, 6*n);
%     A(1:n,1:n)       = A_zc;
%     A(1:n, n+1:2*n)  = A_zc;
%     A(1:n,3*n+1:4*n) = A_yc;
%     A(1:n,4*n+1:5*n) = A_yc;
%     A(n+1:2*n, 2*n+1:3*n) = A_zp;
%     A(n+1:2*n, 5*n+1:6*n) = A_yp;
%     B = A;
% 
%     %reg_params_all corresponds to the penalty for [z_p; z_h; z_n; y_p; y_h; y_n]
%     % in that order
%     reg_params_all =[tau; tau*gamma; tau; tau; tau*gamma; tau]; 
end

% ------------------------- Prepare for algorithm -------------------------
% Set up function handles for computing A and A^T
AT  = @(x)A'*x;
A   = @(x)A*x;

% set maximum number of iterations, tol, and when to print to screen
miniter = 5;
maxiter = 1000;
tolerance = 1e-8;
verbose = 1;
stopcriterion = 3; 
    % 3: Relative changes in iterate
    % 4: relative changes in objective

% Simple initialization:
% initialization of f to start with AT(s) rescaled to a 
% least-squares fit to the mean intensity
f_init = (sum(sum(s_obs)).*numel(AT(s_obs)))...
    ./(sum(sum(AT(s_obs))) .*sum(sum((AT(ones(size(s_obs)))))))...
    .*AT(s_obs);





% =========================================================================
% =                  SPIRAL Method reconstruction                         =
% =========================================================================

    % Run SPIRAL algorithm
    [fhatSPIRAL, iterations_SPIRAL, objective_SPIRAL,...
    reconerror_SPIRAL, cputime_SPIRAL] ...
        = NEBULA(s_obs,A,reg_params, reg_params_all,'poisson',subvectors,...
        'maxiter',maxiter,...
        'Initialization',f_init,...
        'AT',AT,...
        'miniter',miniter,...
        'stopcriterion',stopcriterion,...
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
        'truth',f_true,...
        'verbose',verbose);

    
% =========================================================================
% =                  NEBULA Method reconstruction                         =
% =========================================================================

    % Run NEBULA algorithm
    [fhatNEBULA, iterations_NEBULA, objective_NEBULA,...
    reconerror_NEBULA, cputime_NEBULA] ...
        = NEBULA(s_obs,A,reg_params, reg_params_all,'negative binomial',subvectors,...
        'maxiter',maxiter,...
        'Initialization',f_init,...
        'AT',AT,...
        'miniter',miniter,...
        'stopcriterion',stopcriterion,...
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
        'truth',f_true,...
        'verbose',verbose);
    
    % Separate reconstruction into each signal
    if lower(kind) == 'diploid'
        fhatSPIRAL_p = 2*fhatSPIRAL(1:n)      + fhatSPIRAL(3*n+1:4*n);
        fhatSPIRAL_h = 2*fhatSPIRAL(n+1:2*n)  + fhatSPIRAL(4*n+1:5*n);
        fhatSPIRAL_n = 2*fhatSPIRAL(2*n+1:3*n)+ fhatSPIRAL(5*n+1:6*n);
        fhatSPIRAL_c = fhatSPIRAL_h + fhatSPIRAL_n;
    elseif lower(kind) == 'haploid'
        fhatSPIRAL_p = fhatSPIRAL(1:n);
        fhatSPIRAL_h = fhatSPIRAL(n+1:2*n);
        fhatSPIRAL_n = fhatSPIRAL(2*n+1:3*n);
        fhatSPIRAL_c = fhatSPIRAL_h + fhatSPIRAL_n;
    end

    %Separate reconstruction into each signal
    if lower(kind) == 'diploid'
        fhatNEBULA_p = 2*fhatNEBULA(1:n)      + fhatNEBULA(3*n+1:4*n);
        fhatNEBULA_h = 2*fhatNEBULA(n+1:2*n)  + fhatNEBULA(4*n+1:5*n);
        fhatNEBULA_n = 2*fhatNEBULA(2*n+1:3*n)+ fhatNEBULA(5*n+1:6*n);
        fhatNEBULA_c = fhatNEBULA_h + fhatNEBULA_n;
    elseif lower(kind) == 'haploid'
        fhatNEBULA_p = fhatNEBULA(1:n);
        fhatNEBULA_h = fhatNEBULA(n+1:2*n);
        fhatNEBULA_n = fhatNEBULA(2*n+1:3*n);
        fhatNEBULA_c = fhatNEBULA_h + fhatNEBULA_n;
    end      

    
    % overall plot
    figure
    subplot(1,3,1);
    [X_n,Y_n,T_n,AUC_n] = perfcurve(f_true,fhatNEBULA, 1); 
    [X_s,Y_s,T_s,AUC_s] = perfcurve(f_true,fhatSPIRAL, 1);
    plot(X_n,Y_n, '-r', 'LineWidth',2); hold on
    plot(X_s, Y_s, '--b', 'LineWidth',1.5)
    legend(strcat('NEBULA = ', num2str(AUC_n)), strcat('SPIRAL  = ', num2str(AUC_s)), 'FontSize',14, 'Location', 'southeast');
    xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
    title({'ROC Curves for total reconstruction',['\tau = ' num2str(reg_params(1)),...
        ' \gamma = ' num2str(reg_params(2))]},'FontSize',16)
    
    
    %parent reconstruction
    subplot(1,3,2);
    [X_n,Y_n,T_n,AUC_n] = perfcurve(f_p,fhatNEBULA_p, 1); 
    [X_s,Y_s,T_s,AUC_s] = perfcurve(f_p,fhatSPIRAL_p, 1);
    plot(X_n,Y_n, '-r', 'LineWidth',2); hold on
    plot(X_s, Y_s, '--b', 'LineWidth',1.5)
    legend(strcat('NEBULA = ', num2str(AUC_n)), strcat('SPIRAL  = ', num2str(AUC_s)), 'FontSize',14, 'Location', 'southeast');
    xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
    title({'ROC Curves for parent reconstruction',['\tau = ' num2str(reg_params(1)),...
        ' \gamma = ' num2str(reg_params(2))]},'FontSize',16)
    
    %child reconstruction
    subplot(1,3,3);
    [X_n,Y_n,T_n,AUC_n] = perfcurve(f_c,fhatNEBULA_c, 1); 
    [X_s,Y_s,T_s,AUC_s] = perfcurve(f_c,fhatSPIRAL_c, 1);
    plot(X_n,Y_n, '-r', 'LineWidth',2); hold on
    plot(X_s, Y_s, '--b', 'LineWidth',1.5)
    legend(strcat('NEBULA = ', num2str(AUC_n)), strcat('SPIRAL  = ', num2str(AUC_s)), 'FontSize',14, 'Location', 'southeast');
    xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16);
    title({'ROC Curves for child reconstruction',['\tau = ' num2str(reg_params(1)),...
        ' \gamma = ' num2str(reg_params(2))]},'FontSize',16)
% Confusion Matrix
% Not well labelled ATM, so may be hard to understand
    figure
    subplot(1,2,1)
    threshold = 0.5;
    fhatNEBULA_bin = fhatNEBULA > threshold;
    C = confusionmat(f_true, double(fhatNEBULA_bin), 'Order', [1, 0]);
    cm = confusionchart(C);
    cm.Title = 'Structural Variant Confusion Matrix using NEBULA';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    cm.NormalizedValues;
    subplot(1,2,2)
    fhatSPIRAL_bin = fhatSPIRAL > threshold;
    C = confusionmat(f_true, double(fhatSPIRAL_bin), 'Order', [1, 0]);
    cm = confusionchart(C);
    cm.Title = 'Structural Variant Confusion Matrix using SPIRAL';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    cm.NormalizedValues;

%end