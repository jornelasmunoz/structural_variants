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

% Load Simulated Data 
addpath([genpath('/Users/jocelynornelasmunoz/Desktop/Research/structural_variants/'), ...
         genpath('/Users/jocelynornelas/iCloud Drive (Archive)/Desktop/UC Merced/Research/structural_variants/') ])

filename = 'data/haploid_20pctNovel_10k_100n.mat'; load(filename)

% f_p = f_p_neg_binom;
% f_h = f_c_inh_neg_binom;
% f_n = f_c_nov_neg_binom;
% s_p = y_p_neg_binom;
% s_c = y_c_neg_binom;
% A_c = A_c_neg_binom;
% A_p = A_p_neg_binom;

% Define true signal f, observed signal s
f_true = double([f_p; f_h; f_n]);
s_obs = double([s_p; s_c]);
subvectors = 3;
N = length(f_true);
n = N/subvectors; 

% Set up diagonal block matrix A (coverage matrix)
A = sparse(2*n, 3*n);
A(1:n,1:n) = A_p;
A(n+1:2*n,n+1:2*n) = A_c;
A(n+1:2*n, 2*n+1:3*n) = A_c;
B = A;

% Set up function handles for computing A and A^T
AT  = @(x)A'*x;
A   = @(x)A*x;

% set maximum number of iterations, tol, and when to print to screen
maxiter = 1000;
tolerance = 1e-8;
verbose = 10;


% Simple initialization:
% initialization of f to start with
% AT(s) rescaled to a least-squares fit to the mean intensity
f_init = (sum(sum(s_obs)).*numel(AT(s_obs)))...
    ./(sum(sum(AT(s_obs))) .*sum(sum((AT(ones(size(s_obs)))))))...
    .*AT(s_obs);

% Define parameters tau and gamma 
tauvals= [1];
gamma= 10;

for i=1:length(tauvals)
     t= tauvals(i);
     tau=[t; t*gamma];
    
     

    % =====================================================================
    % =                  SPIRAL Method reconstruction                     =
    % =====================================================================

    % Run SPIRAL algorithm
    [fhat_SPIRAL, iterations_SPIRAL, objective_SPIRAL,...
    reconerror_SPIRAL, cputime_SPIRAL] ...
        = NEBULA(s_obs,A,tau,'poisson',subvectors,...
        'maxiter',maxiter,...
        'Initialization',f_init,...
        'AT',AT,...
        'miniter',5,...
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
        'truth',f_true,...
        'verbose',verbose);

%     % Separate reconstruction into each signal
%     fhat_SPIRAL_p = 2*fhat_SPIRAL(1:n)      + fhat_SPIRAL(3*n+1:4*n);
%     fhat_SPIRAL_h = 2*fhat_SPIRAL(n+1:2*n)  + fhat_SPIRAL(4*n+1:5*n);
%     fhat_SPIRAL_n = 2*fhat_SPIRAL(2*n+1:3*n)+ fhat_SPIRAL(5*n+1:6*n);
%     fhat_SPIRAL_c = fhat_SPIRAL_h + fhat_SPIRAL_n;

    [tpr_SPIRAL,fpr_SPIRAL,thresholds_SPIRAL] = roc(f_true,fhat_SPIRAL);
    % =====================================================================
    % =                  NEBULA Method reconstruction                     =
    % =====================================================================

    % Run NEBULA algorithm
    [fhat_NEBULA, iterations_NEBULA, objective_NEBULA,...
    reconerror_NEBULA, cputime_NEBULA] ...
        = NEBULA(s_obs,A,tau,'negative binomial',subvectors,...
        'maxiter',maxiter,...
        'Initialization',f_init,...
        'AT',AT,...
        'miniter',5,...
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
        'truth',f_true,...
        'verbose',verbose);
    
%     % Separate reconstruction into each signal
%     fhat_NEBULA_p = 2*fhat_NEBULA(1:n)      + fhat_NEBULA(3*n+1:4*n);
%     fhat_NEBULA_h = 2*fhat_NEBULA(n+1:2*n)  + fhat_NEBULA(4*n+1:5*n);
%     fhat_NEBULA_n = 2*fhat_NEBULA(2*n+1:3*n)+ fhat_NEBULA(5*n+1:6*n);
%     fhat_NEBULA_c = fhat_NEBULA_h + fhat_NEBULA_n;
%     
    [X_n,Y_n,T_n,AUC_n] = perfcurve(f_true,fhat_NEBULA, 1); 
    [X_s,Y_s,T_s,AUC_s] = perfcurve(f_true,fhat_SPIRAL, 1);
    plot(X_n,Y_n, '-r', 'LineWidth',2); hold on
    plot(X_s, Y_s, '--b', 'LineWidth',1.5)
    legend(strcat('NEBULA = ', num2str(AUC_n)), strcat('SPIRAL = ', num2str(AUC_s)))
    %ROC_curve
    %save_to_JSON
    figure
    fhat_NEBULA_bin = fhat_NEBULA > 0.5;
    C = confusionmat(f_true, double(fhat_NEBULA_bin), 'Order', [1, 0]);
    cm = confusionchart(C);
    cm.Title = 'Structural Variant Confusion Matrix using NEBULA';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    cm.NormalizedValues;
    figure
    fhat_SPIRAL_bin = fhat_SPIRAL > 0.5;
    C = confusionmat(f_true, double(fhat_SPIRAL_bin), 'Order', [1, 0]);
    cm = confusionchart(C);
    cm.Title = 'Structural Variant Confusion Matrix using SPIRAL';
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
    cm.NormalizedValues;

    %calculate precision and recall
    % look into why ROC curve is not good for imbalanced data
end