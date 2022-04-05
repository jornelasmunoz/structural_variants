% Gen1P1CData_nov_NegBinom_Data.m

% This code generates a parent and a child data and saves them to
% appropriate file.  Parameters to change are n, k, lambda_p, lambda_c,
% etc., but do not edit naming of files.

% comments which start with APL were written by Andrew Peter Lazar
clear all
clc
% APL: dispersion parameter r for negative binomial for the maximum value
% we should have r=1
r=1;

% create parent data

% Initialize the size and number of nonzero elements in parent signals, f_p
n = 3; %10^2;
k = 1; %10; %number of SVs

% Initialize lambda and epsilon values

lambda_c = 4;
lambda_p = 4;
erreps = 0.01;

% Change depending on appropriate loading data
pctNovel = 1;%0.2; % if pctNovel = 1, then child signal does not share any variants with the parent



q    = randperm(n);  

%Give parent the first k variants
f_p_neg_binom = zeros(n,1);
f_p_neg_binom(q(1:k)) = 1;


%%%%%%%%%%%%%%%%%%%%% Create Child Signal %%%%%%%%%%%%%%%%%%%%%%

%Child has 2 types of variants:
%(1) inherited f_c_inh
%(2) Novel (not interited from parent)
%Give child the variant from k*pct, k*pct+k
startVal          = k*(pctNovel); %50
endVal            = startVal + k; %550
f_c_neg_binom = zeros(n,1);
f_c_neg_binom(q(startVal+1:endVal)) = 1; %this gives the child SVs where pctNovel*k of
                                    %them overlap with the parent and the
                                    %other (k+1):endVal are novel

f_c_inh_neg_binom = zeros(n,1);
f_c_nov_neg_binom = zeros(n,1);
f_c_inh_neg_binom(q(startVal+1:k)) = 1;  % the variants inherited from the parent
f_c_nov_neg_binom(q(k+1:endVal)) = 1;  % the variants that are novel in the child



% Create observations based on y_c ~ Poisson ( lambda_c - eps) I*fc + epsi))
% Generate observations and A_c, A_p matrices

A_c_neg_binom = (lambda_c - erreps)*speye(n);

%APL: computing the mean and variance vectors for the child
%mu_c= A_c*(f_c_nov + f_c_inh) + erreps;
mu_c= A_c_neg_binom*(f_c_nov_neg_binom + f_c_inh_neg_binom) + erreps;
var_c= mu_c + (1/r)*(mu_c).^2;


%APL: Computing the negative binomial random variable where the parameters
%represent nbinrnd(number of successes, probability of success) for the
%child

y_c_neg_binom =nbinrnd(mu_c.^2./(var_c - mu_c), mu_c./var_c);

%y_c_notp1 = poissrnd((lambda_c - erreps)*f_c_notp1 + erreps);
%y_c_notp2 = poissrnd((lambda_c - erreps)*f_c_notp2 + erreps);

A_p_neg_binom = (lambda_p - erreps)*speye(n);

%APL: computing the mean and variance vectors for the parent
%mu_p=A_p*f_p + erreps;
mu_p =A_p_neg_binom*f_p_neg_binom + erreps;
var_p = mu_p +(1/r)*(mu_p).^2;

%APL: Computing the negative binomial random variable where the parameters
%represent nbinrnd(number of successes, probability of success) for the
%parent

y_p_neg_binom = nbinrnd(mu_p.^2./(var_p - mu_p), mu_p./var_p);

%y_p = poissrnd((lambda_p - erreps)*f_p + erreps);
%A_p2 = (lambda_p2 - erreps)*speye(n);
%y_p2 = poissrnd((lambda_p2 - erreps)*f_p2 + erreps);


% Clear unnecessary variables
clear temp;

% save data in format p1(lambda_p1)_p2(lambda_p2)_c(lambda_c)_percentsharedpar.mat
datname = sprintf('dummy_data');%sprintf('neg_binom_nov_p%d_c%d_%dperNov',lambda_p,lambda_c,pctNovel*100);

save(datname)