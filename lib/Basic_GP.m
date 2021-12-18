function [x_final] = Basic_GP(A, b, x_hat)

%%% load in saved variables %%%
% A = A.B;
% b = b.b;

% primal objective function
f = @(i) (1/2)*(norm(i-x_hat,2))^2;

% dual objective function
r = A*x_hat - b;
g = @(j) (1/2)*(j'*(A*A')*j)+(j'*r);

%%%%%%% Initaialize variables %%%%%%%
lambda = zeros(length(b),1);
k = 1;
%eps1 = 1e-5;
%eps2 = 1e-3;
eps3 = 1e-5;
x_curr = x_hat;
beta = 1/2;
mu = 0.1;

% calculate value of dual obj function
dual_obj = g(lambda);

% calculate value of primal obj function
primal_obj = f(x_curr);

% duality gap
dual_gap = primal_obj - (-1 * dual_obj);

% calc residual
r = A*x_hat - b;

v = sum(r(r < 0));
% check feasibiity of primal problem
if (any(r < 0))
    infeasibility(k) = v;
else
    infeasibility(k) = 0;
end

% fprintf('===========================================================================================\n')
% fprintf('Itn   Primal Obj.   Dual Obj.   Dual. Gap   Infeasibilty     alpha     ||x_n - x_{n-1}}||_2\n')
% fprintf('===========================================================================================\n')
% fprintf('%3d   %10.6f   %10.6f   %10.6f   %10.6f     N/A          N/A\n', k, primal_obj, -dual_obj, dual_gap, infeasibility);
% tt = 2;

k = 2;
tic
%%%%%%% Proposed method %%%%%%%
while ( (k == 2 || norm(x_curr-x_prev,2) > eps3) && k < 100)
    % norm(x_curr-x_prev,2)
    % abs(dual_gap)
    % abs(infeasibility)
    
    % set current value of x as previous 
    x_prev = x_curr;
    
    % calc gradient
    gr = (A*(A'*lambda) + r);
    
    % check if gradient = 0
    if (norm(gr,2) == 0)
        fprintf('Found minimizer; gradient of g = 0 \n')
        break
    end
    
    % intialize alpha
    alpha = norm(gr, 2)^2 / norm(A'*gr, 2)^2;
    
    % backtracking line search
    while true
        temp = lambda - alpha * gr;
        temp(temp < 0) = 0;
        
        condition1 = g(temp);
        condition2 = g(lambda) - mu * gr'* (lambda - temp);
        
        if (condition1 <= condition2)
            break
        end
        
        alpha  = beta * alpha;
    end
    
    % set new lambda using alpha from above
    lambda = temp;

    % set negative components to zero
    lambda(lambda < 0) = 0;
    
    % save value of dual obj function
    dual_obj = g(lambda);
    
    % update solution x to primal problem
    x_curr = x_hat + A'*lambda;
    
    % save value of primal obj function
    primal_obj = f(x_curr);
    
    % calc duality gap
    dual_gap = primal_obj - (-1 * dual_obj);
    
    res = A*x_curr - b;
    v = sum(res(res < 0));
    % check feasibiity of primal problem
    if (any(res < 0))
        infeasibility = v;
    else
        infeasibility = 0;
    end
    
    
%     if tt == 10
%         fprintf('===========================================================================================\n')
%         fprintf('Itn   Primal Obj.   Dual Obj.   Dual. Gap   Infeasibilty     alpha     ||x_n - x_{n-1}}||_2\n')
%         fprintf('===========================================================================================\n')
%         tt = 0; 
%     end
%     fprintf('%3d   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f\n', k, primal_obj, -dual_obj, dual_gap, infeasibility, alpha, norm(x_curr-x_prev,2));
%     tt = tt+1;
    
    k = k+1; % update counter
    
end

elapsed_time = toc;
x_final = x_curr;
total_iterations = k-2;


%%%%%% Save Variables %%%%%% 
% filename = [pwd '/experiments/Basic_GP_m_' num2str(m) '_n_' num2str(n) '_run_' num2str(run) '.mat'];
% save(filename, 'dual_obj', 'primal_obj', 'infeasibility', 'x', 'lambda', 'elapsed_time');

%%%%% Print Statements %%%%%%
% fprintf('------------------------------------ \n')
% fprintf('Projection Results \n')
% fprintf('Elapsed time is %f seconds. \n', elapsed_time)
% fprintf('Total number of iterations: %d\n', k-1)
% fprintf('------------------------------------ \n')

end
