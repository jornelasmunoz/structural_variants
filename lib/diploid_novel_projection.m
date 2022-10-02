function [f_feas] = diploid_novel_projection(f, subvectors, maxit)
%disp("Entering proj")
%Inputs:
%   f: subsolution reconstruction vector 
%   subvectors: number of vectors which f is composed of
%   maxit: number of maximum iterations to project the vector onto the
%          feasible region
N = length(f);
n = N/subvectors;

% The following assumes the vector f = [zP, zH, zN, yP, yH, yN]
zP = f(    1:  n);
zH = f(  n+1:2*n);
zN = f(2*n+1:3*n);
yP = f(3*n+1:4*n);
yH = f(4*n+1:5*n);
yN = f(5*n+1:6*n);

% proj_mat =  [1, 0, 0, 1, 0, 0;
%             -1, 0, 0,-1, 0, 0;
%              0, 1, 0, 0, 1, 0;
%              0,-1, 0, 0,-1, 0;
%              0, 0, 1, 0, 0, 1;
%              0, 0,-1, 0, 0,-1;
%              0, 1, 1, 0, 1, 1;
%              0,-1,-1, 0,-1,-1;
%              0, 1, 0, 0, 0, 0;
%              1,-1, 0, 1, 0, 0;
%              1, 0, 0, 0, 0, 0;
%             -1, 1, 0, 0, 1, 0;
%             -1, 0,-1,-1, 0,-1];
% b = [0;-1;0;-1;0;-1;0;-1;0;0;0;0;-1];
% 
% for i = 1:n
%     x_hat = [zP(i);zH(i);zN(i);yP(i);yH(i);yN(i)];
%     x_hat = Basic_GP(proj_mat, b, x_hat);
%     zP(i) = x_hat(1);
%     yP(i) = x_hat(2);
%     zH(i) = x_hat(3);
%     yH(i) = x_hat(4);
%     zN(i) = x_hat(5);
%     yN(i) = x_hat(6);
% end

% Check feasibility of solution and initialize projections, if needed
for i = 1:n
    % start iteration counter
    iter = 1;
    
    while (iter < maxit)

            % Check if reconstruction vector element already satisfies the
            % constraints
            
            if ( (0 <= zP(i)) && (zP(i) <= 1) && ...
                 (0 <= zH(i)) && (zH(i) <= 1) && ...
                 (0 <= zN(i)) && (zN(i) <= 1) && ...
                 (0 <= yP(i)) && (yP(i) <= 1) && ...
                 (0 <= yH(i)) && (yH(i) <= 1) && ...
                 (0 <= yN(i)) && (yN(i) <= 1) && ...
                 (0 <= zP(i) + yP(i)) && (zP(i) + yP(i) <= 1) && ...
                 (0 <= zH(i) + yH(i)) && (zH(i) + yH(i) <= 1) && ...
                 (0 <= zN(i) + yN(i)) && (zN(i) + yN(i) <= 1) && ...
                 (0 <= zH(i) + yH(i) + zN(i) + yN(i)) && ...
                 (zH(i) + yH(i) + zN(i) + yN(i) <= 1) && ...
                 (zH(i) <= (zP(i) + yP(i))) && ...
                 (zP(i) <= (zH(i) + yH(i))) && ...
                 (zN(i) + yN(i) <= 1 - (zP(i) + yP(i))) && ...
                 (1 - (zP(i) + yP(i)) <= 1) )
             
            break;
            
            else
            % If constraints are not satisfied, start projections onto
            % feasible regions
            
            % Step 0: Initialize indicator variables for block coord descent
            zP_pos = sort([0, zP(i),1]); z_p = zP_pos(2);
            zH_pos = sort([0, zH(i),1]); z_h = zH_pos(2);
            zN_pos = sort([0, zN(i),1]); z_n = zN_pos(2);
            yP_pos = sort([0, yP(i),1]); y_p = yP_pos(2);
            yH_pos = sort([0, yH(i),1]); y_h = yH_pos(2);
            yN_pos = sort([0, yN(i),1]); y_n = yN_pos(2);
             
            %z_p = zP(i); y_p = yP(i);
            
            % Enforce constraint: homog. + heter. SV <= 1
                if (z_h + y_h > 1)
                    z_h = 0.5;
                    y_h = 0.5;
                end
                if (z_n + y_n > 1)
                    z_n = 0.5;
                    y_n = 0.5;
                end
                
            end
           
            
            % Continue projections with previous iterate
        %===================================================
        %=                 PROJECTIONS                     =
        %===================================================
        % Step 1: Fix P,N, and vary H 
        [zH_new, yH_new] = fix_one_indiv(z_p,y_p,z_h,y_h,z_n,y_n, 'H');
        z_h = zH_new;
        y_h = yH_new;
        
        % Step 2: Fix P,H, and vary N
        [zN_new, yN_new] = fix_one_indiv(z_p,y_p,z_h,y_h,z_n,y_n, 'N');
        z_n = zN_new;
        y_n = yN_new;
        
        % Step 3: Fix H, N, and vary P
        [zP_new, yP_new] = fix_one_indiv(z_p,y_p,z_h,y_h,z_n,y_n, 'P'); 
        z_p = zP_new;
        y_p = yP_new;
        
        % re-assign values to vector
        zP(i) = z_p;
        yP(i) = y_p;
        zH(i) = z_h;
        yH(i) = y_h;
        zN(i) = z_n;
        yN(i) = y_n;

        iter = iter + 1;        
    end
end
        
        
        
   f_feas = [zP;zH;zN;yP;yH;yN];  

end

