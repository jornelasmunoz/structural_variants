function [f_feas] = diploid_novel_projection(f, subvectors, maxit)
%Inputs:
%   f: subsolution reconstruction vector 
%   subvectors: number of vectors which f is composed of
%   maxit: number of maximum iterations to project the vector onto the
%   feasible region
N = length(f);
n = N/subvectors;

% The following assumes the vector f = [zP, zH, zN, yP, yH, yN]
zP = f(    1:  n);
zH = f(  n+1:2*n);
zN = f(2*n+1:3*n);
yP = f(3*n+1:4*n);
yH = f(4*n+1:5*n);
yN = f(5*n+1:6*n);


% Check feasibility of solution and initialize projections, if needed

for i = 1:n
    
    % start iteration counter
    iter = 1;
    
    while (iter < maxit)
        if (iter == 1)
            %disp('iter ==1')
            % Check if reconstruction vector already satisfies the
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
             iter = maxit;
             f_feas = [zP;zH;zN;yP;yH;yN];
             %disp('constraints feasible initially')
            return;
            
            else
            % If constraints are not satisfied, start projections onto
            % feasible regions
            %disp('constraints NOT feasible initially')
            
            % Ensure elements are nonnegative
            %zP_pos = sort([0, zP(i),1]); z_p = zP_pos(2);
            zH_pos = sort([0, zH(i),1]); z_h = zH_pos(2);
            zN_pos = sort([0, zN(i),1]); z_n = zN_pos(2);
            %yP_pos = sort([0, yP(i),1]); y_p = yP_pos(2);
            yH_pos = sort([0, yH(i),1]); y_h = yH_pos(2);
            yN_pos = sort([0, yN(i),1]); y_n = yN_pos(2);
            
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
           
            
        else 
            % Continue projections with previous iterate
            z_h = zH_new;
            y_h = yH_new;
            z_n = zN_new;
            y_n = yN_new;
        end
        
        %===================================================
        %=                 PROJECTIONS                     =
        %===================================================
        
        % Fix H, N, and vary P
        [zP_new, yP_new] = fix_one_indiv(zP(i),yP(i),z_h,y_h,z_n,y_n, 'P'); 
        z_p = zP_new;
        y_p = yP_new;
        % Fix P,N, and vary H 
        [zH_new, yH_new] = fix_one_indiv(z_p,y_p,zH(i),yH(i),z_n,y_n, 'H');
        z_h = zH_new;
        y_h = yH_new;
        
        % Fix P,H, and vary N
        [zN_new, yN_new] = fix_one_indiv(z_p,y_p,z_h,y_h,zN(i),yN(i), 'N');
        %z_n = zN_new;
        %y_n = yN_new;
        
        iter = iter + 1;
    end
    zP(i) = zP_new;
    yP(i) = yP_new;
    zH(i) = zH_new;
    yH(i) = yH_new;
    zN(i) = zN_new;
    yN(i) = yN_new;
end

f_feas = [zP;zH;zN;yP;yH;yN];