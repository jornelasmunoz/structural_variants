function [z_proj,y_proj] = fix_one_indiv(zP,yP,zH,yH,zN,yN,vary)
%fix_one_indiv takes in the elements of a vector in 6 dimensions and
%returns the projection onto a 2d plane (z_proj, y_proj)
% Inputs:
%   zP,yP,zH,yH,zN,yN
%   vary: string in ['P', 'H', 'N'] which determines the plane we project
%         onto


% we check whether the point is within the feasible region. If it is not,
% we return a point within the feasible region closest to our input point


% initialize projection variables
z_proj = 0;
y_proj = 0;

switch upper(vary)
    case {"P","H"}
    % vary P, fix H and N 
    % (zP,yP) is a point in the plane
    % zH,yH,zN,yN are fixed values (i.e. constants)
    % Similarly if we vary H and fix P,N
    
    if (vary == 'P')
        x1 = zP; x2 = yP;
        zF = zH; yF = yH;
        %disp('vary P')
    elseif (vary == 'H')
        x1 = zH; x2 = yH;
        zF = zP; yF = yP;
        %disp('vary H')
    end
    
    % Feasible Region
    if (x1 >= 0 && x2 >= 0 && x2 >= zF -x1 && ...
            x2 <= 1 - (zN + yN) - x1 && ...
            x1 <= (zF + yF) )
        z_proj = x1;
        y_proj = x2;
        %disp('Feasible Region')
        return;
        
    % Region 1 
    elseif (x2 >= 1 - (zN + yN) - x1 && ...
            x2 <= x1 + 1 - (zN + yN) && ...
            x2 >= x1 + 1 - 2*(zF + yF)- (zN + yN))
        %disp('Region 1')
        z_proj = (1/2)*(1 + x1 - x2 - (zN + yN)) ;
        y_proj = 1 - (zN + yN) - z_proj;
        return;
    
    % Region 2
    elseif (x1 >= (zF + yF) && ...
            x2 >= 1 - (zN + yN) - (zF + yF) && ...
            x2 <= x1 + 1 - 2*(zF + yF) - (zN + yN))
        %disp('Region 2')
        z_proj = zF + yF;
        y_proj = 1 - (zN + yN) - z_proj;
        return;
    
    % Region 3
    elseif (x1 >= (zF + yF) && x2 >= 0 && ...
            x2 <= 1 - (zN + yN) - (zF + yF))
        %disp('Region 3')
        z_proj = zF + yF;
        y_proj = x2;
        return;
    
    % Region 4
    elseif (x1 >= (zF + yF) && x2 <= 0)
        %disp('Region 4')
        z_proj = zF + yF;
        y_proj = 0;
        return;
   
    % Region 5    
    elseif (x1 >= zF && x2 <= 0 && x1 <= (zF + yF))
        %disp('Region 5')
        z_proj = x1;
        y_proj = 0;
        return;
        
    % Region 6
    elseif (x2 <= x1 - zF && x1 <= zF && x2 <= 0)
        %disp('Region 6')
        z_proj = zF;
        y_proj = 0;
        return;
        
    % Region 7 
    elseif (x2 <= x1 + zF && ...
            x2 <= zF - x1 && ...
            x2 >= x1 - zF)
        %disp('Region 7')
        z_proj = (1/2) * (zF + x1 - x2);
        y_proj = (1/2) * (zF + x2 - x1);
        return;
        
    % Region 8
    elseif (x1 <= 0 && x2 <= zF && x2 >= x1 + zF)
        %disp('Region 8')
        z_proj = 0;
        y_proj = zF;
        return;
        
    % Region 9
    elseif (x1 <= 0 && x2 >= zF && x2 <= 1- (zN + yN))
       % disp('Region 9')
        z_proj = 0;
        y_proj = x2;
        return;
        
    % Region 10
    elseif (x2 >= x1 + 1 && x2 >= 1)
        %disp('Region 10')
        z_proj = 0;
        y_proj = 1- (zN + yN);
        return;
    %return;
    end
    
    
    
    
    case "N"
    % vary N, fix P and H 
    % (zN,yN) is a point in the plane
    % zP,yP,zH,yH are fixed values (i.e. constants)
    %disp('vary N')
    min_of_P_and_H = min([1- (zP+yP), 1- (zH+yH)]);
    
    % Feasible Region
    if (zN >= 0 && yN >= 0 && yN <= min_of_P_and_H -zN)
        z_proj = zN;
        y_proj = yN;
        return;
    
    % Region 1
    elseif (yN >= min_of_P_and_H - zN  && ...
            yN >= zN - min_of_P_and_H  && ...
            yN <= zN + min_of_P_and_H )
        z_proj = 1/2 * (zN - yN + min_of_P_and_H);
        y_proj = min_of_P_and_H - z_proj;
        return;
    
    % Region 2
    elseif (yN <= zN - min_of_P_and_H && ...
            zN >= min_of_P_and_H)
        z_proj = min_of_P_and_H;
        y_proj = 0;
        return;
    
    % Region 3
    elseif (yN <= 0 && zN >= 0 && ...
            zN <= min_of_P_and_H)
        z_proj = zN;
        y_proj = 0;
        return;
    
    % Region 4
    elseif (yN <= 0 && zN <= 0)
        z_proj = 0;
        y_proj = 0;
        return;
    
    
    % Region 5
    elseif (yN >= 0 && zN <= 0 && ...
            yN >= min_of_P_and_H)
        z_proj = 0;
        y_proj = yN;
        return;
    
    % Region 6
    elseif (yN >= zN + min_of_P_and_H && ...
            yN >= min_of_P_and_H)
        z_proj = 0;
        y_proj = min_of_P_and_H;
        return;
    end
    
%     otherwise
%         disp('vary must be the string P, H, or N')
end 
end