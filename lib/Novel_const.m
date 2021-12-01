function [f_subs,region]= Novel_const(f_vec)

% The following function enforces constraints of our problem within 
% the feasible region for the minimization subproblem of the form:
%
%  phi(f_c, f_n, f_p) = 1/2*((f_c - a)^2 + (f_n - b)^2 + (f_p - c)^2)
%
% where a= s_c - mu, b= s_n - mu_n, and c = s_p - mu

N= length(f_vec);
m=N/3;
c_vec=f_vec(1:m);
n_vec=f_vec(m+1:2*m);
p_vec=f_vec((2*m+1):3*m);


%Initialize vectors r, s, and t

r= zeros(m,1);
s= zeros(m,1);
t= zeros(m,1);
region=zeros(m,1);


% This ....
for i=1:m
    c=c_vec(i);
    n=n_vec(i);
    p=p_vec(i);


%--------------------------------------------------
% Inside Region
%--------------------------------------------------

    if (c-p)<= 0 && c >=0 && n >= 0 && (n+p-1)<= 0

            r(i)=c;
            s(i)=n;
            t(i)=p;
            region(i)=1;
        
%--------------------------------------------------
% Vertices
%--------------------------------------------------


% Vertex (1) --> (0,0,0)

    elseif c<=(-p) && n < 0 && p<= 0
        r(i)=0;
        s(i)=0;
        t(i)=0;
        region(i)=2;
   
% Vertex (2) --> (0,0,1)

    elseif c<0 && n<(p-1) && p>= 1
        r(i)=0;
        s(i)=0;
        t(i)=1;
        region(i)=3;
 
% Vertex (3) --> (0,1,0)

    elseif c<(n-p-1) && n>= 1 && p<(n-1)
        r(i)=0;
        s(i)=1;
        t(i)=0;
        region(i)=4;
        
% Vertex (4) --> (1,0,1)

    elseif c>= 1 && n<(c+p-2) && p>2-c
        r(i)=1;
        s(i)=0;
        t(i)=1;
        region(i)=5;
%-----------------------------------------
% Edge
%-----------------------------------------

%Edge R_(0,n,0)

    elseif c<-p && 0<n && n<1 && p<0
        r(i)=0;
        s(i)=n;
        t(i)=0;
        region(i)=6;
        
%Edge R_(r_1,s_1,t_1)

    elseif c> 2-(2*n)-p && n<1+c+p && p<(2*c)-1+n && p<2-c+n
        r(i)=(1+c-n+p)/3;
        s(i)=(2-c+n-p)/3;
        t(i)=(1+c-n+p)/3;
        region(i)=7;
        
%Edge R_(0,0,p)

    elseif c<0 && n<0 && 0<p && p<1
        r(i)=0;
        s(i)=0;
        t(i)=p;
        region(i)=8;
        
%Edge R_(0,s_2,t_2)

    elseif c<0 && n<1+p && 1-n<p && p<n+1
        r(i)=0;
        s(i)=(1+n-p)/2;
        t(i)=(1-n+p)/2;
        region(i)=9;

%Edge R_(r_3,0,t_3)

    elseif p<c && c<2-p && n<0 && p>(-c)
        r(i)=(c+p)/2;
        s(i)=0;
        t(i)=(c+p)/2;
        region(i)=10;

%Edge R_(c,0,1)

    elseif 0<c && c<1 && n<p-1 && p>1
        r(i)=c;
        s(i)=0;
        t(i)=1;
        region(i)=11;

%----------------------------------------------
% Surfaces
%----------------------------------------------


%Face R_(c,0,p)

    elseif c>=0 && c<=p && n<=0 && p<=1
        r(i)=c;
        s(i)=0;
        t(i)=p;
        region(i)=12;

%Face R_(0,n,p)

    elseif c<=0 && n>=0 && p<=1-n && p>=0
        r(i)=0;
        s(i)=n;
        t(i)=p;
        region(i)=13;

%Face R_(c,s_4,t_4)

    elseif c>=0 && n<=1-(2*c)+p && 1-n<=p && p<=n+1
        r(i)=c;
        s(i)=(1-p+n)/2;
        t(i)=(1+p-n)/2;
        region(i)=14; 
        
%Face R_(r_5,n,t_5)

    elseif c>=p && c>=-p && 0<=n && n<=1 && p<=-c-(2*n)+2
        r(i)=(c+p)/2;
        s(i)=n;
        t(i)=(c+p)/2;
        region(i)=15;
   
    else
        region(i)=-1;
   
    end

end
	f_subs=[r;s;t];

end


    
    

