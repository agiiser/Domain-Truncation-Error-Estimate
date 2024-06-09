clear all
clc
close all
l=999;          % Number of grid points in each axis
M=zeros(l,l);   % Intializing M
N=M; O=M;       % Intializing N, and O
ds=20/(l+1);    % step size for stock price axis
s=ds*(1:1:l);   % Stock price axis grid values
sm=max(s)+ds;   % Boundary of the stock price axis
T=1;            % Terminal time
dt=T/(l+1);     % step size for time axis
t=dt*(1:1:l);   % time axis grid values
for j =1:length(s)
    for i=1:length(t)
       [a,b,c,D]=est(t(i),s(j), sm, T);
       M(i,j)=exp(a);   % Our estimate
       N(i,j)=exp(b);   % Kangro estimate
       O(i,j)=c;        % Indicator function of the sub-domain
    end
end
%Figure-1
contour(s,t,N-M,'ShowText','on', LineWidth=2)
axis([0 sm 0 T])
grid on
xlabel('stock price'), ylabel('time')
 %Figure-2 surface plot uses the column and row indices of the elements in Z as the x- and y-coordinates.
 surf(s,t,M,'FaceColor','g','EdgeColor','none') 
 hold on
 surf(s,t,N,'FaceColor','flat','EdgeColor','none')
 surf(s,t, O,'FaceColor','r','EdgeColor','none')
 text(1,1,1.2,"Indicator function of "+"𝒟")
 text(18,0.9,1.4,"\Psi_1")
 text(18,0.1,0.6,"\Psi_1"), text(18.2,0.1,0.7,"-")
 grid on
 axis([0 sm 0 T 0 2])
 xlabel('stock price'), ylabel('time'), zlabel('Estimate value')
 alpha 0.9

function [f1, f2, f3, D]=est(t1,s1, sm, T)
sig=0.1;  % 0.40; %  
r  =0.2;  % 0.01; %  
D=sig^2-2*r; Dp=max(D,0); am=sig^2;
f1= (-log(sm/s1)*(2+ (Dp/am)*log(sm/s1))+ (am+abs(D))*(T-t1))/(2*Dp*(T-t1)+ am/(am+Dp));  %Our paper
f2=(-log(sm/s1)*(log(sm/s1)+min(0,D)*(T-t1))/(2*am*(T-t1))); %Kangro Theorem 4 equation (22)
f3=1;
    if log(sm/s1)+D*(T-t1)<0
        f3=0;
    end
end