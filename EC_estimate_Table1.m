%%%%%%%%%%%%%%%%%%%   START for Table 1 of our paper
clc
clear all
format long
K=10;smax=2*K;r=0.12;sig=0.1;t=0;T=2;dl=((sig^2)-2*r);Dlp=max(dl,0);
s=[7,8,9,10,11,12,13,14,15];
for i=1:length(s)
lhs=log(smax/s(i)); rhs=-(T-t)*((sig^2)-2*r);%% ||v-w||_{\infty}<K for both works
fprintf('The error estimate for s= %d.\n',s(i));
if lhs<rhs
  fprintf('\n Kangro paper provides no error estimate for s= %d.\n',s(i));
else
%%% error by Kangro
err_estimate_in_kangro_work = K*exp(-(log(smax/s(i))*(log(smax/s(i))+min(0,dl)*(T-t)))/(2*(sig^2)*(T-t)))
end
%%% Derived estimate
derived_err_estimate= K*exp(-(log(smax/s(i))*((Dlp/(sig^2))*(log(smax/s(i)))+2)+(((sig^2)+abs(dl))*(T-t)))/(2*((Dlp*(T-t))+((sig^2)/((sig^2)+Dlp)))))
end