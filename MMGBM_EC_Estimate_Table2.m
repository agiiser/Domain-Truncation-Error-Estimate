%%%%%%%%%%%%%%%%%%   START for Table 2 of our paper
clc;clear;
format long
K=10;smax=3*K;r1=0.15;r2=0.15;sig1=0.14;sig2=0.20;T=3;t=0;
A=[(sig1^2)-(2*r1),(sig2^2)-(2*r2)];dl=min(A);
B=[sig1,sig2];sig=max(B);Dlp=max(dl,0);
s=[7,8,9,10,11,12,13,14,15]; %%% ||v-w||_{\infty}<K for both works
for i=1:length(s)
lhs=log(smax/s(i)); rhs=-(T-t)*dl;
fprintf('The error estimate for s= %d.\n',s(i));
if lhs<rhs
   fprintf('\n Kangro paper provides no error estimate for s= %d.\n',s(i));
else
%%% Estimate in Kangro’s paper
 err_estimate_in_kangro_work = K*exp(-(log(smax/s(i))*(log(smax/s(i))+min(0,dl)*(T-t)))/(2*(sig^2)*(T-t)))
end
%%% Derived estimate
derived_err_estimate= K*exp(-(log(smax/s(i))*((Dlp/(sig^2))*(log(smax/s(i)))+2)+(((sig^2)+abs(dl))*(T-t)))/(2*((Dlp*(T-t))+((sig^2)/((sig^2)+Dlp)))))
end
%%%%%%%%%%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%%%%%%%%%