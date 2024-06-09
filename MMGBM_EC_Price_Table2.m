%Antithetic Simulation Method for Pricing EC under MMGBM
clc
clear all
%_Parameters and Variables____________________
TPM=[0,1;1,0];  %Transition Probability Matrix of arbitrary number of states
CPM=transpose(cumsum(transpose(TPM)));
L=[0.5 0.5];    %input('\n Enter the rate vector:')
sigma=[0.14, 0.20];
r=0.15;
K=10;           %Strike price
Sini=[7,8,9,10,11,12,13,14,15];  %Initial Stock Prices
Xini=2;         %Initial State
Ft=3;           %Time Horizon
%______________________________________________
dt=0.004;       %1 day in year unit, the time increment
sdt=sqrt(dt);
l=500000;       % No. (Even Number) of paths
n=50;           % No. of estimators required
m=ceil(Ft/dt);
nini=length(Sini);
F=zeros(n,nini); Estimate= zeros(1,nini); errorF=zeros(1,nini);
%Simulation starts here________________________
for k=1:n
 k
 Mka=zeros(l,nini);
 for i=1:l/2
     Ca=zeros(1,m); Da=zeros(1,m);
     Ciant=0;       Diant=0;
     T=zeros(1,m);  X=zeros(1,m);
     X(1)=Xini;
     a=0;
     j=0;
     while a<Ft
         j=j+1;
         T(j)=exprnd(1/L(X(j)));
         b=min(a+T(j),Ft);
            Ca(j)=(r-0.5*(sigma(X(j))^2))*(b-a);
            Da(j)=sigma(X(j))*normrnd(0,1)*sqrt(b-a);
         a=b;
         cumdist=CPM(X(j),:);
         X(j+1)=find(cumdist>rand(),1);
     end
     Caint=sum(Ca(:));
     Daint=sum(Da(:));
     % Application of Antithetic Simulation Method
     Mka(i,:)     =max(0,exp(Caint+Daint)*Sini-K); 
     Mka(i+l/2,:) =max(0,exp(Caint-Daint)*Sini-K);   
 end
 for ini=1:nini
     F(k,ini)=(exp(-r*Ft))*mean(Mka(:,ini));
 end 
end
for ini=1:nini
       Estimate(ini)=mean(F(:,ini));
       errorF(ini)=std(F(:,ini))/sqrt(n);
end
sprintf('Estimated option prices at regime %d for initial stocks %g, %g etc are', Xini, Sini(1),Sini(2))
Estimate
errorF