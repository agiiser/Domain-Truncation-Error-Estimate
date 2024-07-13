clc;clear;%format long
%%%%%%%%%%%%%%%%%%%   START for Table 1 of our paper
K=10;smax=2*K;r=0.12;sig=0.1;t=0;T=2;B=2*K;
s=[7,8,9,10,11,12,13,14,15];w=zeros(length(s),1);v=zeros(length(s),1);

%%%%%%%%%%  Columns 2 and 3 of Table 1
for i=1:length(s)
[call,put]=blsprice(s(i),K,r,T,sig);
v(i)=call;
w(i)=barrier_exact_solution(s(i),T,K,B,r,sig);
end
disp('s, v and w are');
[s',v,w]

%%%%%%%%%% The following function is used above

function[eta]=barrier_exact_solution(s,TT,st,B,r,sig)
       bk=2*log(st/B);
       rp= r+ 0.5 *sig^2;
       rm= r- 0.5 *sig^2;
       ep= 1+2*r/sig^2;
       em=-1+2*r/sig^2;
       tm= TT;%     !   tm is time to expiry:=T-t
       dn= sig*sqrt(tm);%        ! dn is denominator
 
               sk = log(s/st);
               sb = log(s/B);
               term1= normcdf((sk+rp*tm)/dn)- normcdf((sb+rp*tm)/dn);
               term2= normcdf((sb-rp*tm)/dn)- normcdf((sk+bk-rp*tm)/dn);
               term3= normcdf((sk+rm*tm)/dn)- normcdf((sb+rm*tm)/dn);
               term4= normcdf((sb-rm*tm)/dn)- normcdf((sk+bk-rm*tm)/dn);
               term= st*exp(-r*tm)*(term3-((B/s)^em)*term4);
               eta= s*(term1-((B/s)^ep)*term2)- term;
end