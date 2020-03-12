clc
clear all

a=50;
b=10;
L=a+b;
h=input('Enter discretization size in cm: ');
N=L/h;

%%Constants for core region%%
sigmaC1=0.02935;
sigmaC2=0.10490;
nusigmafC1=0.000242;
nusigmafC2=0.155618;
DC1=1.4380;
DC2=0.3976;
sigmasC1=0;
sigmasC2=0.01563;

%%Constant for reflector region%%
sigmaR1=0.035411;
sigmaR2=0.031579;
nusigmafR1=0;
nusigmafR2=0;
DR1=1.871420;
DR2=0.283409;
sigmasR1=0;
sigmasR2=0.034340;



%%%%% THERMAL GROUP %%%%%

M2=zeros(N+1);
M2(1,1)=(DC2/h)+(sigmaC2*h/2);
M2(1,2)=-DC2/h;

for i=2:N
   if i*h<50+h
       M2(i,i-1)=-DC2/h;
       M2(i,i)=(2*DC2/h)+(sigmaC2*h);
       M2(i,i+1)=-DC2/h;
   end
   if i*h==50+h
       M2(i,i-1)=-DC2/h;
       M2(i,i)=((sigmaR2+sigmaC2)*h/2)+(DC2/h)+(DR2/h);
       M2(i,i+1)=-DR2/h;
   end
   if i*h>50+h
       M2(i,i-1)=-DR2/h;
       M2(i,i)=(2*DR2/h)+(sigmaR2*h);
       M2(i,i+1)=-DR2/h;
   end
end

M2(N+1,N)=-DR2/h;
M2(N+1,N+1)=((DR2/h)+(1/2)+(sigmaR2*h/2));



F2=zeros(N+1);
F2(1,1)=nusigmafC2*h/2;

for i=2:N+1
   if i*h<50+h
       F2(i,i)=nusigmafC2*h;
   end
   if i*h==50+h
       F2(i,i)=nusigmafC2*h/2;
   end
   if i*h>50+h
       F2(i,i)=0;
   end
end

S=zeros(N+1);
S(1,1)=sigmasC2*h/2;
for i=2:N
    if i*h<50+h
       S(i,i)=sigmasC2*h; 
    end
    if i*h==50+h
       S(i,i)=(sigmasC2+sigmasR2)*h/2; 
    end
    if i*h>50+h
       S(i,i)=sigmasR2*h;
    end
end
S(N+1,N+1)=sigmasR2*h/2;


%%%%% FAST GROUP %%%%%

M1=zeros(N+1);
M1(1,1)=(DC1/h)+(sigmaC1*h/2);
M1(1,2)=-DC1/h;

for i=2:N
   if i*h<50+h
       M1(i,i-1)=-DC1/h;
       M1(i,i)=(2*DC1/h)+(sigmaC1*h);
       M1(i,i+1)=-DC1/h;
   end
   if i*h==50+h
       M1(i,i-1)=-DC1/h;
       M1(i,i)=((sigmaR1+sigmaC1)*h/2)+(DC1/h)+(DR1/h);
       M1(i,i+1)=-DR1/h;
   end
   if i*h>50+h
       M1(i,i-1)=-DR1/h;
       M1(i,i)=(2*DR1/h)+(sigmaR1*h);
       M1(i,i+1)=-DR1/h;
   end
end

M1(N+1,N)=-DR1/h;
M1(N+1,N+1)=((DR1/h)+(1/2)+(sigmaR1*h/2));



F1=zeros(N+1);
F1(1,1)=nusigmafC1*h/2;

for i=2:N+1
   if i*h<50+h
       F1(i,i)=nusigmafC1*h;
   end
   if i*h==50+h
       F1(i,i)=nusigmafC1*h/2;
   end
   if i*h>50+h
       F1(i,i)=0;
   end
end


phiold1=ones(N+1,1);
phiold2=ones(N+1,1);
k_old=1;
for i=1:inf
    phinew1=inv(M1)*(1/k_old)*((F1*phiold1)+(F2*phiold2));
    phinew2=inv(M2)*(S*phinew1);
    phinew1core=phinew1(1:a/h);
    phinew2core=phinew2(1:a/h);
    phiold1core=phiold1(1:a/h);
    phiold2core=phiold2(1:a/h);
    k_new=(sum(nusigmafC1*phinew1core)+sum(nusigmafC2*phinew2core))/((1/k_old)*(sum(nusigmafC1*phiold1core)+sum(nusigmafC2*phiold2core)));
    phinew1=phinew1./max(phinew1);
    phinew2=phinew2./max(phinew2);
    if (abs(sum(phiold1)-sum(phinew1))> 10^-5) && (abs(sum(phiold2)-sum(phinew2))> 10^-5) && (abs((k_old-k_new))> 10^-5)
        phiold1=phinew1;
        phiold2=phinew2;
        k_old=k_new;
    else
        break;
    end
end

x=-60:h:60;
phinew1_plot=[flipud(phinew1(2:N+1,1));phinew1];
phinew2_plot=[flipud(phinew2(2:N+1,1));phinew2];
plot(x,phinew1_plot)
hold on
plot(x,phinew2_plot, '--')
title(sprintf('Flux vs x for %d cm mesh size',h))
ylabel('Flux (n/cm^2.s)')
xlabel('x (cm)')
legend('Fast Flux', 'Thermal Flux', 'Location', 'northwest')