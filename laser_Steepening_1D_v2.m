clear
close all
clc


k0=15;
conver=2*pi/k0;


sigma=6*conver;

xi_max=100*conver;
T_max=3000*conver;

dxi=0.04;
dt =0.02;

nt_out=20;

nxi=floor(xi_max/dxi);
nt =floor(T_max/dt/nt_out);


Xi=[0:nxi-1]*dxi;
TT=[0:nt-1]*dt*nt_out;


A   =zeros(nxi,nt);
Chi =ones(nxi,nt);
Psi=zeros(nxi,nt);
Ez=zeros(nxi,nt);
ne=ones(nxi,nt);
gama=ones(nxi,nt);

A(:,1)=4*exp(-(Xi-sigma*8).^2/sigma^2);



% time loop
try
    close(h)
catch
end

h=waitbar(0,'Running...');
tt=1;
for tt1=1:(nt*nt_out-1)
waitbar(tt/nt)



[Chi(:,tt),Psi(:,tt),Ez(:,tt),ne(:,tt),gama(:,tt)]=getChi(A(:,tt),nxi,dxi);

[Anew]=push_laser(Chi(:,tt),A(:,tt),nxi,dxi,dt*0.5,k0);
[Chi(:,tt),Psi(:,tt),Ez(:,tt),ne(:,tt),gama(:,tt)]=getChi(Anew,nxi,dxi);
[Anew]=push_laser(Chi(:,tt),A(:,tt),nxi,dxi,dt,k0);

if(mod(tt1,nt_out)==0)
tt=tt+1; 
end

A(:,tt)=Anew;



% dd0=1.5;dd1=0;dd2=0;
%     
%     for k=1:nxi-1
%     Aold=A(k,tt);
%     SS=(Chi(k,tt)*0.25+1i*k0/dt-dd0/dxi/dt)*Aold+(-2*dd1+0.5d0*dd2)/dxi/dt;
%     Anew=SS/(1i*k0/dt-Chi(k,tt)*0.25-dd0/dxi/dt);
%     dd2=dd1;
%     dd1=Anew-Aold;
%     A(k,tt+1)=Anew;
%     
%     end
    

end
close(h)


for i=1:nt
 energy(i) =  sum(abs((A(1:nxi-1,i)+A(2:nxi,i))*0.5*1i*k0-(A(2:nxi,i)-A(1:nxi-1,i))/dxi).^2)*dxi;
 Echange(i)= -sum((Chi(1:nxi-1,i)+Chi(2:nxi,i))*0.5.*(abs(A(2:nxi,i)).^2-abs(A(1:nxi-1,i)).^2)/dxi)*dxi/2;
end


for i=1:nt
    Bl(2:nxi,i)=(A(1:nxi-1,i)+A(2:nxi,i))*0.5*1i*k0-(A(2:nxi,i)-A(1:nxi-1,i))/dxi;
    Bl(1,i)=A(2,i)*1i*k0;
end

figure
nn=0;
for ntt=[1,150,300,450]*2
nn=nn+1;   
subplot(2,2,nn)
hold on
plot(Xi,abs(Bl(:,ntt))/k0,'r','linewidth',3)
%plot(Xi,Psi(:,1),'k','linewidth',3)
plot(Xi,Chi(:,ntt),'k','linewidth',3)
plot(Xi,Ez(:,ntt),'b','linewidth',3)
set(gca,'linewidth',2);
set(gca,'fontsize',28);
box on;
set(gca,'TickDir','out');
set(gca,'color','w');
set(gcf,'color','w');
xlabel('k_p\xi');
xlim([0,xi_max])
ylim([min(Ez(:,ntt))*1.1,max(abs(Bl(:,ntt)))*1.1/k0])
title(['\omega_pt=',sprintf('%05.1f',ntt*nt_out*dt)],'fontsize',18,'fontweight','normal');
end

figure
plot(TT,gradient(energy,dt*nt_out));
set(gca,'linewidth',2);
set(gca,'fontsize',28);
box on;
set(gca,'TickDir','out');
set(gca,'color','w');
set(gcf,'color','w');
xlabel('\omega_pt');

hold on
plot(TT, Echange)
