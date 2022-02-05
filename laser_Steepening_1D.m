clear
%close all
clc


k0=10;
conver=2*pi/k0;

sigma=6*conver;

xi_max=100*conver;
T_max=2000*conver;

dxi=0.01;
dt =0.2;

nxi=floor(xi_max/dxi);
nt =floor(T_max/dt);


Xi=[0:nxi-1]*dxi;
TT=[0:nt-1]*dt;


A   =zeros(nxi,nt);
Chi =ones(nxi,nt);

A(:,1)=6*exp(-(Xi-sigma*4).^2/sigma^2);

% time loop
for tt=1:nt-1
        
    %  xi loop
    for k=1:nxi-1
    %RK4----1
    gama(k,tt)=0.5*(1+(1+Psi(k,tt))^2+2*abs(A(k,tt))^2)/(1+Psi(k,tt));
    ne(k,tt)=gama(k,tt)/(1+Psi(k,tt));
    kn1=(ne(k,tt)-1)*dxi;
    
    if(k==1); ke1=0; else; ke1=Ez(k,tt)*dxi; end
    Ez1 =Ez(k,tt)+kn1*0.5;
    Ps1=Psi(k,tt)+ke1*0.5;
    
    %RK4----2
    gama(k,tt)=0.5*(1+(1+Ps1)^2+2*abs(A(k,tt)*0.5+A(k+1,tt)*0.5)^2)/(1+Ps1);
    ne(k,tt)=gama(k,tt)/(1+Ps1);
    kn2=(ne(k,tt)-1)*dxi;
    ke2=Ez1*dxi;
    
    
    Ez2 =Ez(k,tt)+kn2*0.5;
    Ps2=Psi(k,tt)+ke2*0.5;
    
    %RK4----3
    gama(k,tt)=0.5*(1+(1+Ps2)^2+2*abs(A(k,tt)*0.5+A(k+1,tt)*0.5)^2)/(1+Ps2);
    ne(k,tt)=gama(k,tt)/(1+Ps2);
    kn3=(ne(k,tt)-1)*dxi;
    ke3=Ez2*dxi;
    
    
    Ez3 =Ez(k,tt)+kn3;
    Ps3=Psi(k,tt)+ke3;
    
    %RK4----4
    
    
    gama(k,tt)=0.5*(1+(1+Ps3)^2+2*abs(A(k+1,tt))^2)/(1+Ps3);
    ne(k,tt)=gama(k,tt)/(1+Ps3);
    kn4=(ne(k,tt)-1)*dxi;
    ke4=Ez3*dxi;
    
    
    Psi(k+1,tt)=Psi(k,tt)+(ke1+2*ke2+2*ke3+ke4)/6;
    Ez(k+1,tt) = Ez(k,tt)+(kn1+2*kn2+2*kn3+kn4)/6;

    
    end
    
    Chi(:,tt)=1./(1+Psi(:,tt));
  
    dd0=1.5;dd1=0;dd2=0;
    
    for k=1:nxi-1
    Aold=A(k,tt);
    SS=(Chi(k,tt)*0.25+1i*k0/dt-dd0/dxi/dt)*Aold+(-2*dd1+0.5d0*dd2)/dxi/dt;
    Anew=SS/(1i*k0/dt-Chi(k,tt)*0.25-dd0/dxi/dt);
    dd2=dd1;
    dd1=Anew-Aold;
    A(k,tt+1)=Anew;
    
    end
    

end




figure
nn=0;
for nt=[1,500,1000,1500]
nn=nn+1;   
subplot(2,2,nn)
hold on
plot(Xi,abs(A(:,nt)),'r','linewidth',3)
%plot(Xi,Psi(:,1),'k','linewidth',3)
plot(Xi,ne(:,nt),'k','linewidth',3)
plot(Xi,Ez(:,nt),'b','linewidth',3)
set(gca,'linewidth',2);
set(gca,'fontsize',28);
box on;
set(gca,'TickDir','out');
set(gca,'color','w');
set(gcf,'color','w');
xlabel('k_p\xi');
xlim([0,xi_max])
ylim([min(Ez(:,nt))*1.1,max(abs(A(:,nt)))*1.1])

end



figure;plot(sum(abs(A(:,:)).^2));