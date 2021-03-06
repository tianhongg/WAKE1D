function [Chi,Psi,Ez,ne,gama]=getChi(A,nxi,dxi)

gama=ones(nxi,1);
ne  =ones(nxi,1);
Psi=zeros(nxi,1);
Ez=zeros(nxi,1);


        
    %  xi loop
    for k=1:nxi-1
    %RK4----1
    gama(k)=0.5*(1+(1+Psi(k))^2+2*abs(A(k))^2)/(1+Psi(k));
    ne(k)=gama(k)/(1+Psi(k));
    kn1=(ne(k)-1)*dxi;
    
    if(k==1); ke1=0; else; ke1=Ez(k)*dxi; end
    Ez1 =Ez(k)+kn1*0.5;
    Ps1=Psi(k)+ke1*0.5;
    
    %RK4----2
    gama(k)=0.5*(1+(1+Ps1)^2+2*abs(A(k)*0.5+A(k+1)*0.5)^2)/(1+Ps1);
    ne(k)=gama(k)/(1+Ps1);
    kn2=(ne(k)-1)*dxi;
    ke2=Ez1*dxi;
    
    
    Ez2 =Ez(k)+kn2*0.5;
    Ps2=Psi(k)+ke2*0.5;
    
    %RK4----3
    gama(k)=0.5*(1+(1+Ps2)^2+2*abs(A(k)*0.5+A(k+1)*0.5)^2)/(1+Ps2);
    ne(k)=gama(k)/(1+Ps2);
    kn3=(ne(k)-1)*dxi;
    ke3=Ez2*dxi;
    
    
    Ez3 =Ez(k)+kn3;
    Ps3=Psi(k)+ke3;
    
    %RK4----4
    
    
    gama(k)=0.5*(1+(1+Ps3)^2+2*abs(A(k+1))^2)/(1+Ps3);
    ne(k)=gama(k)/(1+Ps3);
    kn4=(ne(k)-1)*dxi;
    ke4=Ez3*dxi;
    
    
    Psi(k+1)=Psi(k)+(ke1+2*ke2+2*ke3+ke4)/6;
    Ez(k+1) = Ez(k)+(kn1+2*kn2+2*kn3+kn4)/6;

    
    end
    
    Chi=1./(1+Psi(:));


end