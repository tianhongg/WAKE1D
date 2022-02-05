function [Anew]=push_laser(Chi,A,nxi,dxi,dt,k0)

bb=zeros(1,nxi);
dd=zeros(1,nxi);

aa=+0.5/dxi/dt*ones(1,nxi);
cc=-0.5/dxi/dt*ones(1,nxi);
aa(1)=0;cc(nxi)=0;

for k=1:nxi
    
bb(k)=1i*k0/dt-Chi(k)*0.25;

if(k<nxi&&k>1)
dd(k)=(Chi(k)*0.25+1i*k0/dt)*A(k)-(A(k+1)-A(k-1))/2/dxi/dt;
end

end

dd(1)=(Chi(1)*0.25+1i*k0/dt)*A(1)-(A(2))/2/dxi/dt;
dd(nxi)=(Chi(nxi)*0.25+1i*k0/dt)*A(nxi)-(-A(nxi-1))/2/dxi/dt;

[Anew(:)]=solve_tridiag(aa,bb,cc,dd,nxi);

end