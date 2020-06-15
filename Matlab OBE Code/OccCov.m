function x=OccCov(i,j,V,D)
    n=size(D,1);
    Vin=inv(V);
    D0=zeros(n);
    D0(1,1)=1;
    Eta=V*D0*Vin;
    x=0;
    for k=2:n
       Dk=zeros(n); 
       Dk(k,k)=1;
       Uk=V*Dk*Vin;
       x=x-(Eta(i,i)*Uk(i,j)+Eta(j,j)*Uk(j,i))/D(k,k);
    end       
    x=simplify(x);
end