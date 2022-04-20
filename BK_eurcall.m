function [Val, V, S, tau]=BK_eurcall(r,sigma,S0,K,T,Nt)
    Smax = 2*K;
    Smin = 0;
    Ns = Nt/10;

    dt = T/Nt;
    ds = (Smax-Smin)/Ns;

    V(1:Ns+1, 1:Nt+1)=0.0;

    S=Smin+(0:Ns)*ds;
    tau=(0:Nt)*dt;
    V(1:Ns+1,1)=max(S-K,0);
    for i=1:Ns+1
        for j=2:Nt+1
            d1=(log(S(i)/K)+(r+sigma^2/2)*tau(j))/(sigma*sqrt(tau(j)));
            d2=d1-sigma*sqrt(tau(j));
            V(i,j)=S(i)*normcdf(d1)-K*exp(-r*tau(j))*normcdf(d2);
        end
    end
    ind=find(S==S0);
    Val=V(ind,Nt+1);
end


