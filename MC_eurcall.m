function Val = MC_eurcall( r, sigma, S0, K, T, Nt, Np)
    dt=T/Nt;
    S=zeros(Nt+1,Np);
    S(1,:)=S0*ones(Np,1)';
    figure();
    for i=1:Np
        for j=1:Nt
            S(j+1,i)=S(j,i)*exp((r-0.5*sigma^2)*dt+sigma*normrnd(0,1)*sqrt(dt));
        end
        plot((0:Nt)*dt, S(:,i));
        hold on
    end
    ST=S(Nt+1,:)';
    Phi=max(ST-K,0);
    Val=exp(-r*T)*mean(Phi);
end