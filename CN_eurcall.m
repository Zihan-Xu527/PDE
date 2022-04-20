function [Val, V, S, tau]=CN_eurcall(r,sigma,S0,K,T,Nt)

    Smax = 2*K;
    Smin = 0;
    Ns = Nt/10;
%     Ns = ceil(sqrt(0.9*Nt/(sigma*sigma*T)));

    dt = T/Nt;
    ds = (Smax-Smin)/Ns;

    V(1:Ns+1, 1:Nt+1)=0.0;

    S=Smin+(0:Ns)*ds;
    tau=(0:Nt)*dt;

    V(1:Ns+1,1)=max(S-K,0); % European Call contract function

    V(1,1:Nt+1)=0;
    V(Ns+1,1:Nt+1)=Smax-K*exp(-r*tau);


    A=zeros(Ns-1,1);
    B1=zeros(Ns-1,1);
    B2=zeros(Ns-1,1);
    C=zeros(Ns-1,1);
    
    for i=1:Ns-1
            A(i)=0.25*dt*(r*i-sigma*sigma*i*i);
            
            B1(i)=1+0.5*dt*(sigma*sigma*i*i+r);
            B2(i)=1-0.5*dt*(sigma*sigma*i*i+r);
            
            C(i)=-0.25*dt*(sigma*sigma*i*i+r*i);
    end
    M1=spdiags([C B1 A], [-1 0 1], Ns-1, Ns-1);
    M1=M1';
    M2=spdiags([-A B2 -C], [0 1 2], Ns-1, Ns+1);
    g=zeros(Ns-1,1);

    for j=1:Nt
        RHS=M2*V(:,j);
        g(1)=A(1)*V(1,j+1);
        g(Ns-1)=C(Ns-1)*V(Ns+1,j+1);
        RHS=RHS-g;
        V(2:Ns,j+1)=M1\RHS;
    end

    ind=find(S==S0);
    Val=V(ind,Nt+1);
end
    