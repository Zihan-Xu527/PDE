% Evaluate an European Call option by using an explicit method
clear all; close all; clc;
%Parameters of the problem:
r=0.02;
sigma = 0.25;
S0 = 10;
K = 10;
T = 1.; 

Nt = [200, 400, 800,1600]; % number of time step
h = T./Nt;
n=length(Nt);
k_ex = zeros(size(Nt));
k_im = zeros(size(Nt));
k_CN = zeros(size(Nt));
Np = 200;

[Val_true, V_true, S_true, tau_true]=BK_eurcall(r, sigma,S0, K, T, Nt(n));
val_MC = MC_eurcall(r, sigma, S0, K, T, Nt(1), Np);


err_ex = zeros(size(Nt));
err_im = zeros(size(Nt));
err_CN = zeros(size(Nt));

for i=1:n
    
    [Val_ex, V_ex, S_ex, tau_ex]=explicit_eurcall(r, sigma,S0, K, T, Nt(i));
    k_ex(i) = S_ex(2);
    err_ex(i) = abs(Val_ex-Val_true);

    [Val_im, V_im, S_im, tau_im]=implicit_eurcall(r, sigma,S0, K, T, Nt(i));
    k_im(i) = S_im(2);
    err_im(i) = abs(Val_im-Val_true);

    [Val_CN, V_CN, S_CN, tau_CN]=CN_eurcall(r, sigma,S0, K, T, Nt(i));
    k_CN(i) = S_CN(2);
    err_CN(i) = abs(Val_CN-Val_true);

end

figure()
subplot(2,2,1)
plot(S_true,V_true(:,1),'r-',S_true,V_true(:,round(Nt(n)/2)),'g-',S_true,V_true(:,Nt(n)+1),'b-');
xlabel('S');
ylabel('V(S,tau)');
title('Exact Solution');
legend('V(T,S(T))', 'V(T/2,S(T/2))','V(0,S(0))', 'location', 'Northwest')
subplot(2,2,2)
plot(S_ex,V_ex(:,1),'r-',S_ex,V_ex(:,round(Nt(n)/2)),'g-',S_ex,V_ex(:,Nt(n)+1),'b-');
xlabel('S');
ylabel('V(S,tau)');
title('Explicit FD');
legend('V(T,S(T))', 'V(T/2,S(T/2))','V(0,S(0))', 'location', 'Northwest')
subplot(2,2,3)
plot(S_im,V_im(:,1),'r-',S_im,V_im(:,round(Nt(n)/2)),'g-',S_im,V_im(:,Nt(n)+1),'b-');
xlabel('S');
ylabel('V(S,tau)');
title('Fully implicit FD');
legend('V(T,S(T))', 'V(T/2,S(T/2))','V(0,S(0))', 'location', 'Northwest')
subplot(2,2,4)
plot(S_CN,V_CN(:,1),'r-',S_CN,V_CN(:,round(Nt(n)/2)),'g-',S_CN,V_CN(:,Nt(n)+1),'b-');
xlabel('S');
ylabel('V(S,tau)');
title('Crank Nicolson');
legend('V(T,S(T))', 'V(T/2,S(T/2))','V(0,S(0))', 'location', 'Northwest')
suptitle('European Call Option value, at t=0, T/2 and T')

figure()
subplot(2,2,1)
mesh(tau_true, S_true, V_true);
title('Exact solution')
xlabel('tau')
ylabel('S')
subplot(2,2,2)
mesh(tau_ex, S_ex, V_ex);
title('Explicit Method')
xlabel('tau')
ylabel('S')
subplot(2,2,3)
mesh(tau_im, S_im, V_im);
title('Implicit Method')
xlabel('tau')
ylabel('S')
subplot(2,2,4)
mesh(tau_CN, S_CN, V_CN);
title('Crank Nicolson Method')
xlabel('tau')
ylabel('S')
suptitle('European Call Option value, V(S, tau)')

figure()
loglog(h,err_ex,'o-', h,err_im,'x-', h,err_CN,'*-' )
% axis([.5*min(h) 1.5*max(h) ])
legend('Explicit Method', 'Implicit Method', 'Crank Nicolson Method')
title('log-log plot of errors vs. time step')

figure()
loglog(k_ex,err_ex,'o-', k_im,err_im,'x-', k_CN,err_CN,'*-' )
% axis([.5*min(h) 1.5*max(h) ])
legend('Explicit Method', 'Implicit Method', 'Crank Nicolson Method')
title('log-log plot of errors vs. spatial step')

error_table(h,err_ex)
error_table(k_ex,err_ex)
% error_loglog(h,err_ex)
% 
error_table(h,err_im)
error_table(k_im,err_im)
% error_loglog(h,err_im)
% 
error_table(h,err_CN)
error_table(k_CN,err_CN)
% error_loglog(h,err_CN)


