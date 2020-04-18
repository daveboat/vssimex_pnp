% minimal version of running semi implicit, for profiling purposes
% this is the code which runs for figure 4 of the pnp paper

clear;
epsilon = .01;
parameters;

bc=1;


N1=30+1; 
N2=30;
N3=30;
N=N1+N2+N3;
L1=1/3;
L2=1/3;
L3=1/3;

[x,dx]=makemesh(N1,N2,N3,L1,L2,L3);
t_end = 1;
dt=1e-4; 
t_save=0.1; 
tol=1e-6;
range=tol/3;
dt_max=1;
dt_min=1e-12;
c_max=25;
eta_max=1.1;
eta_min=0.9;
c=0;

n=1;
t(n)=0;
for i=1:N
    cp(i,n)=1+0.1*sin(2*pi*i/N);
    cm(i,n)=1+0.1*sin(2*pi*i/N);
end
phi_x_right(n)=0;
phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
err_save(n)=0;
run