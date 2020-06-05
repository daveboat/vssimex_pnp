% if dt=0.00005;
%
% Using the ghost-point method in concentration, the ratios are
% Ratios: 4.1467, 4.0752, 4.0383, 4.0193, 4.0097
%
% Using the direct method with a two-point stencil in concentration.m, the
% ratios are
% Ratios: 6.1557, 0.7042, 1.5263, 1.7928, 1.9027
%
% Using the direct method with a three-ponit stencil in concentration.m,
% the ratios are
% Ratios: 18.7968, 0.6775, 1.2326, 1.6843, 1.8554

% if dt=0.000005;
%
% Using the ghost-point method in concentration, the ratios are
% Ratios: 3.9394, 3.9702, 3.9852, 3.9926, 3.9963
%
% Using the direct method with a two-point stencil in concentration.m, the
% ratios are
% Ratios: 1.5975, 1.8185, 1.9136, 1.9578, 1.9791
%
% Using the direct method with a three-ponit stencil in concentration.m,
% the ratios are
% Ratios: 1.3676, 1.7280, 1.8729, 1.9385, 1.9697

% if dt=0.0000005;
%
% Using the ghost-point method in concentration, the ratios are
% Ratios: 3.9866, 3.9933, 3.9967, 3.9984, 3.9995
%
% Using the direct method with a two-point stencil in concentration.m, the
% ratios are
% Ratios: 1.9560, 1.9783, 1.9892, 1.9946, 1.9973
%
% Using the direct method with a three-ponit stencil in concentration.m,
% the ratios are
% Ratios: 1.9509, 1.9758, 1.9880, 1.9940, 1.9970
%
% NOTE: to use the three methods you have to first modify concentration.m to have
% the correct boundary conditions and then run this script.

clear;

start = tic;
dt=0.00005;
dt=0.000005;
dt=0.0000005;
t_end = 20*dt;


%physical (i.e. not numerical) parameters script
epsilon=0.01;
parameters;

%type of boundary condition. bc=0 for current (set current in cur.m)
%                            bc=1 for voltage (set voltage in voltage.m)
bc=0;

%set up geometry
N1=10; % Dave had N1=N2=N3 but this doesn't give a uniform mesh.
N1=10+1; % if I'm going to have a uniform mesh then I need L1/(N1-1)=L2/N2=L3/N3.
N2=10;
N3=10;
N=N1+N2+N3;
L1=1/3;
L2=1/3;
L3=1/3;
%call my meshing function
[x,dx]=makemesh(N1,N2,N3,L1,L2,L3);

%set up timescale
n=1;
t(n)=0;
% t_end=.001;
% dt=1e-6; %initial dt
% dt = .0001;
t_save=0.01; %save every multiple of this parameter
%set up error tolerances for adaptive stepping
tol=1e-6;
range=tol/3;
range=Inf;
dt_max=1; %biggest i'm willing to let dt get
%dt_max=Inf;
dt_min=1e-8; %smallest i'm willing to let dt get
c_max=25; %maximum number of convergence steps before i take dt=dt_min
eta_max=1.1; %eta_min dt_old < dt_new < eta_max dt_old
eta_min=0.9;
c=0; %convergence counter

%set up file name for saving
formatOut = 'mm_dd_yy';
date=datestr(now,formatOut);
file_name=sprintf('%s_eps=%g_delta=%g_D=%g_k=%g_endtime=%g.mat',date,epsilon_1,delta,D_p,kc_right,t_end);

debug1=zeros(N,1);
%set up initial conditions
for i=1:N
     cp(i,n)=1+0.1*sin(2*pi*i/N);
%     cp(i,n)=1;
    cm(i,n)=1+0.1*sin(2*pi*i/N);
    debug1(i,n)=1;
end
phi_x_right(n)=0;
phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
err_save(n)=0;
run;
cp_1=cp;
phi_x_right_1=phi_x_right;
t_1=t;

dt=dt/2;
n=1;
run;
cp_2=cp;
phi_x_right_2=phi_x_right;
t_2=t;

dt=dt/2;
n=1;
run;
cp_3=cp;
phi_x_right_3=phi_x_right;
t_3=t;

dt=dt/2;
n=1;
run;
cp_4=cp;
phi_x_right_4=phi_x_right;
t_4=t;

dt=dt/2;
n=1;
run;
cp_5=cp;
phi_x_right_5=phi_x_right;
t_5=t;

dt=dt/2;
n=1;
run;
cp_6=cp;
phi_x_right_6=phi_x_right;
t_6=t;

dt=dt/2;
n=1;
run;
cp_7=cp;
phi_x_right_7=phi_x_right;
t_7=t;

r_1=cat(1,cp_1(:,21),phi_x_right_1(21));
r_2=cat(1,cp_2(:,41),phi_x_right_2(41));
r_3=cat(1,cp_3(:,81),phi_x_right_3(81));
r_4=cat(1,cp_4(:,161),phi_x_right_4(161));
r_5=cat(1,cp_5(:,321),phi_x_right_5(321));
r_6=cat(1,cp_6(:,641),phi_x_right_6(641));
r_7=cat(1,cp_7(:,1281),phi_x_right_7(1281));

rat1=norm(r_1-r_2)/norm(r_2-r_3);
rat2=norm(r_2-r_3)/norm(r_3-r_4);
rat3=norm(r_3-r_4)/norm(r_4-r_5);
rat4=norm(r_4-r_5)/norm(r_5-r_6);
rat5=norm(r_5-r_6)/norm(r_6-r_7);

display(sprintf('Ratios: %.4f, %.4f, %.4f, %.4f, %.4f', rat1, rat2, rat3,rat4,rat5));


