% this is to find steady states
clear;


start=tic;

% de = (.1360-.1336)/50;
% EPSILON = .1336:de:.1360;
% 
% de = (.14-.1)/50;
% EPSILON = .1:de:.14;
%
% de = (.35-.01)/50;
% EPSILON = .01:de:.35;

de = (.2-.05)/50;
EPSILON = .05:de:.2;

% first find steady state using uninformed initial data; do this for a
% large value of epsilon
jj = length(EPSILON);
epsilon = EPSILON(jj);
%physical (i.e. not numerical) parameters script
parameters;


%type of boundary condition. bc=0 for current (set current in cur.m)
%                            bc=1 for voltage (set voltage in voltage.m)
bc=1;

%set up geometry
N1=30; % Dave had N1=N2=N3 but this doesn't give a uniform mesh.
N1=30+1; % if I'm going to have a uniform mesh then I need L1/(N1-1)=L2/N2=L3/N3.
N2=30;
N3=30;
N=N1+N2+N3;
L1=1/3;
L2=1/3;
L3=1/3;
%call my meshing function
[x,dx]=makemesh(N1,N2,N3,L1,L2,L3);

%set up timescale
n=1;
t(n)=0;
t_end=10;
dt=1e-4; %initial dt
% dt = .0001;
t_save=0.1; %save every multiple of this parameter
%set up error tolerances for adaptive stepping
tol=1e-6;
range=tol/3;
% range=Inf;
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
run_no_save

% for ii=1:length(t)-1
%     dp(:,ii) = cp(:,ii)-cp(:,end);
%     dm(:,ii) = cm(:,ii)-cm(:,end);
%     dphi(:,ii) = phi(:,ii)-phi(:,end);
% end
% figure(1)
% clf
% hold on
% plot(t(1:end-1),log10(max(abs(dp))+eps));
% plot(t(1:end-1),log10(max(abs(dm))+eps),'r');
% plot(t(1:end-1),log10(max(abs(dphi))+eps),'g');
% figure(1)
%
% figure(2)
% plot(t(1:end-1),diff(t),'.-')

% turn off the adaptive time-stepper
range = inf;
% set the time-step to be half of the final time-step
% DT = diff(t);
% dt = mean(DT(end-50:end))/2; % use this if used run.m, which keeps the entire history of t...
dt = dt/2;
% set the initial data
cp_init = cp(:,end);
cm_init = cm(:,end);
phi_x_right_init = phi_x_right(end);
phi_init = phi(:,end);
err_save_init = err_save(end);
clear cp cm phi phi_x_right j j_d j_f t err_save
n=1;
cp(:,n) = cp_init;
cm(:,n) = cm_init;
phi_x_right(n) = phi_x_right_init;
phi(:,n) = phi_init;
err_save(n) = err_save_init;
t(n) = 0;

run_occasional_save

cp_ss(:,jj) = cp(:,end);
cm_ss(:,jj) = cm(:,end);
phi_ss(:,jj) = phi(:,end);
phi_x_right_ss(:,jj) = phi_x_right(end);
j_ss(jj) = j(end);
j_d_ss(jj) = j_d(end);
j_f_ss(jj) = j_f(end);
epsilon
jj


% clear dp dm dphi
% for ii=1:length(t)-1
%     dp(:,ii) = cp(:,ii)-cp(:,end);
%     dm(:,ii) = cm(:,ii)-cm(:,end);
%     dphi(:,ii) = phi(:,ii)-phi(:,end);
% end
% figure(1)
% clf
% hold on
% plot(t(1:end-1),log10(max(abs(dp))+eps));
% plot(t(1:end-1),log10(max(abs(dm))+eps),'r');
% plot(t(1:end-1),log10(max(abs(dphi))+eps),'g');
% figure(1)
% pause(1)

% finally, now that I have a good guess for the steady state solution, I
% perturb it slightly and use the adaptive time-stepper so that I can get a
% good estimate for dt_end(1).

range = tol/3; % turn on the adaptive time-stepper
t_end = 100; % run for a longer time

clear cp cm phi phi_x_right j j_d j_f debug err_save t
cp(:,1) = cp_ss(:,jj) + +.01*sin(2*pi*x).^3;
cm(:,1) = cm_ss(:,jj) - +.01*sin(2*pi*x).^3;
debug(:,1) = ones(size(x));
phi_x_right(1) = phi_x_right_ss(jj);
n=1;
t(n) = 0;
phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
err_save(n)=0;
run

DT = diff(t);
dt_end(jj) = mean(DT(end-100:end));
LTE_vss(jj) = mean(err_save(end-100:end));

% figure(2)
% clf
% plot(t(1:end-1),DT);
% hold on
% plot(t(end-100:end),DT(end-100:end),'o')
% figure(2)
% pause(1)



for jj=length(EPSILON)-1:-1:1
    clear cp cm phi phi_x_right j j_d j_f debug err_save t
    tic
    epsilon = EPSILON(jj);
    %physical (i.e. not numerical) parameters script
    parameters;
    
    
    %type of boundary condition. bc=0 for current (set current in cur.m)
    %                            bc=1 for voltage (set voltage in voltage.m)
    bc=1;
    
    %set up geometry
    N1=30; % Dave had N1=N2=N3 but this doesn't give a uniform mesh.
    N1=30+1; % if I'm going to have a uniform mesh then I need L1/(N1-1)=L2/N2=L3/N3.
    N2=30;
    N3=30;
    N=N1+N2+N3;
    L1=1/3;
    L2=1/3;
    L3=1/3;
    %call my meshing function
    [x,dx]=makemesh(N1,N2,N3,L1,L2,L3);
    
    %set up timescale
    n=1;
    t(n)=0;
    t_end=10;
    dt=1e-4; %initial dt
    % dt = .0001;
    t_save=0.1; %save every multiple of this parameter
    %set up error tolerances for adaptive stepping
    tol=1e-6;
    range=tol/3;
    % range=Inf;
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
    
    %set up initial conditions by taking previous steady state
    cp(:,1) = cp_ss(:,jj+1);
    cm(:,1) = cm_ss(:,jj+1);
    debug(:,1) = ones(size(x));
    phi_x_right(1) = phi_x_right_ss(jj+1);
    phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
    err_save(n)=0;
    run_no_save
    
    % for ii=1:length(t)-1
    %     dp(:,ii) = cp(:,ii)-cp(:,end);
    %     dm(:,ii) = cm(:,ii)-cm(:,end);
    %     dphi(:,ii) = phi(:,ii)-phi(:,end);
    % end
    % figure(1)
    % clf
    % hold on
    % plot(t(1:end-1),log10(max(abs(dp))+eps));
    % plot(t(1:end-1),log10(max(abs(dm))+eps),'r');
    % plot(t(1:end-1),log10(max(abs(dphi))+eps),'g');
    % figure(1)
    %
    % figure(2)
    % plot(t(1:end-1),diff(t),'.-')
    
    % turn off the adaptive time-stepper
    range = inf;
    % set the time-step to be half of the final time-step
    % DT = diff(t);
    % dt = mean(DT(end-50:end))/2; % use this if used run.m, which keeps the entire history of t...
    dt = dt/2;
    % set the initial data
    cp_init = cp(:,end);
    cm_init = cm(:,end);
    phi_x_right_init = phi_x_right(end);
    phi_init = phi(:,end);
    err_save_init = err_save(end);
    clear cp cm phi phi_x_right j j_d j_f debug err_save t

    n=1;
    cp(:,n) = cp_init;
    cm(:,n) = cm_init;
    phi_x_right(n) = phi_x_right_init;
    phi(:,n) = phi_init;
    err_save(n) = err_save_init;
    t(n) = 0;
    
    run_occasional_save
    
    cp_ss(:,jj) = cp(:,end);
    cm_ss(:,jj) = cm(:,end);
    phi_ss(:,jj) = phi(:,end);
    phi_x_right_ss(:,jj) = phi_x_right(end);
    j_ss(jj) = j(end);
    j_d_ss(jj) = j_d(end);
    j_f_ss(jj) = j_f(end);
    run_time(jj) = toc;
    DT(jj) = dt;
    
    epsilon
    jj
    dt
    
    
    %     clear dp dm dphi
    %     for ii=1:length(t)-1
    %         dp(:,ii) = cp(:,ii)-cp(:,end);
    %         dm(:,ii) = cm(:,ii)-cm(:,end);
    %         dphi(:,ii) = phi(:,ii)-phi(:,end);
    %     end
    %     figure(1)
    %     clf
    %     hold on
    %     plot(t(1:end-1),log10(max(abs(dp))+eps));
    %     plot(t(1:end-1),log10(max(abs(dm))+eps),'r');
    %     plot(t(1:end-1),log10(max(abs(dphi))+eps),'g');
    %     figure(1)
    %     pause(1)
    
    % Now recompute each problem using an adaptive time-stepper; this will give
    % us an estimate of dt_thresh (we hope)  For the initial data, I'm taking
    % the steady state with a perturbation.  The perturbation is chosen so that
    % it's an odd function (and so the integral of the initial data is the same
    % as the integral of the initial data) and so that the perturbation and the
    % derivative of the perturbation are zero at x=0 and x=1.  This means that
    % the initial data satisfies the boundary conditions and I don't need to
    % take sudden short time-steps to resolve a fast transient in time right
    % after t=0.
    range = tol/3; % turn on the adaptive time-stepper
    t_end = 100; % run for a longer time
    
    clear cp cm phi phi_x_right j j_d j_f debug err_save t
    cp(:,1) = cp_ss(:,jj) + +.01*sin(2*pi*x).^3;
    cm(:,1) = cm_ss(:,jj) - +.01*sin(2*pi*x).^3;
    debug(:,1) = ones(size(x));
    phi_x_right(1) = phi_x_right_ss(jj);
    n=1;
    t(n) = 0;
    phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
    err_save(n)=0;
    run
    
    DT = diff(t);
    dt_end(jj) = mean(DT(end-100:end));
    LTE_vss(jj) = mean(err_save(end-100:end));
    
    figure(2)
    clf
    plot(t(1:end-1),DT);
    hold on
    plot(t(end-100:end),DT(end-100:end),'o')
    figure(2)
    pause(1)    
end

save steady_state_data.mat
