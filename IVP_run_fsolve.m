% t_end = 20; using fast version of concentration_implicit_RHS.m
% Newton
% Runtime: 3 minutes and 5.686863 seconds
% Elapsed time is 185.689118 seconds.
% fsolve
% Runtime: 3 minutes and 25.714496 seconds
% Elapsed time is 205.716222 seconds.
% semi-implicit
% Runtime: 0 minutes and 1.430940 seconds
% Elapsed time is 1.432333 seconds.

dt_max = 5
%
% Newton:
% Runtime: 85 minutes and 34.821204 seconds
% Elapsed time is 5134.821394 seconds.
% fsolve:
% Runtime: 90 minutes and 27.056331 seconds
% Elapsed time is 5427.056477 seconds.
% semi-implicit:
% Runtime: 0 minutes and 11.737650 seconds
% Elapsed time is 11.737809 seconds.
% 
% max(NF_implicit(2:end))
% ans =
%    2.3384e-10
% min(NF_implicit(2:end))
% ans =
%    3.7880e-12
% mean(NF_implicit(2:end))
% ans =
%    1.1899e-11
% max(NF_implicit_fsolve(2:end))
% ans =
%    2.2844e-10
% min(NF_implicit_fsolve(2:end))
% ans =
%    4.7343e-12
% mean(NF_implicit_fsolve(2:end))
% ans =
%    1.0928e-11
%
% When running to final time 100, N = 91, tol = 1.e-6
% Newton takes 114 minutes and 22.942542 seconds
% fsolve takes 126 minutes and 56.126083 seconds
% vssbdf2 takes 17.367670 seconds

% dt_max = 1
%
% Runtime: 56 minutes and 51.379575 seconds
% Elapsed time is 3411.379744 seconds.
% Runtime: 53 minutes and 38.222709 seconds
% Elapsed time is 3218.222877 seconds.
% Runtime: 0 minutes and 12.659449 seconds
% Elapsed time is 12.659593 seconds.

clear;



epsilon = .5;
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
t_end=100;
t_end = 20;
dt=1e-4; %initial dt
% dt = .0001;
t_save=0.1; %save every multiple of this parameter
%set up error tolerances for adaptive stepping
tol=1e-6;
% tol = 1.e-4; % make it smaller for the Newton method...
range=tol/3;
% range=Inf;
dt_max=1; %biggest i'm willing to let dt get
% dt_max = .04;
% dt_max = 5;
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

clear t err_save cp cm phi phi_x_right
%set up initial conditions
n=1;
t(n)=0;
for i=1:N
    cp(i,n)=1+0.1*sin(2*pi*i/N);
    %     cp(i,n)=1;
    cm(i,n)=1+0.1*sin(2*pi*i/N);
    debug1(i,n)=1;
end
phi_x_right(n)=0;
phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
err_save(n)=0;
NF_run(n)=0;
start=tic;
run_fully_implicit
toc(start)
t_implicit = t;
err_save_implicit = err_save;
cp_implicit = cp;
cm_implicit = cm;
phi_implicit = phi;
phi_x_right_implicit = phi_x_right;
NF_implicit = NF_run;

clear t err_save cp cm phi phi_x_right NF_run
%set up initial conditions
n=1;
t(n)=0;
for i=1:N
    cp(i,n)=1+0.1*sin(2*pi*i/N);
    %     cp(i,n)=1;
    cm(i,n)=1+0.1*sin(2*pi*i/N);
    debug1(i,n)=1;
end
phi_x_right(n)=0;
phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
err_save(n)=0;
start=tic;
run_fully_implicit_fsolve
toc(start)
t_implicit_fsolve = t;
err_save_implicit_fsolve = err_save;
cp_implicit_fsolve = cp;
cm_implicit_fsolve = cm;
phi_implicit_fsolve = phi;
phi_x_right_implicit_fsolve = phi_x_right;
NF_implicit_fsolve = NF_run;



clear t err_save cp cm phi phi_x_right
%set up initial conditions
n=1;
t(n)=0;
for i=1:N
    cp(i,n)=1+0.1*sin(2*pi*i/N);
    %     cp(i,n)=1;
    cm(i,n)=1+0.1*sin(2*pi*i/N);
    debug1(i,n)=1;
end
phi_x_right(n)=0;
phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
err_save(n)=0;
start=tic;
run
% toc(start)

% save compare_three_data_dtmax_1.mat


% figure(1)
% clf
% hold on
% plot(t(1:end-1),log10(diff(t)),'k')
% plot(t_implicit(1:end-1),log10(diff(t_implicit)),'r')
% plot(t_implicit_fsolve(1:end-1),log10(diff(t_implicit_fsolve)),'g')
% plot(t,log10(dt_max)*ones(size(t)),'--')
% dt_inf_0 = 0.024991129014342;
% dt_inf_1 = 0.022323241029796;
% dt_inf_2 = 0.019050813425106;
% dt_inf_3 = 0.015845235369629
% plot(t,log10(dt_inf_0)*ones(size(t)),':')
% plot(t,log10(dt_inf_1)*ones(size(t)),':')
% plot(t,log10(dt_inf_2)*ones(size(t)),':')
% plot(t,log10(dt_inf_3)*ones(size(t)),':')
% axis([0 60 -7 1])
% print -dps dt_time_varying_voltage.ps
% 
% figure(2)
% clf
% hold on
% plot(x,cp(:,end),'b')
% plot(x,cm(:,end),'b')
% plot(x,phi(:,end),'b')
% plot(x,cp_implicit(:,end),'r')
% plot(x,cm_implicit(:,end),'r')
% plot(x,phi_implicit(:,end),'r')
% plot(x,cp_implicit_fsolve(:,end),'g')
% plot(x,cm_implicit_fsolve(:,end),'g')
% plot(x,phi_implicit_fsolve(:,end),'g')
% 
% figure(3)
% clf
% hold on
% plot(t,log10(err_save),'b')
% plot(t_implicit,log10(err_save_implicit),'r')
% plot(t_implicit_fsolve,log10(err_save_implicit_fsolve),'g')
% plot(t,log10((tol+range)*ones(size(t))),'--')
% plot(t,log10((tol-range)*ones(size(t))),'--')
% 
% figure(4)
% clf
% hold on
% plot(t_implicit,log10(NF_implicit),'r')
% plot(t_implicit_fsolve,log10(NF_implicit_fsolve),'g')



% figure(1)
% clf
% hold on
% plot(t_implicit(1:end-1),diff(t_implicit));
% plot(t(1:end-1),diff(t),'r')
% % I see that the original scheme and the scheme where cp and cm are
% % computed implicitly (but phi is extrapolated forward and depends on cp
% % and cm) have the same thresholding in dt.
% %
% % plot(t,dt_max*ones(size(t)),'--');
% 
% figure(2)
% clf
% hold on
% plot(t_implicit,mean(cp_implicit))
% plot(t,mean(cp),'r')
% % the schemes seem to agree as far as the means of cp are concerned
% 
% figure(3)
% clf
% hold on
% plot(t_implicit,mean(cm_implicit))
% plot(t,mean(cm),'r')
% % the schemes seem to agree as far as the means of cp are concerned
% 
% cp_implicit_ss = cp_implicit(:,end);
% cm_implicit_ss = cm_implicit(:,end);
% for ii=1:length(t_implicit)-1
%     dp_implicit(:,ii) = cp_implicit(:,ii)-cp_implicit_ss;
%     dm_implicit(:,ii) = cm_implicit(:,ii)-cm_implicit_ss;
% end
% figure(4)
% clf
% hold on
% plot(t_implicit(1:end-1),log10(max(abs(dp_implicit))));
% plot(t_implicit(1:end-1),log10(max(abs(dm_implicit))),'r');
% 
% figure(5)
% clf
% hold on
% plot(x,cp_implicit(:,end))
% plot(x,cm_implicit(:,end),'r')
% plot(x,phi_implicit(:,end),'g')
% plot(x,cp(:,end),'o')
% plot(x,cm(:,end),'ro')
% plot(x,phi(:,end),'go')
% 
% % let's compute to find the steady-state profile for both schemes, down to
% % round-off.
% t_end = 10;
% dt_end = mean(diff(t(end-100:end)));
% dt = .8*dt_end;
% range = inf;
% cp_ss_approx = cp(:,end);
% cm_ss_approx = cm(:,end);
% 
% clear t err_save cp cm phi dp dm dp_implicit dm_implicit
% %set up initial conditions
% n=1;
% t(n)=0;
% cp(:,n) = cp_ss_approx + .001*sin(2*pi*x);
% cm(:,n) = cm_ss_approx + .001*sin(2*pi*x);
% % for i=1:N
% %     cp(i,n)=1+0.1*sin(2*pi*i/N);
% %     %     cp(i,n)=1;
% %     cm(i,n)=1+0.1*sin(2*pi*i/N);
% %     debug1(i,n)=1;
% % end
% phi_x_right(n)=0;
% phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
% err_save(n)=0;
% start=tic;
% run
% toc(start)
% 
% cp_ss = cp(:,end);
% cm_ss = cm(:,end);
% phi_ss = phi(:,end);
% phi_x_right_ss = phi_x_right(end);
% for ii=1:length(t)-1
%     dp(:,ii) = cp(:,ii)-cp_ss;
%     dm(:,ii) = cm(:,ii)-cm_ss;
% end
% figure(6)
% clf
% hold on
% plot(t(1:end-1),log10(max(abs(dp))));
% plot(t(1:end-1),log10(max(abs(dm))),'r');
% 
% clear t err_save cp cm phi dp dm dp_implicit dm_implicit
% %set up initial conditions
% n=1;
% t(n)=0;
% cp(:,n) = cp_ss_approx + .001*sin(2*pi*x);
% cm(:,n) = cm_ss_approx + .001*sin(2*pi*x);
% % for i=1:N
% %     cp(i,n)=1+0.1*sin(2*pi*i/N);
% %     %     cp(i,n)=1;
% %     cm(i,n)=1+0.1*sin(2*pi*i/N);
% %     debug1(i,n)=1;
% % end
% phi_x_right(n)=0;
% phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
% err_save(n)=0;
% start=tic;
% run_implicit
% toc(start)
% t_implicit = t;
% err_save_implicit = err_save;
% cp_implicit = cp;
% cm_implicit = cm;
% phi_implicit = phi;
% phi_x_right_implicit = phi_x_right;
% 
% 


% cp_implicit_ss = cp_implicit(:,end);
% cm_implicit_ss = cm_implicit(:,end);
% phi_implicit_ss = phi_implicit(:,end);
% phi_x_right_implicit_ss = phi_x_right(end);
% for ii=1:length(t)-1
%     dp_implicit(:,ii) = cp_implicit(:,ii)-cp_implicit_ss;
%     dm_implicit(:,ii) = cm_implicit(:,ii)-cm_implicit_ss;
% end
% figure(6)
% hold on
% plot(t_implicit(1:end-1),log10(max(abs(dp_implicit))),'--');
% plot(t_implicit(1:end-1),log10(max(abs(dm_implicit))),'r--');
% 
% % now that I have the steady-state profiles, compute with the adaptive mesh
% % one last time so that I can save all the data...
% t_end=6;
% range=tol/3;
% 
% clear t err_save cp cm phi phi_x_right
% %set up initial conditions
% n=1;
% t(n)=0;
% for i=1:N
%     cp(i,n)=1+0.1*sin(2*pi*i/N);
%     %     cp(i,n)=1;
%     cm(i,n)=1+0.1*sin(2*pi*i/N);
%     debug1(i,n)=1;
% end
% phi_x_right(n)=0;
% phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
% err_save(n)=0;
% start=tic;
% run_implicit
% toc(start)
% t_implicit = t;
% err_save_implicit = err_save;
% cp_implicit = cp;
% cm_implicit = cm;
% phi_implicit = phi;
% phi_x_right_implicit = phi_x_right;
% clear t err_save cp cm phi phi_x_right
% %set up initial conditions
% n=1;
% t(n)=0;
% for i=1:N
%     cp(i,n)=1+0.1*sin(2*pi*i/N);
%     %     cp(i,n)=1;
%     cm(i,n)=1+0.1*sin(2*pi*i/N);
%     debug1(i,n)=1;
% end
% phi_x_right(n)=0;
% phi(:,n)=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp(:,n),cm(:,n),z_cp,z_cm, phi_x_right(n),voltage(t(n)) );
% err_save(n)=0;
% start=tic;
% run
% toc(start)
% 
% save data.mat
% 
% % compare the steady states; they agree up to 10^(-11)
% figure(7)
% clf 
% hold on
% plot(x,(abs(cp_ss-cp_implicit_ss)));
% plot(x,(abs(cm_ss-cm_implicit_ss)),'r')
% plot(x,(abs(phi_ss-phi_implicit_ss)),'g')
% % do we have the same eigenmodes?  They have different eigenmodes.  I
% % wonder how the stability domain changes...
% figure(8)
% clf
% hold on
% y = cp(:,end)-cp_ss;
% y = y/y(end);
% plot(x,y)
% y = cp_implicit(:,end)-cp_implicit_ss;
% y = y/y(end);
% plot(x,y,'--')
% y = cm(:,end)-cm_ss;
% y = y/y(end);
% plot(x,y,'r')
% y = cm_implicit(:,end)-cm_implicit_ss;
% y = y/y(end);
% plot(x,y,'r--')
% 
% 
% 







