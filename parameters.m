%defines all constants as a script called in main.m

%nondimensionalized species transport coefficients
D_0=1; %Nondimensionalized base diffusivity. This should be left at 1.
D_p=1; %Relative diffusivity of positive species. D_pos = D_p/D_0
D_m=1; %Relative diffusivity of negative species. D_neg = D_p/D_0

kc_left=1; %butler-volmer charge transfer constants
jr_left=1;
kc_right=1;
jr_right=1;

%nondimensionalized electrical and polarization coefficients
delta=1; % ratio between the stern layer width and the debye length
%epsilon=0.01 % commented out here because it is being set in the run scripts. Epsilon needs to be set before this script is run.
lambda_s=delta*epsilon; % nondimensionalized stern layer width
epsilon_1=epsilon; % epsilon_{1, 2, 3, s} were intended to be used to have different values of epsilon in different parts of the domain, but for this work, they are all equal.
epsilon_2=epsilon;
epsilon_3=epsilon;
epsilon_s_left=epsilon;
epsilon_s_right=epsilon;

z_cp = 1; %charge number for cp
z_cm = -1; %charge number for cm