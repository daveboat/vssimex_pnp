function [ j ] = cur( t )
%cur.m
%called by step.m for current boundary conditions
%takes in current time, t, and returns the imposed current, j

%      j=0.1*(tanh(100*(t-1.5))/2 + 1/2);

% j=0.5;
j = 0;

% j=-(0.1*(tanh(500*(t-2))/2 + 1/2) + 0.1*(tanh(500*(t-3))/2 + 1/2)+0.1*(tanh(500*(t-4))/2 + 1/2));
end
