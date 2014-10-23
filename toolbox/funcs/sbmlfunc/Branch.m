function xdot = Branch(x_values, time,SM)
% function Branch takes
%
% either	1) no arguments
%       	    and returns a vector of the initial values
%
% or    	2) time - the elapsed time since the beginning of the reactions
%       	   x_values    - vector of the current values of the variables
%       	    and returns a vector of the rate of change of value of each of the variables
%
% Branch can be used with MATLABs odeN functions as 
%
%	[t,x] = ode23(@Branch, [0, t_end], Branch)
%
%			where  t_end is the end time of the simulation
%
%--------------------------------------------------------
% output vector


%--------------------------------------------------------
% compartment values

compartmentOne = 1;

%--------------------------------------------------------
% parameter values

k1 = 0;
k2 = 0;
k3 = 0;

%--------------------------------------------------------
% initial values of variables - these may be overridden by assignment rules
% NOTE: any use of initialAssignments has been considered in calculating the initial values

	% floating variable values
	S1 = x_values(1);
	X0 = x_values(2);
	X1 = x_values(3);
	X2 = x_values(4);

%--------------------------------------------------------
% assignment rules

%--------------------------------------------------------
% algebraic rules

%--------------------------------------------------------
% calculate concentration values


	% rate equations
	V(1) = (k1*X0);
	V(2) = (k2*S1);
	V(3) = (k3*S1);
	xdot=SM*V;
end
