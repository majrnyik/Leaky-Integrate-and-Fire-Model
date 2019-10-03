%													
%	Leaky integrate-and-fire model without noise	
%		

% resting potential -> -70 mv to +30 mv
% threshold potential/stimulus -> 50 mv less negative than resting pot
% membrane time constant = r_m * c_m -> around 10ms
% the neuron that just fired goes to a refractory state in which it doesn't respond
%refractory period -> limits the firing rate
% I (t) = Cm(dVm/dt)
% I (t) = - (Vm(t)/Rm) = Cm(dVm/dt)
% RI(t) = - Vm (t) = Omega or R*C (dVm/dt) where omega is R.Cm 

%% CLEANUP
clear
clc
clf

%constants
V_th = -0.02;                           % spike threshold -> -50 mv
tau = 0.02;                             % taum -> R*C -> 10 ms or 0.030;
R_m = 3e7;                              % membrane resistance -> 10 Mega ohms or 3e7
C = 1;                                  % membrane capacitance

%time
dt = 0.0001;                           % time in msec or 0.001
dur = 0.2;                             % duration of simulation in msec or 0.5
n_iter = floor(dur ./ dt)+1;           % number of iterations

%potential 
V = zeros(1, n_iter + 1);               % array to store all potential during simulation                          
V_0 = -0.07;                            % initial voltage
EL = V_0;                               % resting potencial
V_reset = V_0 ;                         % reset value
V (1) = V_0;                            % store V_0 in the first index of V

%injected current
I_inj = zeros(1, n_iter + 1);								 
I_inj(1) = 0;                        				% electrode current input --> 2e-9
I_noise = 2e-9 .* random('Normal', 1, dt+1, [1, n_iter]);	% initial current
%I_noise = I_inj(1) .* randn(1, n_iter + 1); 			% se n tem mean eh white noise (gaussian normal distr)

%spike 
time = linspace(0, dur, n_iter + 1);    % time vector of simulation
t_spike = 0;                            % time of a spike occurrence
n_spikes = 0;                           % number of spikes
r_period = 0.002;                       % refractory period

% using Euler's method for diferential equations we have the following
%1. define f(x,y)
%2. input t_0 and y_0
%3. input step size h and the number of steps n
%4. from j = 1 to n do
%	a) m = f(x_0, y_0)
%	b) y_1 = y_0 + h * m;
%	c) t_1 = t_0 + h
%	d) print t_1 and y_1
%	e) t_0 = t_1
%	f) y_0 = y_1
%4. end

for i = 1 : n_iter
    
    dV = ((EL - V(i) + I_inj(1).*R_m + I_noise(i).*R_m) .* dt) ./ tau;       	% LIF equation w/ noise
	V (i + 1) = V (i) + dV;                                             	% update
    if (V (i + 1) > V_th)							% check if potential is high enough for spike to occurr
        if (n_spikes > 0)							% check if spike has already occurred
            if (time (i) >= (t_spike + r_period))		                % check if in refractory period							
                V (i + 1) = V_reset;					        % reset voltage
                t_spike = time (i);						% mark spike time
                V(i) = 0;
                n_spikes = n_spikes+1;				                % increment number of spikes 
            end          
        else 									% no spike has occured 
        	V (i + 1) = V_reset;						% reset potential
            V(i) = 0;
            t_spike = time (i);							% mark spike time
        	n_spikes = n_spikes + 1;					% increment number of spikes
        end 
    end
end

% plotting the graph %
plot(time, V, 'color', 'blue');
title('Leaky integrate-and-fire model');
xlabel('Time (ms)')
ylabel('Voltage (mV)')
grid on
