%													
%	Leaky integrate-and-fire model without noise	
%		

% CLEANUP
clear

%constants
V_th = -0.02;                           % spike threshold -> -50 mv
tau = 0.03;                             % taum -> R*C -> 10 ms or 0.030;
R_m = 3e7;                              % membrane resistance -> 10 Mega ohms or 3e7
C = 1;                                  % membrane capacitance

%time
dt = 0.001;                             % time in msec or 0.001
dur = 0.16;                             % duration of simulation in msec or 0.5
n_iter = floor(dur ./ dt)+1;            % number of iterations

%potential 
V_0 = -0.07;                            % The initial membrane potential -> 0.05mV
V = zeros(1, n_iter + 1);               % array to store all potential during simulation
V (1) = V_0;                            % store V_0 in the first index of V 
EL = -0.07;                             % resting potencial
V_reset = -0.07;                        % reset value

%injection current
I_inj = 2e-9;                        	% electrode current input --> 2e-9

%spike 
time = linspace(0, dur, n_iter + 1);    % time vector of simulation
t_spike = 0;                            % time of a spike occurrence
n_spikes = 0;                           % number of spikes
r_period = 0.02;                        % refractory period

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
    
    dV = (dt./tau) .* (EL - V(i) + R_m .* I_inj);                       % LIF equation
	V (i + 1) = V (i) + dV;                                         % update
    if (V (i + 1) > V_th)						% check if potential is high enough for spike to occurr
        if (n_spikes > 0)						% check if spike has already occurred
            if (time (i) >= (t_spike + r_period))		        % check if in refractory period							
                V (i + 1) = V_reset;					% reset voltage
                t_spike = time (i);					% mark spike time
                V(i) = -0;
                n_spikes = n_spikes+1;				        % increment number of spikes 
            end          
        else 								% no spike has occured 
        	V (i + 1) = V_reset;					% reset potential
            V(i) = -0;
            t_spike = time (i);						% mark spike time
        	n_spikes = n_spikes + 1;				% increment number of spikes
        end 
    end
end

% plotting the graph %
plot(time, V, 'color', 'blue');
title('Leaky integrate-and-fire model');
xlabel('Time (ms)')
ylabel('Voltage (mV)')
grid on
