function run_charge(data_python)
data.As=
data.as1=
data.as3=
data.C_nom=
data.ce_avg=
data.Cp=
data.cs1_max=
data.cs3_max=
data.diam=
data.Ds1_ref=
data.Ds3_ref=
data.Ea_Ds1=
data.Ea_Ds3=
data.Ea_k1=
data.Ea_k3=
data.eps1s=
data.eps3s=
data.F=
data.h=
data.height=
data.k1_ref=
data.k3_ref=
data.L=
data.R=
data.Rc=
data.rho=
data.Rs1=
data.Rs3=
data.SA_V=
data.T_amb=
data.T_ref=
data.thick1=
data.thick2=
data.thick3=
data.Vc=
data.x1_soc0=
data.x1_soc1=
data.y3_soc0=
data.y3_soc1=
% Building the model
N = 6;      % N+1 is the number of Chebyshev nodes used to discretize 
            % the diffusion equation in each particle
nodes           = get_nodes(data,N);        % Create the Chebyshev nodes
matrices_spm    = get_model(data,nodes,N);  % Create the SPM model
% Initial conditions
x1_init = data.x1_soc1;
y3_init = data.y3_soc1;
T_init  = data.T_amb;
initSPM = get_init(x1_init,y3_init,T_init,data,nodes,matrices_spm);
C_rate = 0.5;     % C-rate
I = @(t) C_rate*data.C_nom*ones(size(t));   % Applied current [A]
tspan = 0:10:30000;  	% Simulation time span
V_limit = [3.0 4.2];  	% Minimum and maximum voltage
event = @(t,y) cutOffVoltage(t,y,I,data,matrices_spm,V_limit);
fun = @(t,y) derivs_spm(t,y,I,data,matrices_spm);
opt = odeset('Events',event);
[result.time,result.state] = ode45(fun,tspan,initSPM.y0,opt);
result = get_postproc( result,data,nodes,matrices_spm,I);
% Voltage vs time
figure;
plot(result.time,result.voltage,'.-');
xlabel('Time [s]');
ylabel('Voltage [V]');
grid on;
end