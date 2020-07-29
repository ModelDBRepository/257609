%This script can be used to test the parameters and dynamics for a
%recurrent network, there is no read-out

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = 0;
tEnd = 5000;                      %Simulation in milli-seconds
tStep = 0.1;                      %0.1 millisecond time step

time = [tStart:tStep:tEnd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hardcoded = false; 
if hardcoded
    synfire = true; %do not forget to reduce the external input to exc rec netw neurons
    if synfire
        mult = 0.15;
        numClusters = 120;
    	EneuronNum = 120;
        IneuronNum = 30;
        neuronNum = 150;   
        createWeightMatrixEIF
        weightsEE = diag(800*ones(EneuronNum-1,1),1);
        weightsEE(end,1) = 800;
    else
        mult = 3; 
        EneuronNum  = mult*800;                 %Number of excitatory neurons in the network
        numClusters = mult*10;                  %Number of clusters
        IneuronNum  = round(0.25*EneuronNum);   %Number of inhibitory neurons in the network
        neuronNum   = EneuronNum + IneuronNum;  %Total number of neurons
        createWeightMatrixEIF                   %manually set the connectivities
    end
else
    A = load('balanced_network.mat'); %import FF weight matrix
    %A = load('SSA20b_2400.mat'); %import SSA weight matrix
    A = A.A;
    
    neuronNum = size(A,1);
    EneuronNum = 0.80*size(A,1);
    IneuronNum = 0.20*size(A,1);
    numClusters = EneuronNum/80; %80 neurons/cluster
    
    weightsEE = A(1:EneuronNum,1:EneuronNum);
    weightsEI = A(1:EneuronNum,EneuronNum+1:neuronNum);
    weightsIE = A(EneuronNum+1:neuronNum,1:EneuronNum);
    weightsII = A(EneuronNum+1:neuronNum,EneuronNum+1:neuronNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for both E and I Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vreset = -60;       %Reset for both exc and inh neurons
C = 300;            %capacitance
tau_abs = 5;        %refractory period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the E-Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vthres = 20;        %Spiking threshold for exc neurons
tau_E = 20;         %Membrane Time Constant
V_E = -70;          %resting potential 
DET = 2;            %slope of exponential
E_E = 0;            %reversal potential
V_T = -52;          %threshold potential (the spiking threshold for inh neurons)
A_T = 10;           %post spike threshold potential increase
tau_T = 30;         %adaptive threshold time scale
EVthreshold = V_T*ones(1,EneuronNum); %neuronal threshold vector for all exc neurons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the I-Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_I = 20;         %Membrane Time Constant
V_I = -62;          %resting potential
E_I = -75;          %reversal potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the short term plasticity/adaptation (E only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_w = 100;                            %adaptation time constant
a = 0;                                  %adaptation slope
b = 1000;                               %adaptation amplitude
w = a*(Vreset-V_E)*ones(1,EneuronNum);  %adaptation vector for all excitatory neurons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the synapses and neuronal conductances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tauedecay = 6;      %decay time for e-synapses
tauerise = 1;       %rise time of e-synapses
tauidecay = 2;      %decay time for i-synapses
tauirise = 0.5;     %rise time of i-synapses

xedecay = zeros(1,neuronNum);
xerise = zeros(1,neuronNum);
xidecay = zeros(1,neuronNum);
xirise = zeros(1,neuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% External input to both E and I neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rex = 4.500;        %external rate to E-neurons        
rix = 2.250;        %external rate to I-neurons
Jeex = 1.6;         %weights for ee external input
Jiex = 1.52;        %weights for ie external input

nextx = zeros(1,neuronNum);      %vector containing the next external input spike times
nextx(1,1:EneuronNum) =  exprnd(1,1,EneuronNum)/rex;         
nextx(1,1+EneuronNum:end) = exprnd(1,1,IneuronNum)/rix;
rx = zeros(1,neuronNum);
rx(1,1:EneuronNum) = rex;
rx(1,EneuronNum+1:end) = rix;
forwardInputsEPrev = zeros(1,neuronNum);
forwardInputsIPrev = zeros(1,neuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulating Network with EIF neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rast = zeros(neuronNum,(tEnd - tStart)/tStep + 1);          %Matrix storing spike times for raster plots
rast_binary = zeros(neuronNum,(tEnd - tStart)/tStep + 1);   %same but with binary numbers
lastAP  = -50 * ones(1,neuronNum);                          %last action potential for refractor period calculation (just big number negative put)

memVol = Vreset+(V_T-Vreset)*rand(neuronNum,(tEnd - tStart)/tStep + 1);

for i =2:(tEnd - tStart)/tStep

    forwardInputsE = zeros(1,neuronNum);
    forwardInputsI = zeros(1,neuronNum);

    for j = 1:neuronNum
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %EXTERNAL INPUT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        while i*tStep > nextx(j)
            nextx(j) = nextx(j) + exprnd(1)/rx(j);
            if j <= EneuronNum
                forwardInputsEPrev(j) = forwardInputsEPrev(j) + Jeex;
            else
                forwardInputsEPrev(j) = forwardInputsEPrev(j) + Jiex;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CONNCECTIVITY CALCULATIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        xerise(j) = xerise(j) -tStep*xerise(j)/tauerise + forwardInputsEPrev(j);
        xedecay(j) = xedecay(j) -tStep*xedecay(j)/tauedecay + forwardInputsEPrev(j);
        xirise(j) = xirise(j) -tStep*xirise(j)/tauirise + forwardInputsIPrev(j);
        xidecay(j) = xidecay(j) -tStep*xidecay(j)/tauidecay + forwardInputsIPrev(j);
        
        gE = (xedecay(j) - xerise(j))/(tauedecay - tauerise);
        gI = (xidecay(j) - xirise(j))/(tauidecay - tauirise);
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %EXCITATORY NEURONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(j <= EneuronNum)
                        
            w(j) = w(j) + (tStep/tau_w)*(a*(memVol(j,i-1) - V_E) - w(j));           %adaptation current            
            EVthreshold(j) = EVthreshold(j) + (tStep/tau_T)*(V_T - EVthreshold(j));  %adapting threshold
            
            %cell dynamics
            v = memVol(j,i-1) + (tStep/tau_E)*(-memVol(j,i-1) + V_E + DET*exp((memVol(j,i-1)-EVthreshold(j))/DET)) ...
                + (tStep/C)*(gE*(E_E - memVol(j,i-1)) + gI*(E_I - memVol(j,i-1)) - w(j));
            
            
            if ((lastAP(j) + tau_abs/tStep)>=i)   %Refractory Period
                v = Vreset;
            end
            
            if (v > Vthres)           %Fire if exceed threshold
                v = Vreset;
                lastAP(j) = i;
                rast(j,i) = j;
                rast_binary(j,i) = 1;
                
                forwardInputsE = forwardInputsE + [weightsEE(:,j);weightsIE(:,j)]';
                
                EVthreshold(j) = EVthreshold(j) + A_T;
                w(j) = w(j) + b;
            end
            
            memVol(j,i) = v;
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %INHIBITORY NEURONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (j > EneuronNum)
            
            %cell dynamics
            v = memVol(j,i-1) + (tStep/tau_I)*(-memVol(j,i-1) + V_I) + ...
                (tStep/C)*(gE*(E_E - memVol(j,i-1)) + gI*(E_I - memVol(j,i-1)));
            
            if ((lastAP(j) + tau_abs/tStep)>=i)     %Refractory Period
                v = Vreset;
            end
            
            if (v > V_T)                  %Fire if exceed threshold
                v = Vreset;
                lastAP(j) = i;
                rast(j,i) = j;
                rast_binary(j,i) = 1;
                
                forwardInputsI = forwardInputsI + [weightsEI(:,j-EneuronNum);weightsII(:,j-EneuronNum)]';
            end
            
            memVol(j,i) = v;
        end             
    
    end
    
    forwardInputsEPrev = forwardInputsE;
    forwardInputsIPrev = forwardInputsI;
    
end

%these plotting scripts are taken from Schaub et al. Emergence of slow-switching assemblies in structured neuronal networks (2015)
[lambdas, W]= get_eigenvalues_LIF(weightsEE,weightsIE,weightsEI,weightsII);
figure; plot(real(lambdas),imag(lambdas),'.','color',[0 0.3 0], 'MarkerSize',24)
set(gca,'FontSize',25)
set(gcf,'Color',[1 1 1])
figure;
plotRASTER
