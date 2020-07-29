clearvars -global
clear all

%This script sets the parameters and launches a training script
%Choices: Hardcoded weight matrix (synfire chain) or learned weight matrix
%Run generateSequence script to simulate spontaneous dynamics
%ABCBA supervisor is hardcoded in this script, supervisor of Fig 6 is loaded from data folder
%different lengths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = 0;
tEnd = 12000;                       %Simulation in milli-seconds 
tStep = 0.1;                        %0.1 millisecond time step
    
time = [tStart:tStep:tEnd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT MATRIX OF TEMPORAL BACKBONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hardcoded = false; 
thirty_clusters = true;
if hardcoded
    synfire = true;
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
    if thirty_clusters
        %matrix60b has 30 exc clusters
        A = load('FF60b_2400_nc80.mat'); %import weight matrix
        A = A.A;
    else
        %This matrix has 80 exc clusters, very large simulation: may run
        %into memory errors on local pc. Cluster recommended.
        A = load('FF180_6400b.mat');
        A = A.A;
    end
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
% SUPERVISOR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if thirty_clusters
    actionNeuronNum = 3;    %read-out neurons
    superNeuronNum = 3;     %supervisor neurons
    interNeuronNum = 3;     %interneurons
    
    %75 ms stimulations
    superLength = 400; %[ms]
    supervisor = ones(actionNeuronNum,superLength/tStep+1); %baseline input
    supervisor(1,25/tStep:100/tStep) = 10.;
    supervisor(1,325/tStep+1:400/tStep) = 10.;
    supervisor(2,100/tStep+1:175/tStep) = 10.;
    supervisor(2,250/tStep+1:325/tStep) = 10.;
    supervisor(3,175/tStep+1:250/tStep) = 10.;
else
    supervisor = load('supervisor_600.mat'); %600 ms part of a bird song
    supervisor = supervisor.supervisor;
    actionNeuronNum = size(supervisor,1);
    superNeuronNum = actionNeuronNum;
    interNeuronNum = actionNeuronNum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT MATRIX TO ACTION READ OUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initially all zero read-out weights, all-to-all connected
W_AE = zeros(actionNeuronNum,EneuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for both E and I Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vreset = -60;       %Reset for both exc and inh neurons
C = 300;            %capacitance
tau_abs = 5;        %refractory period
tau_absA = 1;       %refractory period for read-out neurons
tau_absS = 1;       %refractory period for supervisor neurons
tau_absH = 1;       %refractory period for interneurons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the E-Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vthres = 20;        %Spiking threshold for exc neurons
tau_E = 20;         %Membrane time constant
V_E = -70;          %resting potential 
DET = 2;            %slope of exponential
E_E = 0;            %reversal potential
V_T = -52;          %threshold potential (the spiking threshold for inh neurons)
A_T = 10;           %post spike threshold potential increase
tau_T = 30;         %adaptive threshold time scale
EVthreshold = V_T*ones(1,EneuronNum);           %neuronal threshold vector for all E neurons
EVthresholdA = V_T*ones(1,actionNeuronNum);     %neuronal threshold vector for all read-out neurons
EVthresholdS = V_T*ones(1,superNeuronNum);      %neuronal threshold vector for all supervisor neurons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the I-Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_I = 20;         %Membrane time constant
V_I = -62;          %resting potential
E_I = -75;          %reversal potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the adaptation (Exc only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_w = 100;                %adaptation time constant in the recurrent network
tau_wA = 10;                %adaptation time constant in the read-out
tau_wS = 10;                %adaptation time constant in the supervisor
a = 0;                      %adaptation slope in the recurrent network
aA = 0;                     %adaptation slope of read-out neurons
aS = 0;                     %adaptation slope of supervisor neurons
b = 1000;                   %adaptation amplitude in the recurrent network
bA = 0;                     %adaptation amplitude of read-out neurons
bS = 0;                     %adaptation amplitude of supervisor neurons
w = a*(Vreset-V_E)*ones(1,EneuronNum);              %initial adaptation vector for all recurrent network exc neurons
wA = aA*(Vreset-V_E)*ones(1,actionNeuronNum);       %initial adaptation vector for all read-out neurons
wS = aS*(Vreset-V_E)*ones(1,superNeuronNum);        %initial adaptation vector for all supervisor neurons

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

xedecayA = zeros(1,actionNeuronNum);
xeriseA = zeros(1,actionNeuronNum);
xidecayA = zeros(1,actionNeuronNum);
xiriseA = zeros(1,actionNeuronNum);

xedecayS = zeros(1,superNeuronNum);
xeriseS = zeros(1,superNeuronNum);

xedecayH = zeros(1,interNeuronNum);
xeriseH = zeros(1,interNeuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% External input to both E and I clockneurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rex = 4.50;        %external rate to E-neurons, for synfire chain: 2.75     
rix = 2.25;        %external rate to I-neurons
Jeex = 1.6;         %weights for ee external input
Jiex = 1.52;        %weights for ie external input

nextx = zeros(1,neuronNum);      %vector containing the next external input spike times  for recurrent network neurons
nextx(1,1:EneuronNum) = exprnd(1,1,EneuronNum)/rex;         
nextx(1,1+EneuronNum:end) = exprnd(1,1,IneuronNum)/rix;
rx = zeros(1,neuronNum);
rx(1,1:EneuronNum) = rex;
rx(1,EneuronNum+1:end) = rix;
forwardInputsEPrev = zeros(1,neuronNum);
forwardInputsIPrev = zeros(1,neuronNum);

rhx = 1; %external rate to H-neurons 
jHX = 1.78; %external weight
jAH = 200; %H to A weight (non-plastic)
jHA = 200; %A to H weight (non-plastic)
nextxH = exprnd(1,1,interNeuronNum)/rhx;    %initialize vector containing the next external input spike times for interneurons
forwardInputsHPrev = zeros(1,interNeuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input to read-out neurons (from supervisor neurons and E-neurons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jas = 200;     %weight strenght for connections from the supervisor neurons to read-out neurons
forwardInputsAEPrev = zeros(1,actionNeuronNum); %vector containing excitatory input to the read-outs (rec. netw. + supervisor input)
forwardInputsAIPrev = zeros(1,actionNeuronNum); %vector containing inhibitory input to the read-outs (interneuron input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% External input to supervisor neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jsex = 1.78;     %weights for se external target sequence input
nextxS = exprnd(1,1,superNeuronNum); %initialize next incoming external spike to supervisor
forwardInputsSPrev = zeros(1,superNeuronNum); %vector containing excitatory input to the supervisor neuron (external only)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the synaptic plasticity (vSTDP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th_LTP = -49;               %LTP threshold constant 
th_LTD = -70;               %LTD threshold constant

w_max = 25;                 %maximal read-out weight strength

A_LTP = 0.0008;             %LTP amplitude constant
A_LTD = 0.0014;             %LTD amplitude constant

tau_x = 5;                  %time constant of presynaptic low pass filtered spike train 
tau_u = 10;                 %time constant of postsynaptic low pass filtered membrane voltage (LTD) 
tau_vs = 7;                 %time constant of postsynaptic low pass filtered membrane voltage (LTP) 

x = zeros(1,EneuronNum);                                                 %low pass filtered presynaptic spike train
u = Vreset+(V_T-Vreset)*rand(actionNeuronNum,(tEnd - tStart)/tStep + 1); %low pass filtered postsynaptic membrane voltage (LTD)
vs = Vreset+(V_T-Vreset)*rand(actionNeuronNum,(tEnd - tStart)/tStep + 1);%low pass filtered postsynaptic membrane voltage (LTP)
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Training 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


trainer


%plot of read-out weights
figure;
imagesc(W_AE(end:-1:1,:))
xlabel('Recurrent network neuron index', 'FontSize',25)
ylabel('Read-out neuron', 'FontSize',25)
xticks([1 400 800 1200 1600 2000 2400])
yticks([1 2 3])
yticklabels({'C','B','A'})
colorbar

%generateSequence
