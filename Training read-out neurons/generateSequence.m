%This script simulates and plots spontaneous dynamics of recurrent network and
%read-outs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = 0;
tEnd = 2000;                      %Simulation in milli-seconds
tStep = 0.1;                      %0.1 millisecond time step

time = [tStart:tStep:tEnd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xedecay = zeros(1,neuronNum);
xerise = zeros(1,neuronNum);
xidecay = zeros(1,neuronNum);
xirise = zeros(1,neuronNum);

xedecayA = zeros(1,actionNeuronNum);
xeriseA = zeros(1,actionNeuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset input vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nextx = zeros(1,neuronNum);      %vector containing the next external input spike times 
nextx(1,1:EneuronNum) = exprnd(1,1,EneuronNum)/rex;         
nextx(1,1+EneuronNum:end) = exprnd(1,1,IneuronNum)/rix;
forwardInputsEPrev = zeros(1,neuronNum);
forwardInputsIPrev = zeros(1,neuronNum);
forwardInputsAPrev = zeros(1,actionNeuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating Network with EIF neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rast = zeros(neuronNum,(tEnd - tStart)/tStep + 1);              %Matrix storing spike times for raster plots
rast_binary = zeros(neuronNum,(tEnd - tStart)/tStep + 1);       %same but with binary numbers
rast_A = zeros(actionNeuronNum,(tEnd - tStart)/tStep + 1);      %storing spike times of action neurons
lastAP  = -50 * ones(1,neuronNum);                              %for refractory period calculation (just big neg. number)
lastAPA = -50*ones(1,actionNeuronNum);


memVol = Vreset+(V_T-Vreset)*rand(neuronNum,(tEnd - tStart)/tStep + 1);
memVolA = Vreset+(V_T-Vreset)*rand(actionNeuronNum,(tEnd - tStart)/tStep + 1);

%for keeping track of sequential activity
begin = false;
begin_time = [];
finish = false;

%synapse deletion
prob = 1/3;

%for plotting 
current = zeros(actionNeuronNum,tEnd/tStep+1);
W_AE_temp = W_AE;
for i =2:(tEnd - tStart)/tStep
    
    if mod(i,1000)==0
        i/10 %print time every 100 ms
    end
    
%     if i==15000
%         tempvar = weightsEE(40,41);
%         weightsEE(40,41) = 0;
%     end

    if i==10000
        for m=1:numClusters
            W_AE(:,1+(m-1)*EneuronNum/numClusters:m*EneuronNum/numClusters) = ...
                W_AE(:,1+(m-1)*EneuronNum/numClusters:m*EneuronNum/numClusters).*(rand(3,EneuronNum/numClusters)>prob);
        end
    end
    
    forwardInputsE = zeros(1,neuronNum);
    forwardInputsI = zeros(1,neuronNum);
    forwardInputsA = zeros(1,actionNeuronNum);
    
    %%%%%%%%%%%%%%%%%%%
    %%%CLOCK NETWORK%%%
    %%%%%%%%%%%%%%%%%%%
   
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
                forwardInputsA = forwardInputsA + W_AE(:,j)';
                
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
                v = Vthres;
                lastAP(j) = i;
                rast(j,i) = j;
                rast_binary(j,i) = 1;
                
                forwardInputsI = forwardInputsI + [weightsEI(:,j-EneuronNum);weightsII(:,j-EneuronNum)]';
            end
            
            memVol(j,i) = v;
        end
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   READ-OUT NEURONS   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1:actionNeuronNum
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CONNCECTIVITY CALCULATIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        xeriseA(j) = xeriseA(j) -tStep*xeriseA(j)/tauerise + forwardInputsAPrev(j);
        xedecayA(j) = xedecayA(j) -tStep*xedecayA(j)/tauedecay + forwardInputsAPrev(j);
        
        gA = (xedecayA(j) - xeriseA(j))/(tauedecay - tauerise);  
        current(j,i) = gA;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DYNAMICS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        states = zeros(numClusters,1);
        if i>500
            for k = 1:numClusters
                states(k) = sum(sum(rast_binary(1+(k-1)*EneuronNum/numClusters:k*EneuronNum/numClusters,i-100:i)));
            end
            temp = find(states==max(states));
            if ~begin && temp(1)==1
                begin = true;
                finish = false;
                begin_time = [begin_time;i];
            end
            if ~finish && begin && i-begin_time(end)>superLength/tStep*0.95
                finish = true;
                begin = false;
            end
        end
        
        wA(j) = wA(j) + (tStep/tau_wA)*(aA*(memVolA(j,i-1) - V_E) - wA(j));           %adaptation current            
        EVthresholdA(j) = EVthresholdA(j) + (tStep/tau_T)*(V_T - EVthresholdA(j));  %adapting threshold

        v = memVolA(j,i-1) + (tStep/tau_E)*(-memVolA(j,i-1) + V_E + DET*exp((memVolA(j,i-1)-EVthresholdA(j))/DET)) ...
                + (tStep/C)*(gA*(E_E - memVolA(j,i-1)) - wA(j));
            
        if ((lastAPA(j) + tau_absA/tStep)>=i)   %Refractory Period
            v = Vreset;
        end

        if (v > Vthres)           %Fire if exceed threshold
            v = Vthres;
            lastAPA(j) = i;
            rast_A(j,i) = 1;
            
            wA(j) = wA(j) + bA;
            EVthresholdA(j) = EVthresholdA(j) + A_T;
        end

        memVolA(j,i) = v;
            
    end
    
    forwardInputsEPrev = forwardInputsE;
    forwardInputsIPrev = forwardInputsI;
    forwardInputsAPrev = forwardInputsA;
    
end

rast_A = rast_A(end:-1:1,:);

%plots of sequential activity and input to action neurons
figure;
ax1 = subplot(2,1,1);
plotRASTER
box(ax1,'off')

for k=1:size(begin_time,1)
    hold on
    line([begin_time(k)/10000,begin_time(k)/10000],[1,EneuronNum])
end

ax2 = subplot(2,1,2);
plotACTIONRASTER
box(ax2,'off')
for k=1:size(begin_time,1)
    hold on
    line([begin_time(k)/10000,begin_time(k)/10000],[0,actionNeuronNum+1])
end

spikes = rast_A(end:-1:1,:);
c = find(spikes>0);
spikes(c) = 1;

kernel = exp(-(-100:0.1:100).^2/250);
rate = spikes;
for i=1:actionNeuronNum
    rate(i,:) = conv(spikes(i,:),kernel,'same');
end

figure;
ax1 = subplot(2,1,1);
imagesc(supervisor(end:-1:1,:))
box(ax1,'off')
set(gca,'FontSize',25)
ylabel({'Target sequence'})
yticks([1 2 3])
yticklabels({'C','B','A'})
xticks([500  1500 2500 3500])
xticklabels({50, 150, 250, 350})

ax2 = subplot(2,1,2);
imagesc(rate(end:-1:1,begin_time(2)+350:begin_time(2)+superLength/tStep+350))
box(ax2,'off')
set(gca,'FontSize',25)
xlabel('Time [ms]');ylabel('Read-out neuron rate')
xticks([500  1500 2500 3500])
xticklabels({50, 150, 250, 350})
yticks([1 2 3])
yticklabels({'C','B','A'})

%The average difference (error) between the target and learned supervisor,
%normalized between [0,1]
T = (supervisor - min(min(supervisor)))./(max(max(supervisor))-min(min(supervisor)));
A = rate(end:-1:1,begin_time(2)+350:begin_time(2)+superLength/tStep+350);
L = (A - min(min(A)))./(max(max(A))-min(min(A)));
Average_error = mean(mean(abs(T-L)))

%restore
% weightsEE(40,41) = tempvar;