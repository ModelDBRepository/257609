%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script does trains the read-out synapses %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rast = zeros(neuronNum,(tEnd - tStart)/tStep + 1);                  %Matrix storing spike times for raster plots
rast_binary = zeros(neuronNum,(tEnd - tStart)/tStep + 1);           %same but with binary numbers
rast_A = zeros(actionNeuronNum,(tEnd - tStart)/tStep + 1);          %storing spike times of read-out neurons (binary)
rast_S = zeros(superNeuronNum,(tEnd - tStart)/tStep + 1);           %storing spike times of supervisor neurons (binary)
rast_H = zeros(interNeuronNum,(tEnd - tStart)/tStep + 1);           %storing spike times of interneurons (binary)

%for refractory period calculation (just big negative number)
lastAP  = -50 * ones(1,neuronNum);      
lastAPA = -50*ones(1,actionNeuronNum);
lastAPS = -50*ones(1,superNeuronNum);
lastAPH = -50*ones(1,interNeuronNum);

%membrane potential
memVol = Vreset+(V_T-Vreset)*rand(neuronNum,(tEnd - tStart)/tStep + 1);
memVolA = Vreset+(V_T-Vreset)*rand(actionNeuronNum,(tEnd - tStart)/tStep + 1);
memVolS = Vreset +(V_T - Vreset)*rand(superNeuronNum,(tEnd - tStart)/tStep + 1);
memVolH = Vreset +(V_T - Vreset)*rand(interNeuronNum,(tEnd - tStart)/tStep + 1);

%for supervisor, variables to keep track of when stimulation begins and finishes
begin = false;
begin_time = [];
finish = false;

for i =2:(tEnd - tStart)/tStep
    
    if mod(i,1000/tStep)==0
        i/10    %print every second elapsed time in ms
    end
    
    forwardInputsE = zeros(1,neuronNum);
    forwardInputsI = zeros(1,neuronNum);
    forwardInputsAE = zeros(1,actionNeuronNum);
    forwardInputsAI = zeros(1,actionNeuronNum);
    forwardInputsS = zeros(1,superNeuronNum);
    forwardInputsH = zeros(1,interNeuronNum);
    
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
            
            x(j) = x(j) - (tStep/tau_x)*x(j);
            
            w(j) = w(j) + (tStep/tau_w)*(a*(memVol(j,i-1) - V_E) - w(j));            %adaptation current            
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
                forwardInputsAE = forwardInputsAE + W_AE(:,j)';
                
                EVthreshold(j) = EVthreshold(j) + A_T;                
                w(j) = w(j) + b;
                x(j) = x(j) + 1;
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
    
    %variables to keep track of learning
    %external input to supervisor neurons is turned on at the beginning of
    %the sequential activity in the recurrent network
    states = zeros(numClusters,1);
    if i>500 %burn-in time 50ms
        for k = 1:numClusters
            states(k) = sum(sum(rast_binary(1+(k-1)*EneuronNum/numClusters:k*EneuronNum/numClusters,i-100:i)));
        end
        temp = find(states==max(states));
        if ~begin && temp(1)==1 
            begin = true;
            finish = false;
            begin_time = [begin_time;i];
        end
        if ~finish && begin && i-begin_time(end)>superLength/tStep
            finish = true;
            begin = false;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   SUPERVISOR NEURONS   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for j = 1:superNeuronNum
   
        if begin && ~finish %add external input to supervisor
            rsex = supervisor(j,i-begin_time(end)+1);
        else
            rsex = 1; %baseline external input to supervisor neurons
        end
     
        while i*tStep>nextxS(j)
            nextxS(j) = nextxS(j) + exprnd(1)/rsex;
            forwardInputsSPrev(j) = forwardInputsSPrev(j) + Jsex;
        end

        xeriseS(j) = xeriseS(j) -tStep*xeriseS(j)/tauerise + forwardInputsSPrev(j);
        xedecayS(j) = xedecayS(j) -tStep*xedecayS(j)/tauedecay + forwardInputsSPrev(j);
        
        gE = (xedecayS(j) - xeriseS(j))/(tauedecay - tauerise);
        
        wS(j) = wS(j) + (tStep/tau_wS)*(aS*(memVolS(j,i-1) - V_E) - wS(j));           %adaptation current
        EVthresholdS(j) = EVthresholdS(j) + (tStep/tau_T)*(V_T - EVthresholdS(j));    %adapting threshold

        %voltage dynamics
        v = memVolS(j,i-1) + (tStep/tau_E)*(-memVolS(j,i-1) + V_E + DET*exp((memVolS(j,i-1)-EVthresholdS(j))/DET)) ...
            + (tStep/C)*( gE*(E_E - memVolS(j,i-1)) - wS(j));                
            
        if ((lastAPS(j) + tau_absS/tStep)>=i)   %Refractory Period
            v = Vreset;
        end

        if (v > Vthres)           %Fire if exceed threshold
            v = Vreset;
            lastAPS(j) = i;
            rast_S(j,i) = j;
            
            forwardInputsAE(j) = forwardInputsAE(j) + Jas; 
            
            wS(j) = wS(j) + bS;
            EVthresholdS(j) = EVthresholdS(j) + A_T;
        end

        memVolS(j,i) = v;                    
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     INTER NEURON      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:interNeuronNum   
    
        while i*tStep>nextxH(j)
            nextxH(j) = nextxH(j) + exprnd(1)/rhx;
            forwardInputsHPrev(j) = forwardInputsHPrev(j) + jHX;
        end
        
        xeriseH(j) = xeriseH(j) -tStep*xeriseH(j)/tauerise + forwardInputsHPrev(j);
        xedecayH(j) = xedecayH(j) -tStep*xedecayH(j)/tauedecay + forwardInputsHPrev(j);
        
        gE = (xedecayH(j) - xeriseH(j))/(tauedecay - tauerise);
        
        %cell dynamics
        v = memVolH(j,i-1) + (tStep/tau_I)*(-memVolH(j,i-1) + V_I) + ...
            (tStep/C)*( gE*(E_E - memVolH(j,i-1)) );
        
        if ((lastAPH + tau_absH/tStep)>=i)     %Refractory Period
            v = Vreset;
        end
        
        if (v > V_T)                  %Fire if exceed threshold
            v = Vreset;
            lastAPH = i;
            rast_H(j,i) = 1;
            
            forwardInputsAI(j) = forwardInputsAI(j) + jAH;
        end
        
        memVolH(j,i) = v;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   READ_OUT NEURONS   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1:actionNeuronNum
                        
        %vSTDP
        u(j) = u(j) + (memVolA(j,i-1) - u(j))*tStep/tau_u;
        vs(j) = vs(j) + (memVolA(j,i-1) - vs(j))*tStep/tau_vs;
                                                
        xeriseA(j) = xeriseA(j) -tStep*xeriseA(j)/tauerise + forwardInputsAEPrev(j);
        xedecayA(j) = xedecayA(j) -tStep*xedecayA(j)/tauedecay + forwardInputsAEPrev(j);
        xidecayA(j) = xidecayA(j) - tStep*xidecayA(j)/tauidecay + forwardInputsAIPrev(j);
        xiriseA(j) = xiriseA(j) - tStep*xiriseA(j)/tauirise + forwardInputsAIPrev(j);
        
        gE = (xedecayA(j) - xeriseA(j))/(tauedecay - tauerise);
        gI = (xidecayA(j) - xiriseA(j))/(tauidecay - tauirise);
        
        wA(j) = wA(j) + (tStep/tau_wA)*(aA*(memVolA(j,i-1) - V_E) - wA(j));                 %adaptation current
        EVthresholdA(j) = EVthresholdA(j) + (tStep/tau_T)*(V_T - EVthresholdA(j));    %adapting threshold

        %voltage dynamics
        v = memVolA(j,i-1) + (tStep/tau_E)*(-memVolA(j,i-1) + V_E + DET*exp((memVolA(j,i-1)-EVthresholdA(j))/DET)) ...
            + (tStep/C)*(gE*(E_E - memVolA(j,i-1)) + gI*(E_I - memVolA(j,i-1)) - wA(j));                
            
        if ((lastAPA(j) + tau_absA/tStep)>=i)   %Refractory Period
            v = Vreset;
        end

        if (v > Vthres)           %Fire if exceed threshold
            v = Vreset;
            lastAPA(j) = i;
            rast_A(j,i) = j;
            
            forwardInputsH(j) = forwardInputsH(j) + jHA; 
            
            wA(j) = wA(j) + bA;
            EVthresholdA(j) = EVthresholdA(j) + A_T;
        end

        memVolA(j,i) = v;            
    
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%%   PLASTICITY   %%%
        %%%%%%%%%%%%%%%%%%%%%%
        A_LTPCORR = A_LTP*(w_max - W_AE(j,:))/w_max;    
        LTP = A_LTPCORR.*x*max(memVolA(j,i)-th_LTP,0)*max(vs(j)-th_LTD,0);
        LTD = A_LTD*rast_binary(1:EneuronNum,i)'*max(u(j)-th_LTD,0);
        W_AE(j,:) = W_AE(j,:) + tStep*(LTP - LTD);
        idx = find(W_AE(j,:)<0); %minimum weight is zero
        W_AE(j,idx) = 0;
     
    end
    
    forwardInputsEPrev = forwardInputsE;
    forwardInputsIPrev = forwardInputsI;
    forwardInputsAEPrev = forwardInputsAE;
    forwardInputsAIPrev = forwardInputsAI;
    forwardInputsSPrev = forwardInputsS;
    forwardInputsHPrev = forwardInputsH;
end

