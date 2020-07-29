%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE WEIGHT MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WRatio  = 22.5;             %Ratio of Win/Wout (synaptic weight of within group to neurons outside of the group)
REE = 1.0;                  %Ratio of pin/pout (probability of connection withing group to outside the group)
FF_ratio = 9;               %Ratio of W_up/W_down (above and under diagonal connections)

f = 1/sqrt(mult);           %Factor to scale by synaptic weight parameters by network size


%choose network parameters
wEI     = f*108.9;          %Average weight of inhibitory to excitatory cells
wIE     = f*2.83;           %Average weight of excitatory to inhibitory cells *10 for one neuron clusters
wEE     = f*4.32;           %Average weight of excitatory to excitatory cells
wII     = f*36.22;          %Average weight of inhibitory to inhibitory cells



wEEsub = wEE/(1/numClusters+(1-1/numClusters)/WRatio);          %Average weight for sub-clusters
wEE    = wEEsub/WRatio;

p1  = 0.2/(1/numClusters+(1-1/numClusters)/REE);                %Average probability for sub-clusters
pEE = p1/REE;

p = 0.2;                                                        %connection probability

weightsEI = random('binom',1,p,[EneuronNum,IneuronNum]);        %Weight matrix of inhibitory to excitatory LIF cells
weightsEI = wEI.* weightsEI;

weightsIE = random('binom',1,p,[IneuronNum, EneuronNum]);       %Weight matrix of excitatory to inhibitory cells
weightsIE = wIE.* weightsIE;

weightsII = random('binom',1,p,[IneuronNum, IneuronNum]);     %Weight matrix of inhibitory to inhibitory cells
weightsII = wII.* weightsII;

weightsEE = random('binom',1,pEE,[EneuronNum, EneuronNum]);     %Weight matrix of excitatory to excitatory cells
weightsEE = wEE.* weightsEE;


%Create the group weight matrices and update the total weight matrix
for i = 1:numClusters
    weightsEEsub = random('binom',1,p1,[EneuronNum/numClusters, EneuronNum/numClusters]);
    weightsEEsub = wEEsub.* weightsEEsub;
    weightsEE((i-1)*EneuronNum/numClusters+1:i*EneuronNum/numClusters,(i-1)*EneuronNum/numClusters+1:i*EneuronNum/numClusters) = weightsEEsub;
end

%Create off-diagonal connections in weightsEE
for i = 1:numClusters-1
    weightsEE((i-1)*EneuronNum/numClusters+1:i*EneuronNum/numClusters,i*EneuronNum/numClusters+1:(i+1)*EneuronNum/numClusters) = ...
        weightsEE((i-1)*EneuronNum/numClusters+1:i*EneuronNum/numClusters,i*EneuronNum/numClusters+1:(i+1)*EneuronNum/numClusters)*FF_ratio;
end
weightsEE((numClusters-1)*EneuronNum/numClusters+1:EneuronNum,1:EneuronNum/numClusters) = ...
    weightsEE((numClusters-1)*EneuronNum/numClusters+1:EneuronNum,1:EneuronNum/numClusters)*FF_ratio;


    
%Ensure the diagonals are zero
weightsII = weightsII - diag(diag(weightsII));
weightsEE = weightsEE - diag(diag(weightsEE));
