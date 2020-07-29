%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting raster plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
rast3 = rast_A(end:-1:1,:);

%Total number of spikes from all neurons in the simulation time
totalNumSpks = find (rast3 > 0);       
spksAll = zeros(length(totalNumSpks),2);
spikeCountTracker = 0;                 

for i = 1:actionNeuronNum    
    spikesNeuron = find (rast3(i,:) > 0);
    numSpikesNeuron = length(spikesNeuron);
    
    if (numSpikesNeuron > 0)
        %to convert to seconds
        spksAll(spikeCountTracker+1:spikeCountTracker+numSpikesNeuron,2) = spikesNeuron'*tStep*(10^-3);    
        spksAll(spikeCountTracker+1:spikeCountTracker+numSpikesNeuron,1) = i*ones(numSpikesNeuron,1);
    end
    if (i <= EneuronNum)
        spksExcStruct(i) = struct('times',spikesNeuron'*tStep*(10^-3));
    else
        spksInhStruct(i) = struct('times',spikesNeuron'*tStep*(10^-3));
    end
    
    spikeCountTracker = spikeCountTracker + numSpikesNeuron;
end


firstInhNeuron = find(spksAll(:,1) >= EneuronNum+1);
if (length(firstInhNeuron) == 0)
    spksExc = spksAll;
    spksInh = [0 0];
else
    spksExc = spksAll(1:firstInhNeuron(1) - 1,:);
    spksInh = spksAll(firstInhNeuron(1):length(spksAll),:);
end


%% Raster Plot of Data
%figure;
set(gcf,'Color',[1 1 1])
plot(spksExc(:,2),spksExc(:,1),'r.'); axis([0 size(rast_A,2)/10000 1 actionNeuronNum])   %convert to seconds
set(gca,'FontSize',25)
xlabel('Time [s]'); ylabel('Read-out neuron')
set(gca,'FontSize',25)
if actionNeuronNum < 4
    yticks([1 2 3 ])
    yticklabels({'A','B','C'})
end
