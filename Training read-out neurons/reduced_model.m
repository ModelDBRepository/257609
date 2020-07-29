%reduced weight matrix by averaging over clusters
%this justifies more explicitly the linear intuition

A = load('FF60b_2400_nc80.mat'); %import weight matrix
A = A.A;
neuronNum = size(A,1);
EneuronNum = 0.80*size(A,1);
IneuronNum = 0.20*size(A,1);
numClusters = EneuronNum/80; %80 neurons/cluster

weightsEE = A(1:EneuronNum,1:EneuronNum);
weightsEI = A(1:EneuronNum,EneuronNum+1:neuronNum);
weightsIE = A(EneuronNum+1:neuronNum,1:EneuronNum);
weightsII = A(EneuronNum+1:neuronNum,EneuronNum+1:neuronNum);

reduced_weightsEE = zeros(numClusters,numClusters);
for i=1:numClusters
    for j=1:numClusters
        reduced_weightsEE(i,j) = mean(mean(weightsEE(1+(i-1)*80:i*80,1+(j-1)*80:j*80)));
    end
end

reduced_weightsEI = zeros(numClusters,10);
for i=1:numClusters
    for j=1:10
        reduced_weightsEI(i,j) = mean(mean(weightsEI(1+(i-1)*80:i*80,1+(j-1)*60:j*60)));
    end
end

reduced_weightsIE = zeros(10,numClusters);
for i=1:numClusters
    for j=1:10
        reduced_weightsIE(j,i) = mean(mean(weightsIE(1+(j-1)*60:j*60,1+(i-1)*80:i*80)));
    end
end

reduced_weightsII = zeros(10,10);
for i=1:10
    for j=1:10
        reduced_weightsII(j,i) = mean(mean(weightsII(1+(j-1)*60:j*60,1+(i-1)*60:i*60)));
    end
end
[lambdas, W]= get_eigenvalues_LIF(reduced_weightsEE,reduced_weightsIE,reduced_weightsEI,reduced_weightsII);
figure; plot(real(lambdas),imag(lambdas),'.','color',[0 0.3 0], 'MarkerSize',24)

