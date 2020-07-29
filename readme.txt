This .txt file contains information regarding the code used to produce the results for: 'Learning spatiotemporal sequences using a recurrent spiking network that discretizes time' (published on 31 January 2020 in Plos cb).
Author: Amadeus Maes
Contact: ahm17@ic.ac.uk or amadeus.maes@gmail.com

The code contains three folders:
1) Learning sequential dynamics
2) Training read-out neurons
3) Data and plotting

1)
In this julia code a recurrent network is trained. Change boolean variables in the file runsim.jl to simulate and train networks. Change variables in sim.jl to change network size and parameters. Repeated sequential input embeds a feedforward structure in the weight matrix. A new matrix can be created and trained or existing matrices can be uploaded. You can turn off the structured stimulation and simulate spontaneous dynamics. Matrices can be saved in .txt format and loaded into the matlab code. The RNNs used in the main results are saved in .mat format in folder 3).

2)
In this matlab code read-out synapses are trained. Change variables in the file clockActionNetwork.m to choose the type of simulation you want to run. RNN weights matrices can be loaded into the file from folder 3).  

3)
Some code to plot and the main trained RNN weight matrices are saved in this folder.


Hopefully the code is readable and can be of help. If anything is unclear or does not work please feel free to contact me.
