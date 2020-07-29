#this code was taken and edited from litwin-kumar_doiron_formation_2014
#the original code can be found here: http://lk.zuckermaninstitute.columbia.edu/#code

#edits made by Amadeus Maes 
#contact: amadeus.maes@gmail.com

using PyPlot
using PyCall
using JLD
using DataFrames
using CSV
using LinearAlgebra
using DelimitedFiles
using Distributions

include("sim.jl")
include("plotting.jl")

#if true, loads trained network.  if false, generates a new network.
loadtrained = false
#if true, sets up assemblies that may overlap. A group of neurons will not belong to any assembly.
#if false, sets up disjoint assemblies. All neurons belong to exactly one assembly.
overlapping = false


#stimulus matrix is in the following format:
#column 1: index of stimulated population
#columns 2/3: start/stop time of stimulus
#column 4: rate of external drive (in kHz)

#FF (2 minutes/iteration)
stim = zeros(8000,4)
for i=1:8000
	stim[i,1] = i % 30 + 1
	stim[i,2] = 1+15*i
	stim[i,3] = 10+15*i
	stim[i,4] = 18
end


no_stim = zeros(1,4)
no_stim[1,1] = 1


if loadtrained
	#load FF weight matrix
	weights_FF = load("weight_matrix_2400_nc80a.jld")["weights"]
	#load SSA weight matrix
	weights = load("weight_matrix_2400_nc80a.jld")["weights"]
	#load popmembers matrix
	popmembers = load("popmembers_matrix_2400_nc80.jld")["popmembers"]

	#spontaneous dynamics
	structure1_reinf = zeros(31,30);
	structure2_reinf = zeros(31,30);
	for j=1:29
		structure1_reinf[1,j]=mean(weights_FF[1+80*(j-1):80*j,1+80*j:80*(j+1)]);
		structure2_reinf[1,j]=mean(weights_FF[1+80*j:80*(j+1),1+80*(j-1):80*j]);
	end
	structure1_reinf[1,30]=mean(weights_FF[2321:2400,1:80]);
	structure2_reinf[1,30]=mean(weights_FF[1:80,2321:2400]);


	save("structure1_reinf.jld","structure1_reinf",structure1_reinf)
	save("structure2_reinf.jld","structure2_reinf",structure2_reinf)
	for i=1:30
		global times,ns,Ne,Ncells, weights, weights_FF, T

		times,ns,Ne,Ncells, weights, weights_FF, T, _, _ = sim_degradation(no_stim,weights, weights_FF,popmembers, 120000, 1000)
		#plot_exc_weights(popmembers, weights_FF, Ne, (i-1)%30 +1)
		#plot_spectrum(weights_FF, Ne, Ncells, (i-1)%30 +1)
		#5 sec recording of spontaneous activity of neurons
		#times,ns,Ne,Ncells, weights, T, _, _ = sim(no_stim,weights,popmembers, 3000, 3000)
		#plot_spike_raster(popmembers, times, ns, (i-1)%30 +1)

		for j=1:29
			structure1_reinf[i+1,j]=mean(weights_FF[1+80*(j-1):80*j,1+80*j:80*(j+1)]);
			structure2_reinf[i+1,j]=mean(weights_FF[1+80*j:80*(j+1),1+80*(j-1):80*j]);
		end
		structure1_reinf[i+1,30]=mean(weights_FF[2321:2400,1:80]);
		structure2_reinf[i+1,30]=mean(weights_FF[1:80,2321:2400]);


		save("structure1_reinf.jld","structure1_reinf",structure1_reinf)
		save("structure2_reinf.jld","structure2_reinf",structure2_reinf)

		println("Save weight matrix")
		save("weight_matrix_reinforcement.jld","weights",weights_FF)
	end


else
	#set up balanced network
	times,ns,popmembers,Ne,Ncells, weights, T = simnew(no_stim, 5000, 1000, overlapping)
	#save initial matrix
	#save("weights_init.jld","weights", weights)
	#save cluster organization
	save("popmembers_matrix_2400_nc80.jld","popmembers", popmembers)

	#initial plots
	plot_exc_weights(popmembers, weights, Ne, 0)
	#plot_spectrum(weights, Ne, Ncells, 0)

	times,ns,Ne,Ncells, weights, T, _, _ = sim(stim,weights,popmembers, 3000, 0)
	plot_spike_raster(popmembers, times, ns, 0)

	#train the network
	for i=1:30
		global times,ns,Ne,Ncells, weights, T

		#create stimuli by sampling from transition matrix
		times,ns,Ne,Ncells, weights, T, _, _ = sim(stim,weights,popmembers, 120000, 0)
		plot_exc_weights(popmembers, weights, Ne, (i-1)%30 +1)
		plot_spectrum(weights, Ne, Ncells, (i-1)%30 +1)
		#5 sec recording of spontaneous activity of neurons
		times,ns,Ne,Ncells, weights, T, _, _ = sim(no_stim,weights,popmembers, 2000, 2000)
		plot_spike_raster(popmembers, times, ns, (i-1)%30 +1)

		println("Save weight matrix")
		save("weight_matrix_2400_nc80a.jld","weights",weights)
	end

	writedlm("FF60a_2400_nc80.txt",weights)


	#simulate spontaneous dynamics
	for i=1:30
		global times,ns,Ne,Ncells, weights, T

		times,ns,Ne,Ncells, weights, T, _, _ = sim(no_stim,weights,popmembers, 120000, 1000)
		plot_exc_weights(popmembers, weights, Ne, (i-1)%30 +1)
		plot_spectrum(weights, Ne, Ncells, (i-1)%30 +1)
		#5 sec recording of spontaneous activity of neurons
		times,ns,Ne,Ncells, weights, T, _, _ = sim(no_stim,weights,popmembers, 3000, 3000)
		plot_spike_raster(popmembers, times, ns, (i-1)%30 +1)

		println("Save weight matrix")
		save("weight_matrix_2400_nc80b.jld","weights",weights)
	end

	writedlm("FF60b_2400_nc80.txt",weights)


	#simulate trained network for 5 seconds and plot spike raster
	times,ns,Ne,Ncells, weights, T, _, _ = sim(no_stim,weights,popmembers, 5000, 1000)
	plot_spike_raster(popmembers, times, ns, 51)

end


println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):Ncells]/T)," Hz")
