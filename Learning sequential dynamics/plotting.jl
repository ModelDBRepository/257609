#author: Amadeus Maes
#additional code to make plots

using LinearAlgebra

function plot_spectrum(weights, Ne, Ncells, i)
	println(string("Plot eigenvalue spectrum ",string(i)))
	ylabel("Imag")
	xlabel("Real")
	weights[(1+Ne):Ncells,1:Ne] *= -1
	weights[(1+Ne):Ncells,(1+Ne):Ncells] *= -1
	figure()
	ax = PyPlot.axes()
	ax[:tick_params]("both",labelsize=14)
	ax[:spines]["top"][:set_color]("none") # Remove the top axis boundary
	ax[:spines]["right"][:set_color]("none") # Remove the right axis boundary
	xlabel("Re", size=18)
	ylabel("Im", size=18)
	scatter(real(eigvals(weights)),imag(eigvals(weights)),marker="o",linewidths=0)
	savefig(string("spectra/",string("spectrum_2400_nc80_",string(i)),".png"),dpi=150)
	weights[(1+Ne):Ncells,1:Ne] *= -1
	weights[(1+Ne):Ncells,(1+Ne):Ncells] *= -1
end



function plot_spike_raster(popmembers, times, ns, i)
	println(string("Plot spike raster ",string(i)))
	Npop = size(popmembers,1)
	Nmaxmembers = size(popmembers,2)
	figure()
	ax = PyPlot.axes()
	ax[:tick_params]("both",labelsize=18)
	xlim(0,T)
	ylim(0,sum(popmembers.>0))
	ylabel("Neuron id",size=20)
	xlabel("Time [ms]",size=20)
	tight_layout()

	#plot raster with the order of rows determined by population membership
	rowcount = 0
	for pp = 1:Npop
		for cc = 1:Nmaxmembers
			if popmembers[pp,cc] < 1
				break
			end
			rowcount+=1
			ind = popmembers[pp,cc]
			vals = times[ind,1:ns[ind]]
			y = rowcount*ones(length(vals))
			scatter(vals,y,s=.4,c="r",marker="o",linewidths=0)
		end
	end
	#@printf("\rdone creating plot")
	savefig(string("rasters/",string("raster_2400_nc80_",string(i)),".png"),dpi=150)
end


function plot_exc_weights(popmembers, weights, Ne, i)
	println(string("Plot excitatory weights ",string(i)))
	w = weights[1:Ne,1:Ne]
	exc_weights = zeros(size(w))
	Npop = size(popmembers,1)
	Nmaxmembers = size(popmembers,2)

	#restructure weight matrix to put neurons with same cluster membership together
	prev_count = 0
	others = collect(1:Ne)
	processed = []
	for pp = 1:Npop
		count = 0
		new_group = []
		for cc = 1:Nmaxmembers
			if popmembers[pp,cc]<1
				break
			end
			if !(popmembers[pp,cc] in processed)
				count += 1
				push!(new_group,popmembers[pp,cc])
				push!(processed,popmembers[pp,cc])
			end
		end
		exc_weights[1+prev_count:prev_count+count,1+prev_count:prev_count+count] = w[new_group,new_group]
		for cc = 1:count
			filter!(e->e!=new_group[cc],others)
		end
		exc_weights[1+prev_count+count:Ne,1+prev_count:prev_count+count] = w[others,new_group]
		exc_weights[1+prev_count:prev_count+count,1+prev_count+count:Ne] = w[new_group,others]
		prev_count += count
	end

	exc_weights[1+prev_count:Ne,1+prev_count:Ne] = w[others,others]

	figure()
	ax = PyPlot.axes()
	ax[:tick_params]("both",labelsize=18)
	ax[:spines]["top"][:set_color]("none") # Remove the top axis boundary
	ax[:spines]["right"][:set_color]("none") # Remove the right axis boundary
	xlabel("Post neuron id", size=20)
	ylabel("Pre neuron id",size=20)
	tight_layout()
	matshow(exc_weights,cmap="binary",fignum=0)
	savefig(string("exc_weights/",string("weights_2400_nc80_",string(i)),".png"),dpi=150)
end


