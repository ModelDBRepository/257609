#this code was taken and edited from litwin-kumar_doiron_formation_2014
#the original code can be found here: http://lk.zuckermaninstitute.columbia.edu/#code

#edits made by Amadeus Maes 
#contact: amadeus.maes@gmail.com

using Distributions
using Printf

function simnew(stim, T, stdpdelay, overlapping) #generates new weights and populations with unpotentiated synapses, runs simulation
	println("setting up weights")

	Ne,Ni,jee0,jei0,jie,jii,p = weightpars()
	Ncells = Ne+Ni

	#set up weights
	#note: weights are set up so that w[i,j] is weight from presynaptic i to postsynaptic j
	#this is for performance: iterating over presynaptic indices is more important and
	#Julia uses column-major arrays
	weights = zeros(Ncells,Ncells)
	weights[1:Ne,1:Ne] .= jee0
	weights[1:Ne,(1+Ne):Ncells] .= jie
	weights[(1+Ne):Ncells,1:Ne] .= jei0
	weights[(1+Ne):Ncells,(1+Ne):Ncells] .= jii

	weights = weights.*(rand(Ncells,Ncells) .< p)
	#weights=load("weights_init.jld")["weights"];

	for cc = 1:Ncells
		weights[cc,cc] = 0
	end

	#populations
	Npop = 30 #number of assemblies
	pmembership = .1 #probability of belonging to any assembly, if overlapping==true
	Nmaxmembers = 210 #maximum number of neurons in a population (to set size of matrix)

	#set up populations
	popmembers = zeros(Int,Npop,Nmaxmembers)
	if overlapping
		for pp = 1:Npop
			members = find(rand(Ne) .< pmembership)
			popmembers[pp,1:length(members)] = members
			println(length(members))
		end
	else

		for pp = 1:Npop
			members = collect(1+Ne*(pp-1)/Npop:Ne*pp/Npop)
			popmembers[pp,1:length(members)] = members
		end

	end

	times,ns,Ne,Ncells, weights, T, _, _ = sim(stim,weights,popmembers, T, stdpdelay)
	return times,ns,popmembers,Ne,Ncells, weights, T
end


function sim(stim,weights,popmembers, T, stdpdelay) #runs simulation given weight matrix and populations
	println("setting up parameters")

	Ne,Ni,jee0,jei0,jie,jii,p = weightpars()

	#membrane dynamics
	taue = 20 #e membrane time constant
	taui = 20 #i membrane time constant
	vleake = -70 #e resting potential
	vleaki = -62 #i resting potential
	deltathe = 2 #eif slope parameter
	C = 300 #capacitance
	erev = 0 #e synapse reversal potential
	irev = -75 #i synapse reversal potntial
	vth0 = -52 #initial spike voltage threshold
	ath = 10 #increase in threshold post spike
	tauth = 30 #threshold decay timescale
	vre = -60 #reset potential
	taurefrac = 5 #absolute refractory period
	aw_adapt = 0 #adaptation parameter a
	bw_adapt = 1000 #adaptation parameter b
	tauw_adapt = 100#150 #adaptation timescale

	#connectivity
	Ncells = Ne+Ni
	tauerise = 1 #e synapse rise time
	tauedecay = 6 #e synapse decay time
	tauirise = .5 #i synapse rise time
	tauidecay = 2 #i synapse decay time
	rex = 4.5 #external input rate to e (khz)
	rix = 2.25 #external input rate to i (khz)

	jeemin = 1.45  #minimum ee strength
	jeemax = 32.68 #maximum ee strength

	jeimin = 48.7#40 #minimum ei strength
	jeimax = 243#200 #maximum ei strength

	jex = 1.6	#external to e strength
	jix = 1.52	#external to i strength

	#voltage based stdp
	altd = 0.0014
	altp = 0.0008 #factor ten bigger than original paper clopath et al.
	thetaltd = -70 #ltd voltage threshold
	thetaltp = -49 #ltp voltage threshold
	tauu = 10 #timescale for u variable
	tauv = 7 #timescale for v variable
	taux = 3.5 #timescale for x variable

	#inhibitory stdp
	tauy = 20 #width of istdp curve
	eta = 1 #istdp learning rate
	r0 = .003 #target rate (khz)

	#populations
	Npop = size(popmembers,1) #number of assemblies
	Nmaxmembers = size(popmembers,2) #maximum number of neurons in a population

	#simulation
	dt = .1 #integration timestep
	#T = 4000 #simulation time
	Nskip = 1000 #how often (in number of timesteps) to save w_in
	vpeak = 20 #cutoff for voltage.  when crossed, record a spike and reset
	dtnormalize = 20 #how often to normalize columns of ee weights, how does this constant influence learning?
	#stdpdelay = 1000 #time before stdp is activated, to allow transients to die out
	Nspikes = 5000 #maximum number of spikes to record per neuron

	#weight strengths, measured every second
	relative_strengths = zeros(Int(T/1000),Npop)
	#average conductances to 2 clusters, measured every ms
	# first size(popmembers,1) columns contain gE of clusters, 2nd contains gI of clusters
	# third contains external input of clusters
	avg_conductances = zeros(T,size(popmembers,1)*3)

	times = zeros(Ncells,Nspikes)
	ns = zeros(Int,Ncells)

	forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
	no_noise = zeros(Ncells)
	only_noise = zeros(Ncells)
	forwardInputsI = zeros(Ncells)
	forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
	no_noisePrev = zeros(Ncells)
	only_noisePrev = zeros(Ncells)
	forwardInputsIPrev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xerise_noise = zeros(Ncells)
	xedecay_noise = zeros(Ncells)
	xerise_no_noise = zeros(Ncells)
	xedecay_no_noise = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)

	expdist = Exponential()

	v = zeros(Ncells) #membrane voltage
	nextx = zeros(Ncells) #time of next external excitatory input
	sumwee0 = zeros(Ne) #initial summed e weight, for normalization
	Nee = zeros(Int,Ne) #number of e->e inputs, for normalization
	rx = zeros(Ncells) #rate of external input
	for cc = 1:Ncells
		v[cc] = vre + (vth0-vre)*rand()
		if cc < Ne
			rx[cc] = rex
			nextx[cc] = rand(expdist)/rx[cc]
			for dd = 1:Ne
				sumwee0[cc] += weights[dd,cc]#^1.5 #L2 normalization
				if weights[dd,cc] > 0
					Nee[cc] += 1
				end
			end
		else
			rx[cc] = rix
			nextx[cc] = rand(expdist)/rx[cc]
		end
	end

	vth = vth0*ones(Ncells) #adaptive threshold
	wadapt = aw_adapt*(vre-vleake)*ones(Ne) #adaptation current
	lastSpike = -100*ones(Ncells) #last time the neuron spiked
	trace_istdp = zeros(Ncells) #low-pass filtered spike train for istdp
	u_vstdp = vre*zeros(Ne)
	v_vstdp = vre*zeros(Ne)
	x_vstdp = zeros(Ne)

	Nsteps = round(Int,T/dt)
	inormalize = round(Int,dtnormalize/dt)

	println("starting simulation")

	#begin main simulation loop
	for tt = 1:Nsteps
		if mod(tt,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*tt/Nsteps))
		end

		if mod(tt,10000) == 0  #measure W_in/W_out every second
			relative_strengths[Int(tt/10000),:] = rel_weight_strengths(weights, popmembers, Ne)
		end


		t = dt*tt
		forwardInputsE[:] .= 0.
		forwardInputsI[:] .= 0.
		no_noise[:] .= 0

		#check if we have entered or exited a stimulation period
		tprev = dt*(tt-1)
		for ss = 1:size(stim)[1]
			if (tprev<stim[ss,2]) && (t>=stim[ss,2])  #just entered stimulation period
				ipop = round(Int,stim[ss,1])
				for ii = 1:Nmaxmembers
					if (popmembers[ipop,ii] == -1) || (popmembers[ipop,ii] == 0)
						break
					end
					rx[popmembers[ipop,ii]] += stim[ss,4]
				end
			end

			if (tprev<stim[ss,3]) && (t>=stim[ss,3]) #just exited stimulation period
				ipop = round(Int,stim[ss,1])
				for ii = 1:Nmaxmembers
					if (popmembers[ipop,ii] == -1) || (popmembers[ipop,ii] == 0)
						break
					end
					rx[popmembers[ipop,ii]] -= stim[ss,4]
				end
			end
		end #end loop over stimuli


		if mod(tt,inormalize) == 0 #excitatory synaptic normalization
			for cc = 1:Ne
				sumwee = 0.
				for dd = 1:Ne
					sumwee += weights[dd,cc]#^1.5 #L2 normalization
				end

				for dd = 1:Ne
					if weights[dd,cc] > 0.
						weights[dd,cc] -= (sumwee-sumwee0[cc])/Nee[cc]#*0.000143 #L1 soft normalization dt=1ms tau=70s creates some in betweeen values
						if weights[dd,cc] < jeemin
							weights[dd,cc] = jeemin
						elseif weights[dd,cc] > jeemax
							weights[dd,cc] = jeemax
						end
					end
				end
			end
		end #end normalization


		#update single cells
		spiked = zeros(Bool,Ncells)
		for cc = 1:Ncells
			trace_istdp[cc] -= dt*trace_istdp[cc]/tauy

			while(t > nextx[cc]) #external input
				nextx[cc] += rand(expdist)/rx[cc]
				if cc < Ne
					only_noisePrev[cc] += jex
					forwardInputsEPrev[cc] += jex
					if stdpdelay == 0 && rx[cc] == rex
						forwardInputsEPrev[cc] -= 1.5*jex
					end
				else
					only_noisePrev[cc] += jix
					forwardInputsEPrev[cc] += jix
				end
			end

			#both noise and recurrent input
			xerise[cc] += -dt*xerise[cc]/tauerise + forwardInputsEPrev[cc]
			xedecay[cc] += -dt*xedecay[cc]/tauedecay + forwardInputsEPrev[cc]
			xirise[cc] += -dt*xirise[cc]/tauirise + forwardInputsIPrev[cc]
			xidecay[cc] += -dt*xidecay[cc]/tauidecay + forwardInputsIPrev[cc]

			#only noise
			xerise_noise[cc] += -dt*xerise_noise[cc]/tauerise + only_noisePrev[cc]
			xedecay_noise[cc] += -dt*xedecay_noise[cc]/tauedecay + only_noisePrev[cc]

			#no noise
			xerise_no_noise[cc] += -dt*xerise_no_noise[cc]/tauerise + no_noisePrev[cc]
			xedecay_no_noise[cc] += -dt*xedecay_no_noise[cc]/tauedecay + no_noisePrev[cc]


			if cc < Ne
				vth[cc] += dt*(vth0 - vth[cc])/tauth;
				wadapt[cc] += dt*(aw_adapt*(v[cc]-vleake) - wadapt[cc])/tauw_adapt;
				u_vstdp[cc] += dt*(v[cc] - u_vstdp[cc])/tauu;
				v_vstdp[cc] += dt*(v[cc] - v_vstdp[cc])/tauv;
				x_vstdp[cc] -= dt*x_vstdp[cc]/taux;
			end

			# update membrane voltage
			#total conductances
			ge = (xedecay[cc] - xerise[cc])/(tauedecay - tauerise);
			gi = (xidecay[cc] - xirise[cc])/(tauidecay - tauirise);

			#only noise
			ge_noise = (xedecay_noise[cc] - xerise_noise[cc])/(tauedecay - tauerise);

			#no noise
			ge_no_noise = (xedecay_no_noise[cc] - xerise_no_noise[cc])/(tauedecay - tauerise);


			#not in refractory period, in refr period neuron is clamped to vre
			if t > (lastSpike[cc] + taurefrac)


				if cc < Ne #excitatory neuron (eif), has adaptation
					dv = (vleake - v[cc] + deltathe*exp((v[cc]-vth[cc])/deltathe))/taue + ge*(erev-v[cc])/C + gi*(irev-v[cc])/C - wadapt[cc]/C;
					v[cc] += dt*dv;
					if v[cc] > vpeak
						spiked[cc] = true
						wadapt[cc] += bw_adapt
					end
				else
					dv = (vleaki - v[cc])/taui + ge*(erev-v[cc])/C + gi*(irev-v[cc])/C;
					v[cc] += dt*dv;
					if v[cc] > vth0
						spiked[cc] = true
					end
				end

				if spiked[cc] #spike occurred
					spiked[cc] = true;
					v[cc] = vre;
					lastSpike[cc] = t;
					ns[cc] += 1;
					if ns[cc] <= Nspikes
						times[cc,ns[cc]] = t;
					end
					trace_istdp[cc] += 1.;
					if cc<Ne
						x_vstdp[cc] += 1.;
					end

					if cc < Ne
						vth[cc] = vth0 + ath;
					end

					#loop over synaptic projections
					for dd = 1:Ncells
						if cc < Ne #excitatory synapse
							no_noise[dd] += weights[cc,dd];
							forwardInputsE[dd] += weights[cc,dd];
						else #inhibitory synapse
							forwardInputsI[dd] += weights[cc,dd];
						end
					end

				end #end if(spiked)
			end #end if(not refractory)

			#istdp
			if spiked[cc] && (t > stdpdelay)
				if cc <= Ne #excitatory neuron fired, potentiate i inputs
					for dd = (Ne+1):Ncells
						if weights[dd,cc] == 0.
							continue
						end
						weights[dd,cc] += eta*trace_istdp[dd]
						if weights[dd,cc] > jeimax
							weights[dd,cc] = jeimax
						end
					end
				else #inhibitory neuron fired, modify outputs to e neurons
					for dd = 1:Ne
						if weights[cc,dd] == 0.
							continue
						end
						weights[cc,dd] += eta*(trace_istdp[dd] - 2*r0*tauy)
						if weights[cc,dd] > jeimax
							weights[cc,dd] = jeimax
						elseif weights[cc,dd] < jeimin
							weights[cc,dd] = jeimin
						end
					end
				end
			end #end istdp


			#vstdp, ltd component
			if spiked[cc] && (t > stdpdelay) && (cc < Ne)
				for dd = 1:Ne #depress weights from cc to dd
					if weights[cc,dd] == 0.
						continue
					end

					if u_vstdp[dd] > thetaltd
						weights[cc,dd] -= dt*altd*(u_vstdp[dd]-thetaltd) #changed dt*
						if weights[cc,dd] < jeemin
							weights[cc,dd] = jeemin

						end
					end
				end
			end #end ltd

			#vstdp, ltp component
			if (t > stdpdelay) && (cc < Ne) && (v[cc] > thetaltp) && (v_vstdp[cc] > thetaltd)
				for dd = 1:Ne #potentiate weights from dd to cc
					if weights[dd,cc] == 0.
						continue
					end

					weights[dd,cc] += dt*altp*x_vstdp[dd]*(v[cc] - thetaltp)*(v_vstdp[cc] - thetaltd); #*((jeemax-weights[dd,cc])/(jeemax-jeemin))#weight dependent
					if weights[dd,cc] > jeemax
						weights[dd,cc] = jeemax
					end
				end
			end #end ltp

		end #end loop over cells
		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
		no_noisePrev = copy(no_noise)
		only_noisePrev = copy(only_noise)
	end #end loop over time
	@printf("\r")

	#times = times[:,1:maximum(ns)]


	return times,ns,Ne,Ncells, weights, T, relative_strengths, avg_conductances
end

function weightpars() #parameters needed to generate weight matrix
	mult = 1#8/3
	f = 1/sqrt(mult)
	Ne = Int(mult*2400)
	Ni = Int(mult*600)
	jee0 = f*2.83   #initial ee strength
	jei0 = f*62.87   #initial ei strength
	jie =  f*1.96 #ie strength (not plastic)
	jii =  f*20.91 #ii strength (not plastic)
	p = 0.2
	return Ne,Ni,jee0,jei0,jie,jii,p
end


function sim_degradation(stim,weights, weights_FF,popmembers, T, stdpdelay) #runs simulation given weight matrix and populations
	println("setting up parameters")

	Ne,Ni,jee0,jei0,jie,jii,p = weightpars()

	#membrane dynamics
	taue = 20 #e membrane time constant
	taui = 20 #i membrane time constant
	vleake = -70 #e resting potential
	vleaki = -62 #i resting potential
	deltathe = 2 #eif slope parameter
	C = 300 #capacitance
	erev = 0 #e synapse reversal potential
	irev = -75 #i synapse reversal potntial
	vth0 = -52 #initial spike voltage threshold
	ath = 10 #increase in threshold post spike
	tauth = 30 #threshold decay timescale
	vre = -60 #reset potential
	taurefrac = 5 #absolute refractory period
	aw_adapt = 0 #adaptation parameter a
	bw_adapt = 1000 #adaptation parameter b
	tauw_adapt = 100#150 #adaptation timescale

	#connectivity
	Ncells = Ne+Ni
	tauerise = 1 #e synapse rise time
	tauedecay = 6 #e synapse decay time
	tauirise = .5 #i synapse rise time
	tauidecay = 2 #i synapse decay time
	rex = 4.5 #external input rate to e (khz)
	rix = 2.25 #external input rate to i (khz)

	jeemin = 1.45  #minimum ee strength
	jeemax = 32.68 #maximum ee strength

	jeimin = 48.7#40 #minimum ei strength
	jeimax = 243#200 #maximum ei strength

	jex = 1.6	#external to e strength
	jix = 1.52	#external to i strength

	#voltage based stdp
	altd = 0.0014
	altp = 0.0008 #factor ten bigger than original paper clopath et al.
	thetaltd = -70 #ltd voltage threshold
	thetaltp = -49 #ltp voltage threshold
	tauu = 10 #timescale for u variable
	tauv = 7 #timescale for v variable
	taux = 3.5 #timescale for x variable

	#inhibitory stdp
	tauy = 20 #width of istdp curve
	eta = 1 #istdp learning rate
	r0 = .003 #target rate (khz)

	#populations
	Npop = size(popmembers,1) #number of assemblies
	Nmaxmembers = size(popmembers,2) #maximum number of neurons in a population

	#simulation
	dt = .1 #integration timestep
	#T = 4000 #simulation time
	Nskip = 1000 #how often (in number of timesteps) to save w_in
	vpeak = 20 #cutoff for voltage.  when crossed, record a spike and reset
	dtnormalize = 20 #how often to normalize columns of ee weights, how does this constant influence learning?
	#stdpdelay = 1000 #time before stdp is activated, to allow transients to die out
	Nspikes = 5000 #maximum number of spikes to record per neuron

	#weight strengths, measured every second
	relative_strengths = zeros(Int(T/1000),Npop)
	#average conductances to 2 clusters, measured every ms
	# first size(popmembers,1) columns contain gE of clusters, 2nd contains gI of clusters
	# third contains external input of clusters
	avg_conductances = zeros(T,size(popmembers,1)*3)


	times = zeros(Ncells,Nspikes)
	ns = zeros(Int,Ncells)

	forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
	no_noise = zeros(Ncells)
	only_noise = zeros(Ncells)
	forwardInputsI = zeros(Ncells)
	forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
	no_noisePrev = zeros(Ncells)
	only_noisePrev = zeros(Ncells)
	forwardInputsIPrev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xerise_noise = zeros(Ncells)
	xedecay_noise = zeros(Ncells)
	xerise_no_noise = zeros(Ncells)
	xedecay_no_noise = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)

	expdist = Exponential()

	v = zeros(Ncells) #membrane voltage
	nextx = zeros(Ncells) #time of next external excitatory input
	sumwee0 = zeros(Ne) #initial summed e weight, for normalization
	Nee = zeros(Int,Ne) #number of e->e inputs, for normalization
	rx = zeros(Ncells) #rate of external input
	for cc = 1:Ncells
		v[cc] = vre + (vth0-vre)*rand()
		if cc < Ne
			rx[cc] = rex
			nextx[cc] = rand(expdist)/rx[cc]
			for dd = 1:Ne
				sumwee0[cc] += weights_FF[dd,cc] #L1 normalization
				if weights_FF[dd,cc] > 0
					Nee[cc] += 1
				end
			end
		else
			rx[cc] = rix
			nextx[cc] = rand(expdist)/rx[cc]
		end
	end

	vth = vth0*ones(Ncells) #adaptive threshold
	wadapt = aw_adapt*(vre-vleake)*ones(Ne) #adaptation current
	lastSpike = -100*ones(Ncells) #last time the neuron spiked
	trace_istdp = zeros(Ncells) #low-pass filtered spike train for istdp
	u_vstdp = vre*zeros(Ne)
	v_vstdp = vre*zeros(Ne)
	x_vstdp = zeros(Ne)

	Nsteps = round(Int,T/dt)
	inormalize = round(Int,dtnormalize/dt)

	println("starting simulation")

	#begin main simulation loop
	for tt = 1:Nsteps
		if mod(tt,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*tt/Nsteps))
		end

		if mod(tt,10000) == 0  #measure W_in/W_out every second
			relative_strengths[Int(tt/10000),:] = rel_weight_strengths(weights, popmembers, Ne)
		end


		t = dt*tt
		forwardInputsE[:] .= 0.
		forwardInputsI[:] .= 0.
		no_noise[:] .= 0

		#check if we have entered or exited a stimulation period
		tprev = dt*(tt-1)
		for ss = 1:size(stim)[1]
			if (tprev<stim[ss,2]) && (t>=stim[ss,2])  #just entered stimulation period
				ipop = round(Int,stim[ss,1])
				for ii = 1:Nmaxmembers
					if (popmembers[ipop,ii] == -1) || (popmembers[ipop,ii] == 0)
						break
					end
					rx[popmembers[ipop,ii]] += stim[ss,4]
				end
			end

			if (tprev<stim[ss,3]) && (t>=stim[ss,3]) #just exited stimulation period
				ipop = round(Int,stim[ss,1])
				for ii = 1:Nmaxmembers
					if (popmembers[ipop,ii] == -1) || (popmembers[ipop,ii] == 0)
						break
					end
					rx[popmembers[ipop,ii]] -= stim[ss,4]
				end
			end
		end #end loop over stimuli


		if mod(tt,inormalize) == 0 #excitatory synaptic normalization
			for cc = 1:Ne
				sumwee = 0.
				for dd = 1:Ne
					sumwee += weights_FF[dd,cc] #L1 normalization
				end

				for dd = 1:Ne
					if weights_FF[dd,cc] > 0.
						weights_FF[dd,cc] -= (sumwee-sumwee0[cc])/Nee[cc] #L1 soft normalization dt=1ms tau=70s creates some in betweeen values
						if weights_FF[dd,cc] < jeemin
							weights_FF[dd,cc] = jeemin
						elseif weights_FF[dd,cc] > jeemax
							weights_FF[dd,cc] = jeemax
						end
					end
				end
			end
		end #end normalization


		#update single cells
		spiked = zeros(Bool,Ncells)
		for cc = 1:Ncells
			trace_istdp[cc] -= dt*trace_istdp[cc]/tauy

			while(t > nextx[cc]) #external input
				nextx[cc] += rand(expdist)/rx[cc]
				if cc < Ne
					only_noisePrev[cc] += jex
					forwardInputsEPrev[cc] += jex
					if stdpdelay == 0 && rx[cc] == rex
						forwardInputsEPrev[cc] -= 1.5*jex
					end
				else
					only_noisePrev[cc] += jix
					forwardInputsEPrev[cc] += jix
				end
			end

			#both noise and recurrent input
			xerise[cc] += -dt*xerise[cc]/tauerise + forwardInputsEPrev[cc]
			xedecay[cc] += -dt*xedecay[cc]/tauedecay + forwardInputsEPrev[cc]
			xirise[cc] += -dt*xirise[cc]/tauirise + forwardInputsIPrev[cc]
			xidecay[cc] += -dt*xidecay[cc]/tauidecay + forwardInputsIPrev[cc]

			#only noise
			xerise_noise[cc] += -dt*xerise_noise[cc]/tauerise + only_noisePrev[cc]
			xedecay_noise[cc] += -dt*xedecay_noise[cc]/tauedecay + only_noisePrev[cc]

			#no noise
			xerise_no_noise[cc] += -dt*xerise_no_noise[cc]/tauerise + no_noisePrev[cc]
			xedecay_no_noise[cc] += -dt*xedecay_no_noise[cc]/tauedecay + no_noisePrev[cc]


			if cc < Ne
				vth[cc] += dt*(vth0 - vth[cc])/tauth;
				wadapt[cc] += dt*(aw_adapt*(v[cc]-vleake) - wadapt[cc])/tauw_adapt;
				u_vstdp[cc] += dt*(v[cc] - u_vstdp[cc])/tauu;
				v_vstdp[cc] += dt*(v[cc] - v_vstdp[cc])/tauv;
				x_vstdp[cc] -= dt*x_vstdp[cc]/taux;
			end

			# update membrane voltage
			#total conductances
			ge = (xedecay[cc] - xerise[cc])/(tauedecay - tauerise);
			gi = (xidecay[cc] - xirise[cc])/(tauidecay - tauirise);

			#only noise
			ge_noise = (xedecay_noise[cc] - xerise_noise[cc])/(tauedecay - tauerise);

			#no noise
			ge_no_noise = (xedecay_no_noise[cc] - xerise_no_noise[cc])/(tauedecay - tauerise);


			#not in refractory period, in refr period neuron is clamped to vre
			if t > (lastSpike[cc] + taurefrac)


				if cc < Ne #excitatory neuron (eif), has adaptation
					dv = (vleake - v[cc] + deltathe*exp((v[cc]-vth[cc])/deltathe))/taue + ge*(erev-v[cc])/C + gi*(irev-v[cc])/C - wadapt[cc]/C;
					v[cc] += dt*dv;
					if v[cc] > vpeak
						spiked[cc] = true
						wadapt[cc] += bw_adapt
					end
				else
					dv = (vleaki - v[cc])/taui + ge*(erev-v[cc])/C + gi*(irev-v[cc])/C;
					v[cc] += dt*dv;
					if v[cc] > vth0
						spiked[cc] = true
					end
				end

				if spiked[cc] #spike occurred
					spiked[cc] = true;
					v[cc] = vre;
					lastSpike[cc] = t;
					ns[cc] += 1;
					if ns[cc] <= Nspikes
						times[cc,ns[cc]] = t;
					end
					trace_istdp[cc] += 1.;
					if cc<Ne
						x_vstdp[cc] += 1.;
					end

					if cc < Ne
						vth[cc] = vth0 + ath;
					end

					#loop over synaptic projections
					for dd = 1:Ncells
						if cc < Ne #excitatory synapse
							no_noise[dd] += weights[cc,dd];
							forwardInputsE[dd] += weights[cc,dd];
						else #inhibitory synapse
							forwardInputsI[dd] += weights[cc,dd];
						end
					end

				end #end if(spiked)
			end #end if(not refractory)

			#istdp
			if spiked[cc] && (t > stdpdelay)
				if cc <= Ne #excitatory neuron fired, potentiate i inputs
					for dd = (Ne+1):Ncells
						if weights_FF[dd,cc] == 0.
							continue
						end
						weights_FF[dd,cc] += eta*trace_istdp[dd]
						if weights_FF[dd,cc] > jeimax
							weights_FF[dd,cc] = jeimax
						end
					end
				else #inhibitory neuron fired, modify outputs to e neurons
					for dd = 1:Ne
						if weights_FF[cc,dd] == 0.
							continue
						end
						weights_FF[cc,dd] += eta*(trace_istdp[dd] - 2*r0*tauy)
						if weights_FF[cc,dd] > jeimax
							weights_FF[cc,dd] = jeimax
						elseif weights_FF[cc,dd] < jeimin
							weights_FF[cc,dd] = jeimin
						end
					end
				end
			end #end istdp


			#vstdp, ltd component
			if spiked[cc] && (t > stdpdelay) && (cc < Ne)
				for dd = 1:Ne #depress weights from cc to dd
					if weights_FF[cc,dd] == 0.
						continue
					end

					if u_vstdp[dd] > thetaltd
						weights_FF[cc,dd] -= dt*altd*(u_vstdp[dd]-thetaltd) #changed dt*
						if weights_FF[cc,dd] < jeemin
							weights_FF[cc,dd] = jeemin

						end
					end
				end
			end #end ltd

			#vstdp, ltp component
			if (t > stdpdelay) && (cc < Ne) && (v[cc] > thetaltp) && (v_vstdp[cc] > thetaltd)
				for dd = 1:Ne #potentiate weights from dd to cc
					if weights_FF[dd,cc] == 0.
						continue
					end

					weights_FF[dd,cc] += dt*altp*x_vstdp[dd]*(v[cc] - thetaltp)*(v_vstdp[cc] - thetaltd); #*((jeemax-weights[dd,cc])/(jeemax-jeemin))#weight dependent
					if weights_FF[dd,cc] > jeemax
						weights_FF[dd,cc] = jeemax
					end
				end
			end #end ltp

		end #end loop over cells
		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
		no_noisePrev = copy(no_noise)
		only_noisePrev = copy(only_noise)
	end #end loop over time
	@printf("\r")

	#times = times[:,1:maximum(ns)]


	return times,ns,Ne,Ncells, weights, weights_FF, T, relative_strengths, avg_conductances
end


#function that returns row vector with relative weight strength W_in/W_out for every cluster
#call this function in sim every 1000 iterations (every second) to update the matrix that contains this relative strenghts, sim has to return this matrix to be plotted in runsim.jl
function rel_weight_strengths(weights, popmembers, Ne)

	#populations
	Npop = size(popmembers,1) #number of assemblies
	Nmaxmembers = size(popmembers,2) #maximum number of neurons in a population

	#average strenghts per population
	relative_strengths = zeros(1,Npop)
	in_strengths = zeros(1,Npop)
	out_strengths = zeros(1,Npop)

	#average in_strength computation
	for ipop = 1:Npop
		in_strength = 0.
		N_in = 0.
		for ii = 1:Nmaxmembers-1
			for jj = ii+1:Nmaxmembers
				if (popmembers[ipop,jj] == -1) || (popmembers[ipop,jj] == 0)
					break
				end
				if weights[popmembers[ipop,ii],popmembers[ipop,jj]] == 0 && weights[popmembers[ipop,jj],popmembers[ipop,ii]] != 0
					in_strength += weights[popmembers[ipop,jj],popmembers[ipop,ii]]
					N_in += 1.
				elseif weights[popmembers[ipop,ii],popmembers[ipop,jj]] != 0 && weights[popmembers[ipop,jj],popmembers[ipop,ii]] == 0
					in_strength += weights[popmembers[ipop,ii],popmembers[ipop,jj]]
					N_in += 1.
				elseif weights[popmembers[ipop,ii],popmembers[ipop,jj]] == 0 && weights[popmembers[ipop,jj],popmembers[ipop,ii]] == 0
					in_strength += 0.
				else
					in_strength += weights[popmembers[ipop,ii],popmembers[ipop,jj]] + weights[popmembers[ipop,jj],popmembers[ipop,ii]]
					N_in += 2.
				end

			end
		end
		in_strengths[1,ipop] = in_strength/N_in
	end

	#average out_strenght computation
	for ipop =1:Npop
		out_strength = 0.
		N_out = 0.
		for ii = 1:Nmaxmembers
			if (popmembers[ipop,ii] == -1) || (popmembers[ipop,ii] == 0)
				break
			end
			for jj = 1:Ne
				if any(x->x==jj,popmembers[ipop,ii]) || (weights[popmembers[ipop,ii],jj] == 0 && weights[jj,popmembers[ipop,ii]] == 0)
					out_strength += 0.
				elseif weights[popmembers[ipop,ii],jj] == 0 && weights[jj,popmembers[ipop,ii]] != 0
					out_strength += weights[jj,popmembers[ipop,ii]]
					N_out += 1
				elseif weights[jj,popmembers[ipop,ii]] == 0 && weights[popmembers[ipop,ii],jj] != 0
					out_strength += weights[popmembers[ipop,ii],jj]
					N_out += 1
				else
					out_strength += weights[popmembers[ipop,ii],jj] + weights[jj,popmembers[ipop,ii]]
					N_out += 2
				end
			end
		end
		out_strengths[1,ipop] = out_strength/N_out
	end

	#average relative_strenght computation
	for i=1:Npop
		relative_strengths[1,i] = in_strengths[1,i]/out_strengths[1,i]
	end
	#println(string(relative_strengths))
	return relative_strengths
end


