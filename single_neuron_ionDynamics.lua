--------------------------------------------------------------
-- This script solves the cable equation with HH channels, 	--
-- activating synapses and transmission synapses.			--
--------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"cable_neuron"})

-- init UG
InitUG(3, AlgebraType("CPU", 1))


---------------------------------
-- read command line arguments --
---------------------------------
-- choice of grid
gridName = util.GetParam("-grid", "rat1.ugx")
gridSyn  = string.sub(gridName, 1, string.len(gridName)-4) .. "_syn.ugx"

-- parameters steering simulation
numRefs = util.GetParamNumber("-numRefs", 0)
dt = util.GetParamNumber("-dt", 1e-5)  -- in s
endTime = util.GetParamNumber("-endTime", dt)
nSteps = util.GetParamNumber("-nSteps", endTime/dt)
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")

-- synapse activity parameters
avg_start = util.GetParamNumber("-avgStart", 0.03)
avg_dur = util.GetParamNumber("-avgDur", 2.4e-4)
dev_start = util.GetParamNumber("-devStart", 0.015)
dev_dur = util.GetParamNumber("-devDur", 0.0)
num_synapses = util.GetParamNumber("-nSyn", 140)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")

-- file handling
filename = util.GetParam("-outName", "solution")
filename = filename.."/"


--------------------------
-- biological settings	--
--------------------------
-- settings are according to T. Branco

-- membrane conductances (in units of S/m^2)
g_k_ax = 400.0	-- axon
g_k_so = 200.0	-- soma
g_k_de = 30		-- dendrite

g_na_ax = 3.0e4
g_na_so = 1.5e3
g_na_de = 40.0

g_l_ax = 200.0
g_l_so = 1.0
g_l_de = 1.0

-- specific capacitance (in units of F/m^2)
spec_cap = 1.0e-2

-- resistivity (in units of Ohm m)
spec_res = 1.5

-- reversal potentials (in units of V)
e_k  = -0.09
e_na = 0.06
e_ca = 0.14

-- equilibrium concentrations (in units of mM)
-- comment: these concentrations will not yield Nernst potentials
-- as given above; pumps will have to be introduced to achieve this
-- in the case where Nernst potentials are calculated from concentrations!
k_out  = 4.0
na_out = 150.0
ca_out = 1.5

k_in   = 140.0
na_in  = 10.0
ca_in  = 5e-5

-- equilibrium potential (in units of V)
v_eq = -0.07

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.0e-9
diff_na	= 1.0e-9
diff_ca	= 2.2e-10

-- temperature in units of deg Celsius
temp = 37.0


------------------------------------
-- create domain and approx space --
------------------------------------
-- synapse distribution
synDistr = SynapseDistributor(gridName)
synDistr:clear() -- clear any synapses from grid
synDistr:place_synapses({0.0, 0.5, 0.5}, num_synapses, "AlphaPostSynapse")
if not synDistr:export_grid(gridSyn) then
	print("SynapseDistributor grid export unsuccessful. Aborting.")
	exit()
end

gridName = gridSyn
requiredSubsets = {"soma", "dend", "apic"}
dom = util.CreateDomain(gridName, numRefs, requiredSubsets)

-- check domain is acyclic
isAcyclic = is_acyclic(dom)
if not isAcyclic then
	print("Domain is not acyclic!")
	exit()
end


approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("ca", "Lagrange", 1)

approxSpace:init_levels();
approxSpace:init_surfaces();
approxSpace:init_top_surface();
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true);


--------------------
-- discretization --
--------------------
-- cable equation
CE = CableEquation("soma, dend, apic", true)

CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)

CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)
CE:set_rev_pot_ca(e_ca)

CE:set_k_out(k_out)
CE:set_na_out(na_out)
CE:set_ca_out(ca_out)

CE:set_diff_coeffs({diff_k, diff_na, diff_ca})

CE:set_temperature_celsius(temp)

ss_dend = "dend, apic"

-- Hodgkin and Huxley channels
HH = ChannelHH("v", "soma, dend, apic")
HH:set_conductances(g_k_so, g_na_so, "soma")
HH:set_conductances(g_k_de, g_na_de, ss_dend)

CE:add(HH)

-- leakage (exactly calibrated to achieve zero net current in equilibrium)
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leak = ChannelLeak("v", "soma, dend, apic")
leak:set_cond(g_l_so*tmp_fct, "soma")
leak:set_rev_pot(-0.0592363968739, "soma")
leak:set_cond(g_l_de*tmp_fct, ss_dend)
leak:set_rev_pot(-0.0679474583084, ss_dend)

CE:add(leak)


-- Calcium dynamics
vdcc = VDCC_BG_cable("ca", "soma, dend, apic")
ncx = NCX_cable("ca", "soma, dend, apic")
pmca = PMCA_cable("ca", "soma, dend, apic")
caLeak = IonLeakage("ca", "soma, dend, apic")
leakCaConst = -3.4836065573770491e-9 +	-- single pump PMCA flux density (mol/s/m^2)
			  -1.0135135135135137e-9 +	-- single pump NCX flux (mol/s/m^2)
			  3.3017662162505882e-11
caLeak:set_perm(leakCaConst, ca_in, ca_out, v_eq, 2)

CE:add(ncx)
CE:add(pmca)
CE:add(vdcc)
CE:add(caLeak)

--[[
-- Na-K pump
nak_ax = Na_K_Pump("", "axon")
nak_ax:set_max_flux(2.6481515257588432)	-- mol/(m^2*s)
nak_so = Na_K_Pump("", "soma")
nak_so:set_max_flux(6.05974e-7/4.57658e-06)	-- mol/(m^2*s)
nak_de = Na_K_Pump("", dendSubsets)
nak_de:set_max_flux(1.61593e-8/4.57658e-06)	-- mol/(m^2*s)

CE:add(nak_ax)
CE:add(nak_so)
CE:add(nak_de)


-- ion leakage
kLeak_ax = IonLeakage("k", "axon")
leakKConst_ax = 0.0000040675975261062531 +	-- HH (mol/s/m^2)
				-0.00000010983795579882983	-- Na/K (mol/s/m^2)
kLeak_ax:set_perm(leakKConst_ax, k_in, k_out, v_eq, 1)
kLeak_so = IonLeakage("k", "soma")
leakKConst_so = 2.0338e-06 +				-- HH (mol/s/m^2)
				-(2.0/3.0 * 6.05974e-7)		-- Na/K (mol/s/m^2)
kLeak_so:set_perm(leakKConst_so, k_in, k_out, v_eq, 1)
kLeak_de = IonLeakage("k", "soma")
leakKConst_de = 3.0507e-7 +					-- HH (mol/s/m^2)
				-(2.0/3.0 * 1.61593e-8)		-- Na/K (mol/s/m^2)
kLeak_de:set_perm(leakKConst_de, k_in, k_out, v_eq, 1)

-- TODO: What about Na!?

CE:add(kLeak_ax)
CE:add(kLeak_so)
CE:add(kLeak_de)
--]]


-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing_alpha(
	avg_start,   -- average onset of synaptical activity in [s]
	avg_dur/6.0, -- average tau of activity function in [s]
	dev_start,   -- deviation of onset in [s]
	dev_dur/6.0, -- deviation of tau in [s]
	1.2e-09)      -- peak conductivity in [S]
CE:set_synapse_handler(syn_handler)


--[[
-- electrode stimulation
-- 5nA seem to enervate the pyramidal cell with uniform diameters of 1um
-- (coords for 13-L3pyr-77.CNG.ugx, current given in C/ms)
CE:set_influx(5e-9, 6.54e-05, 2.665e-05, 3.985e-05, 0.0, 0.04)			-- near soma
CE:set_influx(5e-9, 3.955e-06, 1.095e-06, -3.365e-06, 0.001, 0.0025)		-- 1st edge soma to dend
CE:set_influx(0.3e-9, 3.955e-06, 1.095e-06, -3.365e-06, 0.0, 0.03)		-- 1st 1st edge soma to dend
CE:set_influx(0.095e-9, 0.0, 0.0, 0.0, 0.1, 0.1)							-- soma center vertex
CE:set_influx(0.2e-9, 0.0, 0.0, 0.0, 0.005, 0.0005)						-- soma center vertex
CE:set_influx(10.0e-9, 0.000139, 0.00020809, -2.037e-05, 0.005, 0.005)	-- distal apical dendrite vertex v1
CE:set_influx(10.0e-9, -3.96e-06, 0.0002173, -5.431e-05, 0.005, 0.005)	-- distal apical dendrite vertex v2
--]]


-- create domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)

assTuner = domainDisc:ass_tuner()


-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)


-- create operator from discretization
linOp = AssembledLinearOperator(timeDisc)


------------------
-- solver setup	--
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- linear solver --
linConvCheck = CompositeConvCheck(approxSpace, 1, 2e-26, 1e-08)
linConvCheck:set_component_check("v", 1e-21, 1e-12)
linConvCheck:set_verbose(verbose)

ilu = ILU()
linSolver = LinearSolver()
linSolver:set_preconditioner(ilu)
linSolver:set_convergence_check(linConvCheck)
--linSolver:set_debug(dbgWriter)


-------------------
-- time stepping --
-------------------
time = 0.0

-- init solution
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
Interpolate(v_eq, u, "v")
Interpolate(k_in, u, "k");
Interpolate(na_in, u, "na");
Interpolate(ca_in, u, "ca")


-- write start solution
if generateVTKoutput then 
	out = VTKOutput()
	out:print(filename.."vtk/solution", u, 0, time)
end


--[[
-- measurement setup
measFileVm = filename.."measVm.txt"
measFileCa = filename.."measCa.txt"

if ProcRank() == 0 then
	measOutVm = assert(io.open(measFileVm, "a"))
	measOutCa = assert(io.open(measFileCa, "a"))
end
--]]

-- store grid function in vector of  old solutions
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

curr_dt = dt
dtred = 2

lv = 0
maxLv = 10
cb_counter = {}
cb_counter[lv] = 0

while endTime-time > 0.001*curr_dt do
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, curr_dt)
	
	-- reduce time step if cfl < curr_dt
	-- (this needs to be done AFTER prepare_step as channels are updated there)
	dtChanged = false
	cfl = CE:estimate_cfl_cond(solTimeSeries:latest())
	print("estimated CFL condition: dt < " .. cfl)
	while (curr_dt > cfl) do
		curr_dt = curr_dt/dtred
		
		if lv+1 > maxLv then
			print("Time step too small.")
			exit()
		end
		
		lv = lv + 1
		cb_counter[lv] = 0
		print("estimated CFL condition: dt < " .. cfl .. " - reducing time step to " .. curr_dt)
		dtChanged = true
	end
	
	-- increase time step if cfl > curr_dt / dtred (and if time is aligned with new bigger step size)
	while curr_dt*dtred < cfl and lv > 0 and cb_counter[lv] % (dtred) == 0 do
		curr_dt = curr_dt*dtred;
		lv = lv - 1
		cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1]/dtred
		cb_counter[lv+1] = 0
		print ("estimated CFL condition: dt < " .. cfl .. " - increasing time step to " .. curr_dt)
		dtChanged = true
	end
	
	print("++++++ POINT IN TIME " .. math.floor((time+curr_dt)/curr_dt+0.5)*curr_dt .. " BEGIN ++++++")
	
	-- prepare again with new time step size
	if dtChanged == true then 
		timeDisc:prepare_step(solTimeSeries, curr_dt)
	end

	-- assemble linear problem
	matrixIsConst = time ~= 0.0 and dtChanged == false
	assTuner:set_matrix_is_const(matrixIsConst)
	AssembleLinearOperatorRhsAndSolution(linOp, u, b)
	
	-- synchronize (for profiling)
	PclDebugBarrierAll()
	
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(linOp, u, b, linSolver) == false then
		print("Could not apply linear solver.")
		exit()
	end
	
	--[[
	-- log time and vm in Soma
	if ProcRank() == 0 then
		vm_soma  = EvaluateAtClosestVertex(MakeVec(6.9e-07, 3.74e-06, -2.86e-06), 		u, "v", "soma", 		dom:subset_handler())
		vm_axon  = EvaluateAtClosestVertex(MakeVec(-4.05e-06, 6.736e-05, -1.341e-05), 	u, "v", "axon", 		dom:subset_handler())
		vm_dend  = EvaluateAtClosestVertex(MakeVec(-4.631e-05, -0.0001252, 4.62e-06), 	u, "v", "dendrite", 	dom:subset_handler())
		measOutVm:write(time, "\t", vm_soma, "\t", vm_axon, "\t", vm_dend, "\t", -65, "\n")
		ca_soma  = EvaluateAtClosestVertex(MakeVec(6.9e-07, 3.74e-06, -2.86e-06), 		u, "ca", "soma", 		dom:subset_handler())
		ca_axon  = EvaluateAtClosestVertex(MakeVec(-4.05e-06, 6.736e-05, -1.341e-05), 	u, "ca", "axon", 		dom:subset_handler())
		ca_dend  = EvaluateAtClosestVertex(MakeVec(-4.631e-05, -0.0001252, 4.62e-06), 	u, "ca", "dendrite", 	dom:subset_handler())
		measOutCa:write(time, "\t", ca_soma, "\t", ca_axon, "\t", ca_dend, "\t", -65, "\n")
	end
	--]]
	
	-- update to new time
	time = solTimeSeries:time(0) + curr_dt
	
	-- vtk output
	if (generateVTKoutput) then
		if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
			out:print(filename.."vtk/solution", u, math.floor(time/pstep+0.5), time)
		end
	end
	
	-- updte time series (reuse memory)
	oldestSol = solTimeSeries:oldest()
	VecScaleAssign(oldestSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	-- increment check-back counter
	cb_counter[lv] = cb_counter[lv] + 1

	print("++++++ POINT IN TIME " .. math.floor((time)/curr_dt+0.5)*curr_dt .. "  END ++++++")
end

-- end timeseries, produce gathering file
if generateVTKoutput then 
	out:write_time_pvd(filename.."vtk/solution", u) 
end

--[[
-- close measure file
if ProcRank() == 0 then
	measOutVm:close()
	measOutCa:close()
end
--]]

	
