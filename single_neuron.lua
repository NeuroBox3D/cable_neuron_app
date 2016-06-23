--------------------------------------------------------------
-- This script solves the cable equation with HH channels, 	--
-- activating synapses and transmission synapses.			--
--------------------------------------------------------------

print("scrypt start")
-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")


--------------------------------------------------------------------------------
-- Cell specification
--------------------------------------------------------------------------------
cell = util.GetParam("-cellName", "12-L3pyr")
if not cell == "12-L3pyr" or not cell == "31o_pyr" then
	exit("Cell not specified correctly. Type '12-L3pyr' or '31o_pyr'.")
end

print("cell specs works")

--------------------------------------------------------------------------------
-- UG4-Standard-Settings
--------------------------------------------------------------------------------

-- choice of grid
--gridName = util.GetParam("-grid", "grids/test_cell_small_ref_1.ugx")
--gridName = util.GetParam("-grid", "grids/31o_pyramidal19aFI.CNG_with_subsets.ugx")
--gridName = util.GetParam("-grid", "grids/31o_pyramidal19aFI.CNG_with_subsets_and_diams.ugx")
--gridName = util.GetParam("-grid", "grids/31o_pyramidal19aFI.CNG_diams.ugx")
--gridName = util.GetParam("-grid", "grids/13-L3pyr-77.CNG.ugx")

if cell == "12-L3pyr" then
	gridName = util.GetParam("-grid", "../apps/cable/Ca_dyms/grids/13-L3pyr-77.CNG_syn_deg.ugx")
else
	gridName = util.GetParam("-grid", "31o_pyramidal19aFI.CNG_diams_syn.ugx")
end

print(gridName);

-- dimension
dim = 3

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"cable_neuron", "NeuronalTopologyImporter"})

-- parameters steering simulation
numPreRefs	= util.GetParamNumber("-numPreRefs",	0)
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			1e-5) -- in s
endTime		= util.GetParamNumber("-endTime",	  	10.0) -- in s
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)
pstep		= util.GetParamNumber("-pstep",			dt,		"plotting interval")


print(" chosen parameters:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)
print("    pstep       = " .. pstep)

-- Synapse activity parameters
avg_start = util.GetParamNumber("-avgStart"	,  30.0)
avg_dur = util.GetParamNumber(	"-avgDur"	,   2.4)
dev_start = util.GetParamNumber("-devStart"	,  15.0)
dev_dur = util.GetParamNumber(	"-devDur"	,   0.0)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput	= util.HasParamOption("-vtk")

-- file handling
filename = util.GetParam("-outName", "sol_new_clearance_1e-3")
filename = filename.."/"

--------------------------------------------------------------
-- File i/o setup for sample calcium concentration measurement
-------------------------------------------------------------- 
measFileVm = filename.."meas/measVm.txt"
measFileCa = filename.."meas/measCa.txt"

if ProcRank() == 0 then
	measOutVm = assert(io.open(measFileVm, "a"))
	measOutCa = assert(io.open(measFileCa, "a"))
end

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
v_eq = -0.065

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.0e-9
diff_na	= 1.0e-9
diff_ca	= 2.2e-10

-- temperature in units of deg Celsius
temp = 37.0


--------------------------------------------------------------------------------
-- Create, Load, Refine Domain
--------------------------------------------------------------------------------
if cell == "12-L3pyr" then
	neededSubsets = {"soma", "axon", "dendrite", "apical_dendrite"}
else
	neededSubsets = {"soma", "dendrite", "axon"}
end
--dom = util.CreateDomain(gridName, numRefs, neededSubsets)
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, "metis")


--------------------------------------------------------------------------------
-- Synapse distributions via plugin by Lukas Reinhardt
-- (SplitSynapses)
--------------------------------------------------------------------------------
--30 alphasynapses (15 post- and 15 presynapses)
alphasyns = AlphaSynapses(0,15)
alphasyns:set_mean_gMax(3)
alphasyns:set_dev_gMax(1)
alphasyns:set_mean_onset(2)
alphasyns:set_dev_onset(0.5)
alphasyns:set_mean_tau(1)
alphasyns:set_dev_tau(0.3)
alphasyns:set_mean_e(0.5)
alphasyns:set_dev_e(0.1)

--60 exp2synapses
exp2syns = Exp2Synapses(14,30)
exp2syns:set_mean_tau1(3)
exp2syns:set_dev_tau1(1)
exp2syns:set_mean_tau2(4)
exp2syns:set_dev_tau2(1)
exp2syns:set_mean_e(5)
exp2syns:set_dev_e(0.5)
exp2syns:set_mean_w(1)
exp2syns:set_dev_w(0.2)

-- Instantiate SplitSynapseDistributor object and distribute synapses on the grid
synDistr = SplitSynapseDistributor(gridName, gridName.."_out.ugx", false)
synDistr:place_synapses_uniform(alphasyns:get_synapses())
synDistr:place_synapses_uniform(exp2syns:get_synapses())
print(synDistr:export_grid())

--------------------------------------------------------------------------------
-- Distribute Domain
--------------------------------------------------------------------------------
--util.DistributeDomain(dom, "metis")

--print("Saving parallel grid layout")
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 1e-5)


--------------------------------------------------------------------------------
-- create Approximation Space
--------------------------------------------------------------------------------
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("ca", "Lagrange", 1)

approxSpace:init_levels();
approxSpace:init_surfaces();
approxSpace:init_top_surface();
approxSpace:print_layout_statistic()
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true);


--------------------------------------------------------------------------------
--CableEquation constructor creates every needed concentration out of added Channels from Channel list
--------------------------------------------------------------------------------
-- cable equation
if cell == "12-L3pyr" then
	CE = CableEquation("soma, axon, dendrite, apical_dendrite")
else
	CE = CableEquation("soma, dendrite, axon")
end

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


-- Hodgkin and Huxley channels
--HH = ChannelHHNernst("v, k, na", "axon")
HHaxon = ChannelHH("v", "axon")
HHaxon:set_conductances(g_k_ax, g_na_ax)
HHsoma = ChannelHH("v", "soma")
HHsoma:set_conductances(g_k_so, g_na_so)
if cell == "12-L3pyr" then
	HHdend = ChannelHH("v", "dendrite, apical_dendrite")
else
	HHdend = ChannelHH("v", "dendrite")
end
HHdend:set_conductances(g_k_de, g_na_de)

CE:add(HHaxon)
CE:add(HHsoma)
CE:add(HHdend)


--Calcium dynamics
vdcc = VDCC_BG_cable("ca", "dendrite, soma, apical_dendrite")
ncx = NCX_cable("v, ca", "dendrite, soma, apical_dendrite")
pmca = PMCA_cable("v, ca", "dendrite, soma, apical_dendrite")
caLeak = IonLeakage("", "dendrite, soma, apical_dendrite")
caLeak:set_leaking_quantity("ca")
leakCaConst = -3.4836065573770491e-9 +	-- single pump PMCA flux density (mol/s/m^2)
			  -1.0135135135135137e-9 +	-- single pump NCX flux (mol/s/m^2)
			  3.3017662162505882e-11
caLeak:set_perm(leakCaConst, ca_in, ca_out, v_eq)

CE:add(ncx)
CE:add(pmca)
CE:add(vdcc)
CE:add(caLeak)


-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leakAxon = ChannelLeak("v", "axon")
leakAxon:set_cond(g_l_ax*tmp_fct)
leakAxon:set_rev_pot(-0.066148458)
leakSoma = ChannelLeak("v", "soma")
leakSoma:set_cond(g_l_so*tmp_fct)
leakSoma:set_rev_pot(-0.030654022)
if cell == "12-L3pyr" then
	leakDend = ChannelLeak("v", "dendrite, apical_dendrite")
else
	leakDend = ChannelLeak("v", "dendrite")
end
leakDend:set_cond(g_l_de*tmp_fct)
leakDend:set_rev_pot(-0.057803624)

CE:add(leakAxon)
CE:add(leakSoma)
CE:add(leakDend)


-- synapses
syn_handler = SplitSynapseHandler()
--syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing(
	avg_start,	-- average start time of synaptical activity in ms
	avg_dur,	-- average duration of activity in ms (10)
	dev_start,	-- deviation of start time in ms
	dev_dur,	-- deviation of duration in ms
	1.2e-3)		-- peak conductivity in [uS]
CE:set_synapse_handler(syn_handler)

--CE:set_synapse_distributor(sd)



--------------------------------------------------------------------------------
--	ELECTRODE STIMULATION SETUP
--------------------------------------------------------------------------------

--  INFO: coords for 31o_pyramidal19aFI.CNG_with_subsets.ugx --

--	ELECTRODE STIMULATION near soma: 5nA seem to inervate the pyramidal cell with uniform diameters of 1um
--CE:set_influx(5e-9, 6.54e-05, 2.665e-05, 3.985e-05, 0.0, 0.04) -- current given in A

--	!!! DO NOT USE!!! ELECTRODE STIMULATION into soma: 5nA seem to enervate the pyramidal cell with attached diameters (also causing backpropagating APs)
--CE:set_influx(5e-9, 6.9e-06, 3.74e-05, -2.86e-05, 0.0, 0.04) -- current given in A


-- INFO: coords for 13-L3pyr-77.CNG.ugx --current given in C/ms --

--CE:set_influx(5e-9, 3.955e-06, 1.095e-06, -3.365e-06, 0.001, 0.0025) -- 1st edge soma to dend
--CE:set_influx(0.3e-9, 3.955e-06, 1.095e-06, -3.365e-06, 0.0, 0.03) -- 1st 1st edge soma to dend
--CE:set_influx(0.095e-9, 0.0, 0.0, 0.0, 0.1, 0.1) -- soma center vertex
--CE:set_influx(0.2e-9, 0.0, 0.0, 0.0, 0.005, 0.0005) -- soma center vertex
--CE:set_influx(10.0e-9, 0.000139, 0.00020809, -2.037e-05, 0.005, 0.005) -- distal apical dendrite vertex v1
--CE:set_influx(10.0e-9, -3.96e-06, 0.0002173, -5.431e-05, 0.005, 0.005) -- distal apical dendrite vertex v2
-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)

assTuner = domainDisc:ass_tuner()

-------------------------------------------
--  Setup Time Discretization
-------------------------------------------

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)


-------------------------------------------
--  Algebra
-------------------------------------------
-- create operator from discretization
linOp = AssembledLinearOperator(timeDisc)

------------------
-- solver setup	--
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- linear solver --
linConvCheck = CompositeConvCheck(approxSpace, 20, 2e-26, 1e-08)
linConvCheck:set_component_check("v", 1e-21, 1e-12)
linConvCheck:set_verbose(verbose)

ilu = ILU()
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(linConvCheck)
--cgSolver:set_debug(dbgWriter)

----------------------
-- time stepping	--
----------------------

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
if (generateVTKoutput) then 
	out = VTKOutput()
	out:print(filename.."vtk/solution", u, 0, time)
end

-- store grid function in vector of  old solutions
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

curr_dt = dt
dtred = 2

lv = 0
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
	if AssembleLinearOperatorRhsAndSolution(linOp, u, b) == false then 
		print("Could not assemble operator"); exit(); 
	end
	
	-- synchronize (for profiling)
	PclDebugBarrierAll()
	
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(linOp, u, b, cgSolver) == false then
		print("Could not apply linear solver.");
	end
	
	-- log time and vm in Soma
	if ProcRank() == 0 then
		if cell == "12-L3pyr" then
			vm_soma  = EvaluateAtClosestVertex(MakeVec(0.0, 0.0, 0.0), 						u, "v", "soma", 		dom:subset_handler())
			vm_axon  = EvaluateAtClosestVertex(MakeVec(-3.828e-05, -0.00013166, -2.34e-05), u, "v", "axon", 		dom:subset_handler())
			vm_dend  = EvaluateAtClosestVertex(MakeVec(8.304e-05, -1.982e-05, -8.4e-06), 	u, "v", "dendrite", 		dom:subset_handler())
			vm_aDend = EvaluateAtClosestVertex(MakeVec(-3.84e-06, 0.00018561, -3.947e-05), 	u, "v", "apical_dendrite", 	dom:subset_handler())
			measOutVm:write(time, "\t", vm_soma, "\t", vm_axon, "\t", vm_dend, "\t", vm_aDend, "\n")
			ca_soma  = EvaluateAtClosestVertex(MakeVec(0.0, 0.0, 0.0), 						u, "ca", "soma", 		dom:subset_handler())
			ca_axon  = EvaluateAtClosestVertex(MakeVec(-3.828e-05, -0.00013166, -2.34e-05), u, "ca", "axon", 		dom:subset_handler())
			ca_dend  = EvaluateAtClosestVertex(MakeVec(8.304e-05, -1.982e-05, -8.4e-06), 	u, "ca", "dendrite", 		dom:subset_handler())
			ca_aDend = EvaluateAtClosestVertex(MakeVec(-3.84e-06, 0.00018561, -3.947e-05), 	u, "ca", "apical_dendrite", 	dom:subset_handler())
			measOutCa:write(time, "\t", ca_soma, "\t", ca_axon, "\t", ca_dend, "\t", ca_aDend, "\n")
		else	
			vm_soma  = EvaluateAtClosestVertex(MakeVec(6.9e-07, 3.74e-06, -2.86e-06), 		u, "v", "soma", 		dom:subset_handler())
			vm_axon  = EvaluateAtClosestVertex(MakeVec(-4.05e-06, 6.736e-05, -1.341e-05), 	u, "v", "axon", 		dom:subset_handler())
			vm_dend  = EvaluateAtClosestVertex(MakeVec(-4.631e-05, -0.0001252, 4.62e-06), 	u, "v", "dendrite", 	dom:subset_handler())
			measOutVm:write(time, "\t", vm_soma, "\t", vm_axon, "\t", vm_dend, "\t", -65, "\n")
			ca_soma  = EvaluateAtClosestVertex(MakeVec(6.9e-07, 3.74e-06, -2.86e-06), 		u, "ca", "soma", 		dom:subset_handler())
			ca_axon  = EvaluateAtClosestVertex(MakeVec(-4.05e-06, 6.736e-05, -1.341e-05), 	u, "ca", "axon", 		dom:subset_handler())
			ca_dend  = EvaluateAtClosestVertex(MakeVec(-4.631e-05, -0.0001252, 4.62e-06), 	u, "ca", "dendrite", 	dom:subset_handler())
			measOutCa:write(time, "\t", ca_soma, "\t", ca_axon, "\t", ca_dend, "\t", -65, "\n")
		end
	end
	
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
if (generateVTKoutput) then 
	out:write_time_pvd(filename.."vtk/solution", u) 
end

-- close measure file
if ProcRank() == 0 then
	measOutVm:close()
	measOutCa:close()
end
	
