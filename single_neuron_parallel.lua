--------------------------------------------------------------------------------
-- This script solves the cable equation with HH channels, activating         --
-- randomly (uniformly) distributed alpha synapses with randomly (normally)   --
-- distributed activation patterns on a L3 pyramidal cell.                    --
-- It is intended to be used in a parallel setup:                             --
-- Each core will then produce its own results (embarrassing parallelism);    --
-- a global histogram file for the number of evoked APs will be created       --
-- afterwards.                                                                --
--                                                                            --
-- Author:  Markus Breit                                                      --
-- Date:    2015-09-07                                                        --
--------------------------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")


-- dimension
dim = 3

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"cable_neuron"})


--------------------------------------------------------------------------------
-- Settings
--------------------------------------------------------------------------------
cell = util.GetParam("-cellName", "12-L3pyr")
if not cell == "12-L3pyr" or not cell == "31o_pyr" then
	exit("Cell not specified correctly. Type '12-L3pyr' or '31o_pyr'.")
end

rank = ProcRank()
rankAsString = string.format("%02d", rank)

if cell == "12-L3pyr" then
	gridName = "grids/13-L3pyr-77.CNG.ugx"
	gridSyn  = "grids/13-L3pyr-77.CNG_syn_p"..rankAsString..".ugx"
	gridDeg  = "grids/13-L3pyr-77.CNG_syn_p"..rankAsString.."_deg.ugx"
	distro   = {0.0, 0.0, 0.5, 0.5}
	neededSubsets = {"soma", "axon", "dendrite", "apical_dendrite"}
	dendSubsets = "dendrite, apical_dendrite"
else
	gridName = "grids/31o_pyramidal19aFI.CNG.ugx"
	gridSyn  = "grids/31o_pyramidal19aFI.CNG_syn_p"..rankAsString..".ugx"
	gridDeg  = "grids/31o_pyramidal19aFI.CNG_syn_p"..rankAsString.."_deg.ugx"
	distro   = {0.0, 1.0, 0.0}
	neededSubsets = {"soma", "dendrite", "axon"}
	dendSubsets = "dendrite"
end

-- parameters steering simulation
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			1e-5) -- in s
endTime		= util.GetParamNumber("-endTime",	 	0.01) -- in s
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)
pstep		= util.GetParamNumber("-pstep",			dt,		"plotting interval")

-- synapse activity parameters
avg_start = util.GetParamNumber("-avgStart"	,  0.03)
avg_dur = util.GetParamNumber(	"-avgDur"	,  2.4e-3)
dev_start = util.GetParamNumber("-devStart"	,  0.015)
dev_dur = util.GetParamNumber(	"-devDur"	,  0.0)
num_synapses = util.GetParamNumber("-nSyn", 140)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput	= util.HasParamOption("-vtk")

-- file handling
filename = util.GetParam("-outName", ".")
filename = filename.."/"

--------------------------------------------------------------------------------
-- Synapse distributions via plugin by Lukas Reinhardt
--------------------------------------------------------------------------------
---[[
synDistr = SynapseDistributor(gridName)
synDistr:clear() -- clear any synapses from grid
synDistr:place_synapses(distro, num_synapses, "AlphaPostSynapse")
export_succes = synDistr:export_grid(gridSyn)
print("SynapseDistributor grid export successful: " .. tostring(export_succes))
--]]

gridName = gridSyn

--------------------------------------------------------------------------------
-- Synapse degeneration
--------------------------------------------------------------------------------
---[[
deg_factor = util.GetParamNumber("-degFac", 0.5)
deg_factor = deg_factor + 0.5/num_synapses -- rounding instead of floor-ing

synDistr = SynapseDistributor(gridName)
if cell == "12-L3pyr" then
	synDistr:degenerate_uniform(deg_factor, 2) -- first factor means: newNumber = (1-factor)*oldNumber
	synDistr:degenerate_uniform(deg_factor, 3) -- second param is the subset index
else
	synDistr:degenerate_uniform(deg_factor, 1)
end
synDistr:print_status()
synDistr:export_grid(gridDeg)

gridName = gridDeg
--]]

--------------------------------------------------------------
-- File i/o setup for sample calcium concentration measurement
-------------------------------------------------------------- 
measFileVm = filename.."meas/measVm_p"..rankAsString..".txt"
measFileCa = filename.."meas/measCa_p"..rankAsString..".txt"

measOutVm = assert(io.open(measFileVm, "w"))
measOutCa = assert(io.open(measFileCa, "w"))

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
dom = Domain()
LoadDomain(dom, gridName, -2) -- means: load on every proc

-- check domain is acyclic
isAcyclic = is_acyclic(dom)
if not isAcyclic then
	print("Domain is not acyclic!")
	exit()
end

if numRefs > 0 then
	local refiner = GlobalDomainRefiner(dom)
	for i=1,numRefs do
		TerminateAbortedRun()
		refiner:refine()
	end
		
	delete(refiner)
end

--------------------------------------------------------------------------------
-- create Approximation Space
--------------------------------------------------------------------------------
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


-- cable equation
CE = CableEquation("soma, axon, " .. dendSubsets, true)

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
HH = ChannelHH("v", "axon, soma, " .. dendSubsets)
HH:set_conductances(g_k_ax, g_na_ax, "axon")
HH:set_conductances(g_k_so, g_na_so, "soma")
HH:set_conductances(g_k_de, g_na_de, dendSubsets)

CE:add(HH)


-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leak = ChannelLeak("v", "axon, soma, " .. dendSubsets)
leak:set_cond(g_l_ax*tmp_fct, "axon")
leak:set_rev_pot(-0.066148458, "axon")
leak:set_cond(g_l_so*tmp_fct, "soma")
leak:set_rev_pot(-0.030654022, "soma")
leak:set_cond(g_l_de*tmp_fct, dendSubsets)
leak:set_rev_pot(-0.057803624, dendSubsets)

CE:add(leak)


--Calcium dynamics
vdcc = VDCC_BG_cable("ca", "soma, " .. dendSubsets)
ncx = NCX_cable("ca", "soma, " .. dendSubsets)
pmca = PMCA_cable("ca", "soma, " .. dendSubsets)
caLeak = IonLeakage("ca", "soma, " .. dendSubsets)
leakCaConst = -3.4836065573770491e-9 +	-- single pump PMCA flux density (mol/s/m^2)
			  -1.0135135135135137e-9 +	-- single pump NCX flux (mol/s/m^2)
			  3.3017662162505882e-11
caLeak:set_perm(leakCaConst, ca_in, ca_out, v_eq, 2)

CE:add(ncx)
CE:add(pmca)
CE:add(vdcc)
CE:add(caLeak)


-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing_alpha(
	avg_start,	 -- average onset of synaptical activity in [s]
	avg_dur/6.0, -- average tau of activity function in [s]
	dev_start,   -- deviation of onset in [s]
	dev_dur/6.0, -- deviation of tau in [s]
	1.2e-9)		 -- peak conductivity in [S]
CE:set_synapse_handler(syn_handler)


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
-- linear solver --
linConvCheck = CompositeConvCheck(approxSpace, 20, 2e-26, 1e-08)
linConvCheck:set_component_check("v", 1e-21, 1e-12)
linConvCheck:set_verbose(verbose)

ilu = ILU()
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(linConvCheck)

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

min_dt = 1e-10
curr_dt = dt
dtred = 2

apCount = 0
apUp = false

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
	if cfl < min_dt then
		print("Required time step size is lower than admissible. Aborting.")
		exit()
	end
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
		if (generateVTKoutput) then 
			out:write_time_pvd(filename.."vtk/solution", u) 
		end
		exit()
	end
	
	-- log vm and calcium at soma
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
	
	-- count APs
	if vm_soma > 0 and not apUp then
		apCount = apCount + 1;
		apUp = true
	elseif vm_soma < -67 then
		apUp = false
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
measOutVm:close()
measOutCa:close()
	
-- delete geom files
--os.remove("../apps/cable/Ca_dyms/grids/13-L3pyr-77.CNG_syn_p"..rankAsString..".ugx")
os.remove("../apps/cable/Ca_dyms/grids/13-L3pyr-77.CNG_syn_p"..rankAsString.."_deg.ugx")


-- histogram of AP counts
apCountVec = {}
for i=1,50 do
	apCountVec[i] = 0
end
apCountVec[apCount+1] = 1

apCountVec = ParallelVecSum(apCountVec)

if rank == 0 then
	histOut = assert(io.open(filename.."meas/apHist.txt", "w"))
	for i=1,50 do
		histOut:write(i-1, "\t", apCountVec[i], "\n")
	end
	histOut:close()
end

