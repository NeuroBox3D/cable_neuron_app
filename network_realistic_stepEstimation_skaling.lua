------------------------------------------------------
-- This script is intended for testing purposes.	--
-- It solves the cable equation with HH channels,	--
-- activating synapses and transmission synapses.	--
------------------------------------------------------

ug_load_script("ug_util.lua")
--ug_load_script("util/load_balancing_util.lua")

-- choice of grid only local needed
gridName = util.GetParam("-grid", "testNetwork.ugx")

-- dimension
dim = 3
print("Detected dimension "..dim.." in ugx file.\n")

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"SynapseHandler","HH_Kabelnew"})

-- parameters steering simulation
numPreRefs	= util.GetParamNumber("-numPreRefs",	0)
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			1e-5) -- in units of s
endTime		= util.GetParamNumber("-endTime",		0.02) -- in units of s
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)
pstep		= util.GetParamNumber("-pstep",			dt,		"plotting interval")
imbFactor	= util.GetParamNumber("-imb",			1.05,	"imbalance factor")

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- hierarchical distribution?
hDistr		= util.HasParamOption("-hDistr")

-- vertical interfaces?
bVertIntf	= util.HasParamOption("-vIntf")

-- vtk output?
generateVTKoutput	= util.HasParamOption("-vtk")

-- profiling?
doProfiling			= util.HasParamOption("-profile")
SetOutputProfileStats(doProfiling)

-- file handling
fileName = util.GetParam("-outName", "Solvung")


--------------------------
-- biological settings	--
--------------------------

-- settings are according to T. Branco

-- membrane conductances (in units of S/m^2)
g_k_ax = 400.0	-- axon
g_k_so = 200.0	-- soma
g_k_de = 30.0	-- dendrite

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

----------------------------------
-- setup approximation space	--
----------------------------------
-- create, load, refine and distribute domain

write(">> Loading, distributing and refining domain ...")

neededSubsets = {"Axon", "Dendrite", "Soma"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, "metis")

--print("Saving parallel grid layout")
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 1e-5)

-- create Approximation Space
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("ca", "Lagrange", 1)

-- gating functions for HH-Fluxes
approxSpace:init_levels();
approxSpace:init_surfaces();
approxSpace:init_top_surface();
approxSpace:print_layout_statistic()
approxSpace:print_statistic()
order_cuthillmckee(approxSpace);

----------------------
-- setup elem discs	--
----------------------

-- cable equation
CE = CableEquation("Axon, Dendrite, Soma, PreSynapseEdges, PostSynapseEdges")
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
--HH = ChannelHHNernst("v, k, na", "Axon")
HHaxon = ChannelHH("v", "Axon, PreSynapseEdges")
HHaxon:set_conductances(g_k_ax, g_na_ax)
HHsoma = ChannelHH("v", "Soma")
HHsoma:set_conductances(g_k_so, g_na_so)
HHdend = ChannelHH("v", "Dendrite, PostSynapseEdges")
HHdend:set_conductances(g_k_de, g_na_de)

CE:add(HHaxon)
CE:add(HHsoma)
CE:add(HHdend)

-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leakAxon = ChannelLeak("v", "Axon, PreSynapseEdges")
leakAxon:set_cond(g_l_ax*tmp_fct)
leakAxon:set_rev_pot(-0.066148458)
leakSoma = ChannelLeak("v", "Soma")
leakSoma:set_cond(g_l_so*tmp_fct)
leakSoma:set_rev_pot(-0.030654022)
leakDend = ChannelLeak("v", "Dendrite, PostSynapseEdges")
leakDend:set_cond(g_l_de*tmp_fct)
leakDend:set_rev_pot(-0.057803624)

CE:add(leakAxon)
CE:add(leakSoma)
CE:add(leakDend)

-- synapses
syn_handler = NETISynapseHandler()
syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing(
	0.1,	-- average start time of synaptical activity in ms
	5,		-- average duration of activity in ms (10)
	1.0,	-- deviation of start time in ms
	0.5,	-- deviation of duration in ms
	1.2e-3)	-- peak conductivity in units of uS (6e-4)
CE:set_synapse_handler(syn_handler)

-- domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)

assTuner = domainDisc:ass_tuner()

-- speed up assembling
cableAssTuner = CableAssTuner(domainDisc, approxSpace)
cableAssTuner:remove_ghosts_from_assembling_iterator()

-- time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)

-- create instationary operator
linOp = AssembledLinearOperator(timeDisc)

------------------
-- solver setup	--
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- linear solver --
linConvCheck = CompositeConvCheck3dCPU1(approxSpace, 20, 2e-26, 1e-08)
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
	--out:print(fileName .."vtk/Solvung", u, 0, time)
	out:print_subset(fileName .."somatic_signals", u, 2, 0, time)
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
	
	-- update to new time
	time = solTimeSeries:time(0) + curr_dt
	
	-- vtk output
	if (generateVTKoutput) then
		if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
			--out:print(fileName .."vtk/Solvung", u, math.floor(time/pstep+0.5), time)
			out:print_subset(fileName .."somatic_signals", u, 2, math.floor(time/pstep+0.5), time)
		end
	end
	
	-- updte time series (reuse memory)
	oldestSol = solTimeSeries:oldest()
	VecScaleAssign(oldestSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	-- increment check-back counter
	cb_counter[lv] = cb_counter[lv] + 1

	print("++++++ POINT IN TIME " .. math.floor((time)/curr_dt+0.5)*curr_dt .. "  END ++++++")
	
	-- synchronize (for profiling)
	PclDebugBarrierAll()
end

-- end timeseries, produce gathering file
if (generateVTKoutput) then
	--out:write_time_pvd(fileName .."vtk/Solvung", u)
	out:write_time_pvd(fileName .."somatic_signals", u)
end

if doProfiling then
	WriteProfileData(fileName .."pd.pdxml")
end

