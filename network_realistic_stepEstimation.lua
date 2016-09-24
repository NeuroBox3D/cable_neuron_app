------------------------------------------------------
-- This script is intended for testing purposes.	--
-- It solves the cable equation with HH channels,	--
-- activating synapses and transmission synapses.	--
------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- choice of grid
gridName = util.GetParam("-grid", "testNetwork.ugx")

-- dimension
dim = 3

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"cable_neuron"})

-- parameters steering simulation
numPreRefs	= util.GetParamNumber("-numPreRefs",	0)
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			1e-5) -- in units of s
endTime		= util.GetParamNumber("-endTime",		0.1) -- in units of s
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)
pstep		= util.GetParamNumber("-pstep",			dt,		"plotting interval")
imbFactor	= util.GetParamNumber("-imb",			1.05,	"imbalance factor")

-- with simulation of single ion concentrations?
withIons = util.HasParamOption("-ions")

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

-- whether or not subsets are distinguished by layer
subsetsByLayer		= util.HasParamOption("-sslw")

-- file handling
fileName = util.GetParam("-outName", "Solvung")
fileName = fileName.."/"


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

if not subsetsByLayer then
	neededSubsets = {"Axon", "Dendrite", "Soma"}
else
	neededSubsets = {}
end
dom = util.CreateDomain(gridName, numPreRefs, neededSubsets)

-- check domain is acyclic
isAcyclic = is_acyclic(dom)
if not isAcyclic then
	print("Domain is not acyclic!")
	exit()
end

balancer.partitioner = "parmetis"

-- balancer.firstDistLvl = -1 will cause immediate distribution to all procs if redistSteps == 0,
-- but will cause first distribution to occur on level redistSteps otherwise!
-- 0 will distribute to firstDistProcs on (grid-)level 0 and then each proc
-- will redistribute to redistProcs on levels i*redistStep (i=1,2,...)
-- AND ALL THIS ONLY if staticProcHierarchy is set to true!

if hDistr == true then
	balancer.firstDistLvl 		= 0
	balancer.redistSteps 		= 1
	balancer.firstDistProcs		= 256
	balancer.redistProcs		= 256
else
	balancer.firstDistLvl		= -1
	balancer.redistSteps		= 0
end

balancer.imbalanceFactor		= imbFactor
balancer.staticProcHierarchy	= true
balancer.ParseParameters()
balancer.PrintParameters()

-- in parallel environments: use a load balancer to distribute the grid
loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	loadBalancer:enable_vertical_interface_creation(bVertIntf)
	balancer.RefineAndRebalanceDomain(dom, numRefs, loadBalancer)

	print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end

print(dom:domain_info():to_string())


write(">> done\n")

--print("Saving parallel grid layout")
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 1e-5)

-- create Approximation Space
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
if withIons == true then
	approxSpace:add_fct("k", "Lagrange", 1)
	approxSpace:add_fct("na", "Lagrange", 1)
	approxSpace:add_fct("ca", "Lagrange", 1)
end

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

-- collect subsets defined layer-wise
ss_axon = ""
ss_dend = ""
ss_soma = ""
if subsetsByLayer then
	ss_axon = "AXON__L4_STELLATE, AXON__L23_PYRAMIDAL, AXON__L5A_PYRAMIDAL, AXON__L5B_PYRAMIDAL"
	ss_dend = "DEND__L4_STELLATE, DEND__L23_PYRAMIDAL, DEND__L5A_PYRAMIDAL, DEND__L5B_PYRAMIDAL"
	ss_soma = "SOMA__L4_STELLATE, SOMA__L23_PYRAMIDAL, SOMA__L5A_PYRAMIDAL, SOMA__L5B_PYRAMIDAL"
else
	ss_axon = "Axon"
	ss_dend = "Dendrite"
	ss_soma = "Soma"
end

-- cable equation
CE = CableEquation(ss_axon..", "..ss_dend..", "..ss_soma..", PreSynapseEdges, PostSynapseEdges", withIons)
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
if withIons == true then
	HH = ChannelHHNernst("v", ss_axon..", "..ss_dend..", "..ss_soma..", PreSynapseEdges, PostSynapseEdges")
else
	HH = ChannelHH("v", ss_axon..", "..ss_dend..", "..ss_soma..", PreSynapseEdges, PostSynapseEdges")
end
HH:set_conductances(g_k_ax, g_na_ax, ss_axon..", PreSynapseEdges")
HH:set_conductances(g_k_so, g_na_so, ss_soma)
HH:set_conductances(g_k_de, g_na_de, ss_dend..", PostSynapseEdges")

CE:add(HH)

-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leak = ChannelLeak("v", ss_axon..", "..ss_dend..", "..ss_soma..", PreSynapseEdges, PostSynapseEdges")
leak:set_cond(g_l_ax*tmp_fct, ss_axon..", PreSynapseEdges")
leak:set_rev_pot(-0.066148458, ss_axon..", PreSynapseEdges")
leak:set_cond(g_l_so*tmp_fct, ss_soma)
leak:set_rev_pot(-0.030654022, ss_soma)
leak:set_cond(g_l_de*tmp_fct, ss_dend..", PostSynapseEdges")
leak:set_rev_pot(-0.057803624, ss_dend..", PostSynapseEdges")

CE:add(leak)

-- synapses
--[[
syn_handler = NETISynapseHandler()
syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing(
       0.0,    -- average start time of synaptical activity in ms
       2.4,    -- average duration of activity in ms
       0.0,    -- deviation of start time in ms
       0.1,    -- deviation of duration in ms
       1.2e-3) -- peak conductivity in units of uS
CE:set_synapse_handler(syn_handler)
--]]

syn_handler = SplitSynapseHandler()
--syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_ce_object(CE)

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


----------------------------
-- set synapse parameters --
----------------------------
syn_handler:show_status()
it = syn_handler:begin_AlphaPreSynapse()
it_end = syn_handler:end_AlphaPreSynapse()
while it:inequal(it_end) do
	it:get():set_onset(1337)
	it:next()
end

it = syn_handler:begin_AlphaPostSynapse()
it_end = syn_handler:end_AlphaPostSynapse()
while it:inequal(it_end) do
	it:get():set_tau(42)
	it:next()
end

it = syn_handler:begin_Exp2PreSynapse()
it_end = syn_handler:end_Exp2PreSynapse()
while it:inequal(it_end) do
	it:get():set_threshold(-0.42)
	it:next()
end

it = syn_handler:begin_Exp2PostSynapse()
it_end = syn_handler:end_Exp2PostSynapse()
while it:inequal(it_end) do
	it:get():set_tau1(666)
	it:get():set_tau2(420)
	it:next()
end

syn_handler:show_status()



----------------------
-- time stepping	--
----------------------
if generateActivityStats then
	syn_handler:print_synapse_statistics(2)
end

time = 0.0

-- init solution
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
Interpolate(v_eq, u, "v")
if withIons == true then
	Interpolate(k_in, u, "k");
	Interpolate(na_in, u, "na");
	Interpolate(ca_in, u, "ca")
end

-- write start solution
if (generateVTKoutput) then 
	out = VTKOutput()
	out:print(fileName .."vtk/Solvung", u, 0, time)
	--out:print_subsets(fileName .."vtk/somatic_signals", u, ss_soma, 0, time)
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
			out:print(fileName .."vtk/Solvung", u, math.floor(time/pstep+0.5), time)
			--out:print_subsets(fileName .."vtk/somatic_signals", u, ss_soma, math.floor(time/pstep+0.5), time)
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
	out:write_time_pvd(fileName .."vtk/somatic_signals", u)
end

if doProfiling then
	WriteProfileData(fileName .."pd.pdxml")
end

