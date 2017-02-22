-------------------------------------------------------------------------------------
-- This script solves the cable equation on network morphologies created by NeuGen --
-- and converted to ugx by NETI.                                                   --
-- It (re-)configures the primary alpha synapses inserted into the grid by NeuGen  --
-- to a custom activation pattern. Simulation includes membrane potential,         --
-- K, Na, Ca as well. The ions dynamics include ion-specific channels              --
-- and pumps as well as leakage to ensure ion-wise equilibria.                     --
-- Biological parameters are due to T. Branco.                                     --
-------------------------------------------------------------------------------------

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

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- hierarchical distribution?
hDistr		= util.HasParamOption("-hDistr")

-- vertical interfaces?
bVertIntf	= util.HasParamOption("-vIntf")

-- vtk output?
generateVTKoutput	= util.HasParamOption("-vtk")

-- whether or not subsets are distinguished by layer
subsetsByLayer		= util.HasParamOption("-sslw")

-- profiling?
doProfiling			= util.HasParamOption("-profile")
SetOutputProfileStats(doProfiling)

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
CE = CableEquation(ss_axon..", "..ss_dend..", "..ss_soma, true)
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
HH = ChannelHHNernst("v", ss_axon..", "..ss_dend..", "..ss_soma)
HH:set_conductances(g_k_ax, g_na_ax, ss_axon)
HH:set_conductances(g_k_so, g_na_so, ss_soma)
HH:set_conductances(g_k_de, g_na_de, ss_dend)

CE:add(HH)

-- leakage
tmp_fct = math.pow(2.3, (temp-23.0)/10.0)

leak = ChannelLeak("v", ss_axon..", "..ss_dend..", "..ss_soma)
leak:set_cond(g_l_ax*tmp_fct, ss_axon)
leak:set_rev_pot(-0.066210342630746467, ss_axon)
leak:set_cond(g_l_so*tmp_fct, ss_soma)
leak:set_rev_pot(-0.022074360525636, ss_soma)
leak:set_cond(g_l_de*tmp_fct, ss_dend)
leak:set_rev_pot(-0.056314322586687, ss_dend)

CE:add(leak)


-- Na/K Pump
nak_ax = Na_K_Pump("", ss_axon)
nak_ax:set_max_flux(2.6481515257588432)	-- mol/(m^2*s)
nak_so = Na_K_Pump("", ss_soma)
nak_so:set_max_flux(6.05974e-7/4.57658e-06)	-- mol/(m^2*s)
nak_de = Na_K_Pump("", ss_dend)
nak_de:set_max_flux(1.61593e-8/4.57658e-06)	-- mol/(m^2*s)

CE:add(nak_ax)
CE:add(nak_so)
CE:add(nak_de)


-- ion leakage
kLeak_ax = IonLeakage("k", ss_axon)
leakKConst_ax = 0.0000040675975261062531 +  -- HH (mol/s/m^2)
			   -0.00000010983795579882983   -- Na/K (mol/s/m^2)
kLeak_ax:set_perm(leakKConst_ax, k_in, k_out, v_eq, 1)
kLeak_so = IonLeakage("k", ss_soma)
leakKConst_so = 2.0338e-06 +			-- HH (mol/s/m^2)
				-(2.0/3.0 * 6.05974e-7)	-- Na/K (mol/s/m^2)
kLeak_so:set_perm(leakKConst_so, k_in, k_out, v_eq, 1)
kLeak_de = IonLeakage("k", ss_dend)
leakKConst_de = 3.0507e-7 +	            -- HH (mol/s/m^2)
				-(2.0/3.0 * 1.61593e-8)	-- Na/K (mol/s/m^2)
kLeak_de:set_perm(leakKConst_de, k_in, k_out, v_eq, 1)

-- TODO: WHAt about Na!?

CE:add(kLeak_ax)
CE:add(kLeak_so)
CE:add(kLeak_de)


-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing_alpha(
	5e-3,	-- average onset of synaptical activity in [s]
	4e-4,   -- average tau of activity function in [s]
	2.5e-3, -- deviation of onset in [s]
	1.5e-5, -- deviation of tau in [s]
	1.2e-9)	-- peak conductivity in [S]
CE:set_synapse_handler(syn_handler)

-- domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)

assTuner = domainDisc:ass_tuner()

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
linConvCheck:set_adaptive(true)

ilu = ILU()
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(linConvCheck)
--cgSolver:set_debug(dbgWriter)


-------------------------
-- domain distribution --
-------------------------
-- Domain distribution needs to be performed AFTER addition
-- of the synapse handler to the CE object and addition of the
-- CE object to the domain disc (i.e.: when the synapse handler
-- has got access to the grid).
-- The reason is that the synapse handler needs access to the grid
-- to correctly distribute the synapse* attachments.

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

--print("Saving parallel grid layout")
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 1e-5)

-- speed up assembling (if vertical interfaces are present)
if bVertIntf then
	cableAssTuner = CableAssTuner(domainDisc, approxSpace)
	cableAssTuner:remove_ghosts_from_assembling_iterator()
end

order_cuthillmckee(approxSpace);

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
Interpolate(k_in, u, "k");
Interpolate(na_in, u, "na");
Interpolate(ca_in, u, "ca")

-- write start solution
if (generateVTKoutput) then 
	out = VTKOutput()
	--out:print(fileName .."vtk/Solvung", u, 0, time)
	out:print_subset(fileName .."vtk/somatic_signals", u, 2, 0, time)
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
	AssembleLinearOperatorRhsAndSolution(linOp, u, b)
	
	-- synchronize (for profiling)
	PclDebugBarrierAll()
	
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(linOp, u, b, cgSolver) == false then
		print("Could not apply linear solver.");
		exit()
	end
	
	-- update to new time
	time = solTimeSeries:time(0) + curr_dt
	
	-- vtk output
	if (generateVTKoutput) then
		if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
			--out:print(fileName .."vtk/Solvung", u, math.floor(time/pstep+0.5), time)
			out:print_subset(fileName .."vtk/somatic_signals", u, 2, math.floor(time/pstep+0.5), time)
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
