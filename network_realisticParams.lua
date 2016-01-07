------------------------------------------------------
-- This script is intended for testing purposes.	--
-- It solves the cable equation with HH channels,	--
-- activating synapses and transmission synapses.	--
------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- choice of grid
gridName = util.GetParam("-grid", "testNetwork.ugx")

-- dimension
ugxfi = UGXFileInfo()
ugxfi:parse_file(gridName)
dim = ugxfi:physical_grid_dimension(0)
print("Detected dimension "..dim.." in ugx file.\n")

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"SynapseHandler","HH_Kabelnew"})

-- parameters steering simulation
numPreRefs	= util.GetParamNumber("-numPreRefs",	0)
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			0.01) -- in ms
endTime		= util.GetParamNumber("-endTime",		100.0) -- in ms
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)
pstep		= util.GetParamNumber("-pstep",			dt,		"plotting interval")


-- vtk output?
generateVTKoutput	= util.HasParamOption("-vtk")

-- file handling
fileName = util.GetParam("-outName", "Solvung")
fileName = fileName.."/"


--------------------------
-- biological settings	--
--------------------------

-- settings are according to T. Branco

-- membrane conductances (in units of C/m^2/mV/ms = 10^6 S/m^2)
g_k_ax = 4.0e-4	-- axon
g_k_so = 2.0e-4	-- soma
g_k_de = 3.0e-5	-- dendrite

g_na_ax = 3.0e-2
g_na_so = 1.5e-3
g_na_de = 4.0e-5

g_l_ax = 2.0e-4
g_l_so = 1.0e-6
g_l_de = 1.0e-6

-- capacitance (in units of C/mV/m^2 = 10^3 F/m^2)
spec_cap = 1.0e-5

-- resistance (in units of mV ms m / C = 10^-6 Ohm m)
spec_res = 1.5e6

-- reversal potentials (in units of mV)
e_k  = -90.0
e_na = 60.0
e_ca = 140.0

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

-- equilibrium potential (in units of mV)
v_eq = -65.0

-- diffusion coefficients (in units of m^2/ms)
diff_k 	= 1.0e-12
diff_na	= 1.0e-12
diff_ca	= 2.2e-13

-- temperature in units of deg Celsius
temp = 37.0

----------------------------------
-- setup approximation space	--
----------------------------------
-- create, load, refine and distribute domain
--[[
write(">> Loading, distributing and refining domain ...")

neededSubsets = {"Axon", "Dendrite", "Soma"}
dom = util.CreateDomain(gridName, numPreRefs, neededSubsets)

balancer.partitioner = "parmetis"
balancer.firstDistLvl = -1
balancer.redistProcs = 64
balancer.parallelElementThreshold = 24
balancer.maxLvlsWithoutRedist = 0
balancer.ParseParameters()
balancer.PrintParameters()

-- in parallel environments: use a load balancer to distribute the grid
loadBalancer = balancer.CreateLoadBalancer(dom)
balancer.RefineAndRebalanceDomain(dom, numRefs, loadBalancer)

print(dom:domain_info():to_string())

if load_balancer ~= nil then
	loadBalancer:print_quality_records()
end

write(">> done\n")
--]]
---[[
neededSubsets = {"Axon", "Dendrite", "Soma"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, "metis")
--]]

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
OrderCuthillMcKee(approxSpace, true);

----------------------
-- setup elem discs	--
----------------------

-- cable equation
VMD = CableEquation("Axon, Dendrite, Soma, PreSynapseEdges, PostSynapseEdges")
VMD:set_spec_cap(spec_cap)
VMD:set_spec_res(spec_res)

VMD:set_rev_pot_k(e_k)
VMD:set_rev_pot_na(e_na)
VMD:set_rev_pot_ca(e_ca)

VMD:set_k_out(k_out)
VMD:set_na_out(na_out)
VMD:set_ca_out(ca_out)

VMD:set_diff_coeffs({diff_k, diff_na, diff_ca})

VMD:set_temperature_celsius(temp)


-- Hodgkin and Huxley channels
--HH = ChannelHHNernst("v, k, na", "Axon")
HHaxon = ChannelHH("v", "Axon, PreSynapseEdges")
HHaxon:set_conductances(g_k_ax, g_na_ax)
HHsoma = ChannelHH("v", "Soma")
HHsoma:set_conductances(g_k_so, g_na_so)
HHdend = ChannelHH("v", "Dendrite, PostSynapseEdges")
HHdend:set_conductances(g_k_de, g_na_de)

VMD:add(HHaxon)
VMD:add(HHsoma)
VMD:add(HHdend)

-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leakAxon = ChannelLeak("v", "Axon, PreSynapseEdges")
leakAxon:set_cond(g_l_ax*tmp_fct)
leakAxon:set_rev_pot(-66.148458)
leakSoma = ChannelLeak("v", "Soma")
leakSoma:set_cond(g_l_so*tmp_fct)
leakSoma:set_rev_pot(-30.654022)
leakDend = ChannelLeak("v", "Dendrite, PostSynapseEdges")
leakDend:set_cond(g_l_de*tmp_fct)
leakDend:set_rev_pot(-57.803624)

VMD:add(leakAxon)
VMD:add(leakSoma)
VMD:add(leakDend)

-- synapses
syn_handler = NETISynapseHandler()
syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_vmdisc(VMD)
syn_handler:set_activation_timing(
	0.1,	-- average start time of synaptical activity in ms
	5,		-- average duration of activity in ms (10)
	1.0,	-- deviation of start time in ms
	0.5,	-- deviation of duration in ms
	1.2e-3)	-- peak conductivity in units of uS (6e-4)
VMD:set_synapse_handler(syn_handler)

--[[
-- treat unknowns on synapse subset
diri = DirichletBoundary()
diri:add(0.0, "v", "Exp2Syn")
diri:add(0.0, "k", "Exp2Syn")
diri:add(0.0, "na", "Exp2Syn")
diri:add(0.0, "ca", "Exp2Syn")
--]]

-- domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(VMD)
--domainDisc:add(diri)

-- time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)

-- create instationary operator
op = AssembledOperator(timeDisc)
op:init()


------------------
-- solver setup	--
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)


-- GMG --
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)

gmg:set_base_solver(LU())
gmg:set_gathered_base_solver_if_ambiguous(true)

gmg:set_smoother(ILU())
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)

--gmg:set_debug(dbgWriter)


-- linear solver --
linConvCheck = ConvCheck()
linConvCheck:set_maximum_steps(2000)
linConvCheck:set_minimum_defect(1e-50)
linConvCheck:set_reduction(1e-08)
linConvCheck:set_verbose(false)

cgSolver = CG()
cgSolver:set_preconditioner(ILU())
cgSolver:set_convergence_check(linConvCheck)
--cgSolver:set_debug(dbgWriter)

-- non-linear solver --
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 20, 2e-26, 1e-08)
newtonConvCheck:set_component_check("v", 1e-21, 1e-12)
newtonConvCheck:set_verbose(true)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(cgSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)


----------------------
-- time stepping	--
----------------------
newtonSolver:init(op)

time = 0.0

-- init solution
u = GridFunction(approxSpace)
u:set(0.0)
Interpolate(v_eq, u, "v")
Interpolate(k_in, u, "k");
Interpolate(na_in, u, "na");
Interpolate(ca_in, u, "ca")

-- write start solution
if (generateVTKoutput) then 
	out = VTKOutput()
	out:print(fileName .."vtk/Solvung", u, 0, time)
end

-- store grid function in vector of  old solutions
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

for step = 1,nSteps do
	print("++++++ POINT IN TIME " .. math.floor((time+dt)/dt+0.5)*dt .. " BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)

	-- prepare Newton solver
	newtonSolver:prepare(u)
		
	-- apply Newton solver
	if newtonSolver:apply(u) == false then
		print ("Newton solver apply failed at step "..step..".");
		if (generateVTKoutput) then 
			out:write_time_pvd(fileName .."vtk/Solvung", u);
		end
		exit();
	end 
	
	-- update to new time
	time = solTimeSeries:time(0) + dt
	
	-- vtk output
	if (generateVTKoutput) then
		if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
			out:print(fileName .."vtk/Solvung", u, math.floor(time/pstep+0.5), time)
		end
	end
	
	-- updte time series (reuse memory)
	oldestSol = solTimeSeries:oldest()
	VecScaleAssign(oldestSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldestSol, time)

	print("++++++ POINT IN TIME " .. math.floor((time)/dt+0.5)*dt .. "  END ++++++");
end

-- end timeseries, produce gathering file
if (generateVTKoutput) then
	out:write_time_pvd(fileName .."vtk/Solvung", u)
end

