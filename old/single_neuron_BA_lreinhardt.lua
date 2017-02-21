--------------------------------------------------------------
-- This script solves the cable equation with HH channels, 	--
-- activating synapses and transmission synapses.			--
--------------------------------------------------------------


-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")


--------------------------------------------------------------------------------
-- UG4-Standard-Settings
--------------------------------------------------------------------------------

-- choice of grid
--gridName = util.GetParam("-grid", "grids/test_cell_small.ugx")
gridName = util.GetParam("-grid", "grids/31o_pyramidal19aFI.CNG.ugx")
--gridName = util.GetParam("-grid", "grids/31o_pyramidal19aFI.CNG_diams.ugx")

-- dimension
ugxfi = UGXFileInfo()
ugxfi:parse_file(gridName)
dim = ugxfi:physical_grid_dimension(0)
print("Detected dimension "..dim.." in ugx file.\n")

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"SynapseDistributor", "SynapseHandler","HH_Kabelnew"})

-- parameters steering simulation
numPreRefs	= util.GetParamNumber("-numPreRefs",	0)
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			1e-5) -- in units of s
endTime		= util.GetParamNumber("-endTime",		0.015) -- in units of s
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)

print(" chosen parameters:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

-- SynapseDistributor activity parameters
num_synapses = util.GetParamNumber("-nSyn", 1000)
avg_start = util.GetParamNumber("-avgStart", 0.0)
avg_dur = util.GetParamNumber("-avgDur", 15)
dev_start = util.GetParamNumber("-devStart", 0.0)
dev_dur = util.GetParamNumber("-devDur", 0.0)

-- vtk output?
generateVTKoutput	= util.HasParamOption("-vtk")

-- file handling
filename = util.GetParam("-outName", "sol")
filename = filename .. "/sol"


--------------------------
-- biological settings	--
--------------------------
-- membrane conductances
g_Na = 1.2e3		-- in S/m^2
g_K  = 360.0		-- in S/m^2
g_L  = 3.0			-- in S/m^2

-- specific capacitance (in units of F/m^2)
spec_cap = 1.0e-2

-- specific resistance
spec_res = 1.0		-- in Ohm m

-- diameter
diameter = 1.0e-6	-- in m

-- reversal potentials
ena = 0.05	--0.0635129		-- in V
ek  = -0.077	---0.0741266		-- in V

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.0e-9
diff_na	= 1.0e-9
diff_ca	= 2.2e-10


--------------------------------------------------------------------------------
-- Create, Load, Refine Domain
--------------------------------------------------------------------------------
--neededSubsets = {"Inner"}
neededSubsets = {"axon", "dend", "soma"}
dom = util.CreateDomain(gridName, numRefs, neededSubsets)
--dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, "metis")


--------------------------------------------------------------------------------
-- Synapse distributions via plugin by Lukas Reinhardt
--------------------------------------------------------------------------------

sd = SynapseDistributor(dom, "grids/grid_out.ugx", true)

-- PARAMETER SET FOR ACTIVITY TIMING
--[[
15,					--average start time of synaptical activity in ms
1,					--average duration of activity in ms
13, 				--deviation of start time in ms
0.5 				--deviation of duration in ms
--]]

--sd:place_synapses_uniform(100)
--sd:place_synapses_uniform(num_synapses)
sd:place_synapses_uniform(1, num_synapses)
sd:set_activation_timing(avg_start, avg_dur, dev_start, dev_dur)
sd:print_status()


--------------------------------------------------------------------------------
-- Distribute Domain
--------------------------------------------------------------------------------
util.DistributeDomain(dom, "metis")

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
--VMDisc constructor creates every needed concentration out of added Channels from Channel list
--------------------------------------------------------------------------------

-- Hodgkin and Huxley channels
HH = ChannelHHNernst("v, k, na", "axon")
--HH = ChannelHH("v", "axon")

-- leakage
--leakAxon = ChannelLeak("v", "axon")
--leakAxon:set_rev_pot(-0.0544)
--leakDend = ChannelLeak("v", "dend")

-- synapses
--[[
syn_handler = NETISynapseHandler()
syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_ce_object(VMD)
syn_handler:set_activation_timing(
	0.1,	-- average start time of synaptical activity in ms
	5,		-- average duration of activity in ms (10)
	1.0,	-- deviation of start time in ms
	0.5,	-- deviation of duration in ms
	1.2e-3)	-- peak conductivity (6e-4)
VMD:set_synapse_handler(syn_handler)
--]]

-- cable equation
--VMD = VMDisc("Inner")
VMD = VMDisc("axon, dend, soma")
--VMD:set_diameter(diameter)
--VMD:set_diff_coeffs({diff_k, diff_na, diff_ca})
--VMD:set_spec_cap(spec_cap)
--VMD:set_spec_res(spec_res)
--VMD:set_influx_ac(Flux_ac)
VMD:set_synapse_distributor(sd)
VMD:add_channel(HH)
--VMD:add_channel(leakAxon)
--VMD:add_channel(leakDend)


--------------------------------------------------------------------------------
--	ELECTRODE STIMULATION SETUP
--------------------------------------------------------------------------------

--  INFO: coords for 31o_pyramidal19aFI.CNG_with_subsets.ugx

--	ELECTRODE STIMULATION near soma: 5pA seem to inervate the pyramidal cell with uniform diameters of 1um
--VMD:set_influx(5e-9, 6.54e-05, 2.665e-05, 3.985e-05, 0.0, 0.04) -- current given in units of A

--	!!! DO NOT USE!!! ELECTRODE STIMULATION into soma: 5pA seem to enervate the pyramidal cell with attached diameters (also causing backpropagating APs)
--VMD:set_influx(5e-9, 6.9e-06, 3.74e-05, -2.86e-05, 0.0, 0.04) -- current given in units of A


-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(VMD)


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
op = AssembledOperator(timeDisc)
op:init()


-------------------------------------------
-- solver setup
-------------------------------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(1.0)

gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
--bgs = BackwardGaussSeidel()
ilu = ILU()

lu = LU()
--slu = SuperLU()
--ilu:set_debug(dbgWriter)
--ilut = ILUT()

-- create GMG ---

-- Base Solver
--baseConvCheck = ConvCheck()
--baseConvCheck:set_maximum_steps(500)
--baseConvCheck:set_minimum_defect(1e-8)
--baseConvCheck:set_reduction(1e-30)
--baseConvCheck:set_verbose(false)
base = LU()
--base = LinearSolver()
--base:set_convergence_check(baseConvCheck)
--base:set_preconditioner(jac)
--base = sgs
-- Geometric Multi Grid
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_base_solver(base)
gmg:set_smoother(ilu)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
--gmg:set_debug(dbgWriter)

-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(2000)
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-08)
convCheck:set_verbose(false)

-- create Linear Solver
--linSolver = LinearSolver()
--linSolver:set_preconditioner(ILU())
--linSolver:set_convergence_check(convCheck)

-- create ILU Solver
iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilu)
iluSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(ilu)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

-- linear iterative solver
linSolver = LinearSolver()
--linSolver:set_preconditioner(slu)
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)


-- convergence check
--[[
newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(30)
newtonConvCheck:set_minimum_defect(1e-21)
newtonConvCheck:set_reduction(1e-10)
]]--

--newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 20, 2e-26, 1e-08)
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 20, 2e-24, 1e-08)
newtonConvCheck:set_component_check("v", 1e-21, 1e-12)
--newtonConvCheck:set_component_check("Na", 1e-15, 1e-10)
--newtonConvCheck:set_component_check("K", 1e-15, 1e-10)
newtonConvCheck:set_verbose(true)
newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
--newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)
--newtonSolver:set_line_search(newtonLineSearch)


-----------------------
-- Start time stepping
-----------------------
--newtonSolver:set_linear_solver(slu)
newtonSolver:init(op)
time = 0.0
step = 0

------------------------------------------------------------
-- creating grid function
------------------------------------------------------------
-- get grid function
u = GridFunction(approxSpace)
u:set(0.0)

-------------------------------------
--Setting all Startvalues
-------------------------------------

-- set initial value
Interpolate(-0.065, u, "v", time)
Interpolate(54.4, u, "k", time);
Interpolate(10.0, u, "na", time);
Interpolate(5e-5, u, "ca", time)

-- write start solution
if (generateVTKoutput) then 
	out = VTKOutput()
	out:print(filename, u, step, time)
end

-- store grid function in vector of  old solutions
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

for step = 1, nSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)

	-- prepare newton solver
	if newtonSolver:prepare(u) == false then 
		print ("Newton solver prepare failed at step "..step..".")
		if (generateVTKoutput) then 
			out:write_time_pvd(filename, u) 
		end
		exit()
	end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then 
		print ("Newton solver apply failed at step "..step..".")
		if (generateVTKoutput) then 
			out:write_time_pvd(filename, u)
		end
		exit()
	end 
	
	-- update new time
	time = solTimeSeries:time(0) + dt
	
	if (generateVTKoutput) then 
		out:print(filename, u, step, time)
	end
	
	oldestSol = solTimeSeries:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)

	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file
if (generateVTKoutput) then 
	out:write_time_pvd(filename, u) 
end
	
