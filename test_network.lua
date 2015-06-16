------------------------------------------------------
-- This script is intended for testing purposes.	--
-- It solves the cable equation with HH channels,	--
-- activating synapses and transmission synapses.	--
------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")

-- choice of grid
gridName = util.GetParam("-grid", "test.ugx")

-- dimension
ugxfi = UGXFileInfo()
ugxfi:parse_file("/Users/pgottmann/Documents/workspace/ug4/trunk/apps/cable/"..gridName)
dim = ugxfi:physical_grid_dimension(0)
print("Detected dimension "..dim.." in ugx file.\n")

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"SynapseHandler","HH_Kabelnew"})

-- parameters steering simulation
numPreRefs	= util.GetParamNumber("-numPreRefs",	0)
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			0.01) -- in ms
endTime		= util.GetParamNumber("-endTime",		10.0) -- in ms
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)

-- vtk output?
generateVTKoutput	= util.HasParamOption("-vtk")

-- file handling
fileName = util.GetParam("-outName", "Solvung")
fileName = fileName.."/"


--------------------------
-- biological settings	--
--------------------------
-- membrane conductances
g_Na = 1.2e-3		-- in C/m^2/mV/ms = 10^6 S/m^2
g_K  = 0.36e-3		-- in C/m^2/mV/ms = 10^6 S/m^2
g_L  = 0.003e-3		-- in C/m^2/mV/ms = 10^6 S/m^2

-- capacitance
spec_cap = 1.0e-5	-- in C/mV/m^2 = 10^3 F/m^2

-- resistance
spec_res = 1.0e6	-- in mV ms m / C = 10^-6 Ohm m

-- diameter
Diameter = 1.0e-6	-- in m

-- reversal potentials
ena = 63.5129		-- in mV
ek  = -74.1266		-- in mV

-- diffusion coefficients
diff_k 	= 1.0e-12	-- in m^2/ms
diff_na	= 1.0e-12	-- in m^2/ms
diff_ca	= 2.2e-13	-- in m^2/ms

-- accuracy for gating params
ac = 1e-6

----------------------------------
-- setup approximation space	--
----------------------------------

-- Create, Load, Refine and Distribute Domain
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
OrderCuthillMcKee(approxSpace, true);

----------------------
-- setup elem discs	--
----------------------

-- cable equation
VMD = VMDisc("Axon, Dendrite, Soma")
--VMD:set_diameter(Diameter)
--VMD:set_diff_coeffs({diff_k, diff_na, diff_ca})
--VMD:set_spec_cap(spec_cap)
--VMD:set_spec_res(spec_res)
--VMD:set_influx_ac(Flux_ac)

-- Hodgkin and Huxley channels
HH = ChannelHH("v", "Axon")
VMD:add_channel(HH)
--HH = ChannelHHNernst("v, k, na", "axon")

-- synapses
syn_handler = NETISynapseHandler()
syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_vmdisc(VMD)
syn_handler:set_activation_timing(
	0.1,	-- average start time of synaptical activity in ms
	5,		-- average duration of activity in ms (10)
	1.0,	-- deviation of start time in ms
	0.5,	-- deviation of duration in ms
	1.2e-12)	-- peak conductivity (6e-4)
VMD:set_synapse_handler(syn_handler)

-- treat unknowns on synapse subset
diri = DirichletBoundary()
diri:add(0.0, "v", "Exp2Synapses")
diri:add(0.0, "k", "Exp2Synapses")
diri:add(0.0, "na", "Exp2Synapses")
diri:add(0.0, "ca", "Exp2Synapses")


-- domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(VMD)
domainDisc:add(diri)

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
linConvCheck:set_reduction(1e-04)
linConvCheck:set_verbose(false)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(ILU())
bicgstabSolver:set_convergence_check(linConvCheck)
--bicgstabSolver:set_debug(dbgWriter)

-- non-linear solver --
newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 20, 2e-26, 1e-08)
newtonConvCheck:set_component_check("v", 1e-21, 1e-12)
newtonConvCheck:set_verbose(true)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
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
Interpolate(-63.8167, u, "v")
Interpolate(54.4, u, "k");
Interpolate(10.0, u, "na");
Interpolate(5e-5, u, "ca")

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
	timeDisc:prepare_step_elem(solTimeSeries, dt)
	
	-- update presynaptic Vm values (must be done AFTER prep_step_elem)
	syn_handler:update_presyn()

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
		out:print(fileName .."vtk/Solvung", u, step, time)
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

