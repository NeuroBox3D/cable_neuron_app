------------------------------------------------------
-- This script is intended for testing purposes.	--
-- It solves the cable equation with HH channels,	--
-- activating synapses and transmission synapses.	--
------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")

-- choice of grid
gridName = util.GetParam("-grid", "../apps/cable/Ca_dyms/grids/31o_pyramidal19aFI.CNG_diams.ugx")

-- dimension
dim = 3

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"SynapseHandler","HH_Kabelnew"})

-- parameters steering simulation
numPreRefs	= util.GetParamNumber("-numPreRefs",	0)
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			0.01) -- in ms
endTime		= util.GetParamNumber("-endTime",		100.0) -- in ms
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
diameter = 2.0e-7	-- in m

-- reversal potentials
ena = 50.0			-- in mV
ek  = -77.0			-- in mV

-- diffusion coefficients
diff_k 	= 1.0e-12	-- in m^2/ms
diff_na	= 1.0e-12	-- in m^2/ms
diff_ca	= 2.2e-13	-- in m^2/ms

----------------------------------
-- setup approximation space	--
----------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"axon", "dend", "soma"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, "metis")

--print("Saving parallel grid layout")
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 1e-5)

-- create Approximation Space
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("h", "Lagrange", 1)
approxSpace:add_fct("m", "Lagrange", 1)
approxSpace:add_fct("n", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
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
InfluxPlacex = 6.9e-06 -- in m
InfluxPlacey = 3.74e-05 -- in m
InfluxPlacez = -2.86e-05
flux_ac = 1e-5
InfluxValue = 1e-14

function injection3d(t, x, y, z)
	--print("x: "..x.."y: "..y.."z: "..z.." time: "..time)
	if  t >= 0
		and t <= 20
		and math.abs(x - InfluxPlacex) < flux_ac
		and math.abs(y - InfluxPlacey) < flux_ac
		and math.abs(z - InfluxPlacez) < flux_ac
	then
	--print("flux is done")
	return InfluxValue
	else
	return 0.0
	end
end

injection_flux = LuaFunctionNumber()
injection_flux:set_lua_callback("injection"..dim.."d", (1+dim))

HH = ElemDiscHH_Nernst_FV1(approxSpace,"v, h, m, n, na, k", "axon")
HH:set_injection(injection_flux)
HH:set_diameter(diameter)
HH:set_spec_capa(spec_cap)
HH:set_spec_res(spec_res)
HH:set_consts(g_Na, g_K, g_L)
HH:set_accuracy(1e-6)
HH:set_nernst_consts(8.3136, 310.0, 96485.0)
HH:set_diff_Na(diff_na)
HH:set_diff_K(diff_k)

VMD = VMDisc("dend, soma")
--VMD:set_diameter(diameter)

diri = DirichletBoundary()
diri:add(0.0, "h", "dend, soma")
diri:add(0.0, "n", "dend, soma")
diri:add(0.0, "m", "dend, soma")
diri:add(0.0, "ca", "axon")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(VMD)
domainDisc:add(HH)
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


convCheckL = CompositeConvCheck(approxSpace, 2000, 1e-30, 1e-1)
convCheckL:set_component_check("v", 1e-30, 1e-04)
convCheckL:set_component_check("h,m,n", 5e-20, 1e-04)
convCheckL:set_component_check("na,k,ca", 5e-40, 1e-04)
convCheckL:set_verbose(false)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(ILU())
bicgstabSolver:set_convergence_check(convCheckL)
--bicgstabSolver:set_debug(dbgWriter)

newtonConvCheck = CompositeConvCheck(approxSpace, 20, 1e-15, 1e-10)
newtonConvCheck:set_component_check("v", 5e-27, 1e-08)
newtonConvCheck:set_component_check("h,m,n", 5e-15, 1e-10)
newtonConvCheck:set_component_check("na,k,ca", 1e-26, 1e-10)
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
-------------------------------------
--Setting all Startvalues
-------------------------------------

-- set initial value
print("Interpolation start values")

-- set initial values for Na and K
StartValue = -65.0
Interpolate(10, u, "na", time)
Interpolate(140.0, u, "k", time)
Interpolate(5e-5, u, "ca", time)
Interpolate(StartValue, u, "v", time)

-- Startvalues independent from VM
AlphaHn_test = (math.exp(1-0.1*(StartValue+65))-1);
if (math.abs(AlphaHn_test) > 1e-5) then
	AlphaHn = (0.1-0.01*(StartValue+65))/(math.exp(1-0.1*(StartValue+65))-1)
	else
	AlphaHn = 0.1
end
--print("StartAlphaHn: ".. AlphaHn)
BetaHn = 0.125*math.exp(-(StartValue+65)/80);

AlphaHm_test = math.exp(2.5-0.1*(StartValue+65))-1.0;
if (math.abs(AlphaHm_test) > 1e-5) then
	AlphaHm = (2.5 - 0.1*(StartValue+65)) / AlphaHm_test;
else
	AlphaHm = 1.0;
end
--AlphaHm = (2.5 - 0.1*(StartValue+65))/(math.exp(2.5-0.1*(StartValue+65))-1);
BetaHm = 4.0*math.exp(-(StartValue+65)/18);

AlphaHh = 0.07*math.exp((-1*(StartValue+65))/20);
BetaHh = 1/(math.exp(3-0.1*(StartValue+65))+1);

-- initialisiert H-Werte auf Startwert
h = AlphaHh/(AlphaHh+BetaHh);
startValueh = ConstUserNumber(0.596113)--(0.596113)--(h)--
Interpolate(h, u, "h", time)

-- initialisiert M-Werte auf Startwert
m = AlphaHm/(AlphaHm+BetaHm);
startValuem = ConstUserNumber(0.0529348)--(0.0529354)--(m/10)
Interpolate(m, u, "m", time)

-- initialisiert N-Werte auf Startwert
n = AlphaHn/(AlphaHn+BetaHn);
startValuen = ConstUserNumber(0.317681)--(0.31768)--(n)
Interpolate(n, u, "n", time)



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

