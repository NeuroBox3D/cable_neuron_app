-------------------------------------------------------------
-- This script solves the cable equation with HH channels. --
-- The model is fully implicit on the axon, where the HH   --
-- channels are located. Activation through current        --
-- injection at the soma.                                  --
--                                                         --
-- author: mbreit                                          --
-- date:   2019-06-15                                      --
-------------------------------------------------------------
ug_load_script("ug_util.lua")

-- choice of grid
gridName = util.GetParam("-grid", "grids/31o_pyramidal19aFI.CNG_diams.ugx")

-- init UG
InitUG(3, AlgebraType("CPU", 1))
AssertPluginsLoaded({"cable_neuron"})


-- parameters steering simulation
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs", 0)
dt = util.GetParamNumber("-dt", 1e-5)  -- in s
endTime = util.GetParamNumber("-endTime", 0.1)  -- in s
nSteps = util.GetParamNumber("-nSteps", endTime/dt)

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")

-- file handling
fileName = util.GetParam("-outName", "solution")
fileName = fileName.."/"


-------------------------
-- biological settings --
-------------------------
-- membrane conductances
g_Na = 1.2e3		-- in S/m^2
g_K  = 360.0		-- in S/m^2
g_L  = 3.0			-- in S/m^2

-- specific capacitance
spec_cap = 1.0e-2	-- in F/m^2

-- resistivity
spec_res = 1.0		-- in Ohm m

-- reversal potentials
e_na = 0.050		-- in V
e_k = -0.077		-- in V
e_l = -0.0544       -- in V

-- diffusion coefficients
diff_k 	= 1.0e-9	-- in m^2/s
diff_na	= 1.0e-9	-- in m^2/s
diff_ca	= 2.2e-10	-- in m^2/s

-- equilibrium potential
v_eq = -0.065       -- in V

-------------------------------
-- setup approximation space --
-------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"axon", "dend", "soma"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, "metis")

-- create Approximation Space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("h", "Lagrange", 1)
approxSpace:add_fct("m", "Lagrange", 1)
approxSpace:add_fct("n", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
approxSpace:add_fct("ca", "Lagrange", 1)

-- gating functions for HH-Fluxes
approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true)

----------------------
-- setup elem discs --
----------------------
influxCoordsx = 0.0 -- in m
influxCoordsy = 0.0 -- in m
influxCoordsz = 0.0
injectionDensity = 5.0  -- in A/m^2

function injection3d(t, x, y, z)
	if t >= 0 and t <= 0.001
		and math.abs(x - influxCoordsx) < 1e-8
		and math.abs(y - influxCoordsy) < 1e-8
		and math.abs(z - influxCoordsz) < 1e-8
	then
		return injectionDensity
	end
	return 0.0
end

injection_flux = LuaFunctionNumber()
injection_flux:set_lua_callback("injection3d", 4)

HH = ImplicitActiveCableDiscNernst("v, h, m, n, k, na", "axon")
HH:set_injection(injection_flux)
HH:set_spec_cap(spec_cap)
HH:set_spec_res(spec_res)
HH:set_conductances(g_K, g_Na, g_L)
HH:set_diffusion_constants(diff_k, diff_na)
HH:set_rev_pot(e_k, e_na, e_l)

CE = CableEquation("dend, soma", true)

diri = DirichletBoundary()
diri:add(0.0, "h", "dend, soma")
diri:add(0.0, "n", "dend, soma")
diri:add(0.0, "m", "dend, soma")
diri:add(5e-5, "ca", "axon")


---------------------------------
-- setup domain discretization --
---------------------------------
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)
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

-- set initial value
print("Interpolation start values")

-- set initial values for Na and K
Interpolate(12, u, "na", time)
Interpolate(155.0, u, "k", time)
Interpolate(5e-5, u, "ca", time)
Interpolate(v_eq, u, "v", time)

AlphaHh = 70.0*math.exp((-(v_eq+0.065))/0.020)
BetaHh = 1e3/(math.exp(-100.0*(v_eq+0.035)) + 1.0)
h_eq = AlphaHh/(AlphaHh+BetaHh)
Interpolate(h_eq, u, "h", time)

AlphaHm = -1e5*(v_eq+0.040)/(math.exp(-100.0*(v_eq+0.040)) - 1.0)
BetaHm = 4e3*math.exp(-(v_eq+0.065)/0.018)
m_eq = AlphaHm/(AlphaHm+BetaHm)
Interpolate(m_eq, u, "m", time)

AlphaHn = -1e4*(v_eq+0.055)/(math.exp(-100.0*(v_eq+0.055)) - 1.0)
BetaHn = 125.0*math.exp(-(v_eq+0.065)/0.080)
n_eq = AlphaHn/(AlphaHn+BetaHn)
Interpolate(n_eq, u, "n", time)


-- write start solution
if generateVTKoutput then 
	out = VTKOutput()
	out:print(fileName .."vtk/solution", u, 0, time)
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
		if generateVTKoutput then 
			out:write_time_pvd(fileName .."vtk/solution", u);
		end
		exit();
	end 
	
	-- update to new time
	time = solTimeSeries:time(0) + dt
	
	-- vtk output
	if generateVTKoutput then 
		out:print(fileName .."vtk/solution", u, step, time)
	end
	
	
	-- updte time series (reuse memory)
	oldestSol = solTimeSeries:oldest()
	VecScaleAssign(oldestSol, 1.0, u)
	solTimeSeries:push_discard_oldest(oldestSol, time)

	print("++++++ POINT IN TIME " .. math.floor((time)/dt+0.5)*dt .. "  END ++++++");
end

-- end timeseries, produce gathering file
if generateVTKoutput then
	out:write_time_pvd(fileName .."vtk/solution", u)
end

