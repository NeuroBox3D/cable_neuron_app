--------------------------------------------------------------
-- This script is supposed to serve as reference to compare --
-- a full 3d simulation against (for debugging purposes).   --
-- It solves the cable equation with HH channels and two    --
-- electrodes.			                                    --
--                                                          --
-- author: mbreit                                           --
-- date:   2017-10-28                                       --
--------------------------------------------------------------

ug_load_script("ug_util.lua")

-- init UG
InitUG(3, AlgebraType("CPU", 1));
AssertPluginsLoaded({"cable_neuron"})


--------------
-- settings --
--------------
gridName = "../neuron/grids/grid1d.ugx"
neededSubsets = {"soma", "dend", "dend_inj1", "dend_inj2", "axon", "axon_myel"}

-- parameters steering simulation
numRefs = util.GetParamNumber("-numRefs", 0)
dt = util.GetParamNumber("-dt", 1e-5) -- in s
endTime = util.GetParamNumber("-endTime", 1.0)  -- in s
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")

-- file handling
filename = util.GetParam("-outName", "test/neuron1d")
filename = filename .. "/"


-------------------------
-- biological settings --
-------------------------
-- membrane conductances (in units of S/m^2)
g_k = 360.0
g_na = 1.2e3
g_l = 3.0

-- specific capacitance (in units of F/m^2)
spec_cap = 3.54e-3

-- resistivity (in units of Ohm m)
spec_res = 2.0

-- reversal potentials (in units of V)
e_k  = -0.077
e_na = 0.05

-- equilibrium potential (in units of V)
v_eq = -0.065

-- temperature in units of deg Celsius
temp = 37.0


-------------------------
-- create approx space --
-------------------------
dom = util.CreateDomain(gridName, numRefs, neededSubsets)

-- scale domain from um to m
scale_domain(dom, 1e-6) 

-- check domain is acyclic
isAcyclic = is_acyclic(dom)
if not isAcyclic then
	print("Domain is not acyclic!")
	exit()
end

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)

approxSpace:init_levels();
approxSpace:init_top_surface();
OrderCuthillMcKee(approxSpace, true);


--------------------
-- discretization --
--------------------
-- cable equation
CE = CableEquation("soma, dend, dend_inj1, dend_inj2, axon, axon_myel", false)

CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)

CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)

CE:set_temperature_celsius(temp)

-- Hodgkin and Huxley channels
HH = ChannelHH("v", "axon, soma")
HH:set_conductances(g_k, g_na, "axon, soma")
HH:enable_temperature_dependency(false)
CE:add(HH)

-- leakage
leak = ChannelLeak("v", "soma, dend, dend_inj1, dend_inj2, axon")
leak:set_cond(g_l)
leak:set_rev_pot(-0.0546)
CE:add(leak)

-- electrode stimulation
CE:set_influx_subset(4, 1.5, 0.006, 0.0) -- dend_inj1: 1.5 A/(m^2*s)
CE:set_influx_subset(5, 1.5, 0.006, 0.0) -- dend_inj1: 1.5 A/(m^2*s)


-- create domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)

-- get assembling tuner for speedup later
assTuner = domainDisc:ass_tuner()

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)

-- create operator from discretization
linOp = AssembledLinearOperator(timeDisc)


------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- linear solver --
convCheck = CompositeConvCheck(approxSpace, 2, 2e-26, 1e-08)
convCheck:set_component_check("v", 1e-21, 1e-12)
convCheck:set_verbose(verbose)

ilu = ILU()
solver = LinearSolver()
solver:set_preconditioner(ilu)
solver:set_convergence_check(convCheck)
--solver:set_debug(dbgWriter)


-------------------
-- time stepping --
-------------------
time = 0.0

-- init solution
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
Interpolate(v_eq, u, "v")

-- write start solution
if generateVTKoutput then 
	out = VTKOutput()
	out:print(filename.."vtk/solution", u, 0, time)
end

-- store grid function in vector of old solutions
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

curr_dt = dt
dtred = 2

lv = 0
maxLv = 10
cb_counter = {}
cb_counter[lv] = 0

while endTime-time > 0.001*curr_dt do
	-- prepare next step
	timeDisc:prepare_step(solTimeSeries, curr_dt)
	
	-- reduce time step if cfl < curr_dt
	-- (this needs to be done AFTER prepare_step as channels are updated there)
	dtChanged = false
	cfl = CE:estimate_cfl_cond(solTimeSeries:latest())
	print("estimated CFL condition: dt < " .. cfl)
	while (curr_dt > cfl) do
		curr_dt = curr_dt/dtred
		
		if lv+1 > maxLv then
			print("Time step too small.")
			exit()
		end
		
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
		
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(linOp, u, b, solver) == false then
		print("Could not apply linear solver.")
		exit()
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


	
