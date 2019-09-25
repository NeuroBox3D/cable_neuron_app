--------------------------------------------------------------------------------
-- This script solves the cable equation with HH channels and leakage.        --
-- Activation is realized by randomly distributed synapses.                   --
--                                                                            --
-- Authors: Markus Breit, Pascal Gottmann                                     --
-- Date:    2015-08-19                                                        --
--------------------------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"cable_neuron"})

-- init UG
InitUG(3, AlgebraType("CPU", 1))


---------------------------------
-- read command line arguments --
---------------------------------
-- choice of grid
gridName = util.GetParam("-grid", "cable_neuron_app/grids/13-L3pyr-77.CNG.ugx")
gridSyn = string.sub(gridName, 1, string.len(gridName)-4) .. "_syn.ugx"

-- parameters steering simulation
numRefs = util.GetParamNumber("-numRefs", 0)
dt = util.GetParamNumber("-dt", 1e-5)  -- in s
endTime = util.GetParamNumber("-endTime", dt)
nSteps = util.GetParamNumber("-nSteps", endTime/dt)
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")

-- synapse activity parameters
avg_start = util.GetParamNumber("-avgStart", 0.003)
avg_dur = util.GetParamNumber("-avgDur", 2.4e-4)
dev_start = util.GetParamNumber("-devStart", 0.001)
dev_dur = util.GetParamNumber("-devDur", 0.0)
num_synapses = util.GetParamNumber("-nSyn", 750)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")

-- file handling
outPath = util.GetParam("-outName", "solution")
outPath = outPath.."/"


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
v_eq = -0.07

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.0e-9
diff_na	= 1.0e-9
diff_ca	= 2.2e-10

-- temperature in units of deg Celsius
temp = 37.0


------------------------------------
-- create domain and approx space --
------------------------------------
-- synapse distribution
synDistr = SynapseDistributor(gridName)
synDistr:clear() -- clear any synapses from grid

-- place half of the synapses on subets 1 and 2 ("dendrite" and "apical_dendrite") each
synDistr:place_synapses({0.0, 0.0, 0.5, 0.5}, num_synapses, "AlphaPostSynapse")
if not synDistr:export_grid(gridSyn) then
	print("SynapseDistributor grid export unsuccessful. Aborting.")
	exit()
end

gridName = gridSyn


-- collect functional subset groups
-- this has to be adapted according to the geometry used
somaSubsets = {"soma"}
dendSubsets = {"dendrite", "apical_dendrite"}
axonSubsets = {"axon"}


allSubsets = {}
allSubsetsString = ""
for _, v in pairs(somaSubsets) do
    table.insert(allSubsets, v)
    allSubsetsString = allSubsetsString .. ", " .. v
end
for _, v in pairs(dendSubsets) do
    table.insert(allSubsets, v)
    allSubsetsString = allSubsetsString .. ", " .. v
end
for _, v in pairs(axonSubsets) do
    table.insert(allSubsets, v)
    allSubsetsString = allSubsetsString .. ", " .. v
end
if allSubsetsString:len() > 2 then
	allSubsetsString = allSubsetsString:sub(3)
end

dom = util.CreateDomain(gridName, numRefs, allSubsets)

-- check domain is acyclic
isAcyclic = is_acyclic(dom)
if not isAcyclic then
	print("Domain is not acyclic!")
	exit()
end


approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)

approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true)


--------------------
-- discretization --
--------------------
-- cable equation
CE = CableEquation(allSubsetsString, false)

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
HH = ChannelHH("v", allSubsetsString)
if #axonSubsets > 0 then HH:set_conductances(g_k_ax, g_na_ax, axonSubsets) end
if #somaSubsets > 0 then HH:set_conductances(g_k_so, g_na_so, somaSubsets) end
if #dendSubsets > 0 then HH:set_conductances(g_k_de, g_na_de, dendSubsets) end

CE:add(HH)

-- leakage (exactly calibrated to achieve zero net current in equilibrium)
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leak = ChannelLeak("v", allSubsetsString)
if #axonSubsets > 0 then leak:set_cond(g_l_ax*tmp_fct, axonSubsets) end
if #axonSubsets > 0 then leak:set_rev_pot(-0.070212, axonSubsets) end
if #somaSubsets > 0 then leak:set_cond(g_l_so*tmp_fct, somaSubsets) end
if #somaSubsets > 0 then leak:set_rev_pot(-0.059236, somaSubsets) end
if #dendSubsets > 0 then leak:set_cond(g_l_de*tmp_fct, dendSubsets) end
if #dendSubsets > 0 then leak:set_rev_pot(-0.067947, dendSubsets) end

CE:add(leak)


-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing_alpha(
	avg_start,   -- average onset of synaptical activity in [s]
	avg_dur/6.0, -- average tau of activity function in [s]
	dev_start,   -- deviation of onset in [s]
	dev_dur/6.0, -- deviation of tau in [s]
	1.2e-09)      -- peak conductivity in [S]
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
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- linear solver --
linConvCheck = CompositeConvCheck(approxSpace, 1, 2e-26, 1e-08)
linConvCheck:set_component_check("v", 1e-21, 1e-12)
linConvCheck:set_verbose(verbose)

ilu = ILU()
linSolver = LinearSolver()
linSolver:set_preconditioner(ilu)
linSolver:set_convergence_check(linConvCheck)
--linSolver:set_debug(dbgWriter)


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
	out:print(outPath.."vtk/solution", u, 0, time)
end

--[[
-- measurement setup
measFileVm = outPath.."measVm.txt"
if ProcRank() == 0 then
	measOutVm = assert(io.open(measFileVm, "a"))
end
--]]

-- store grid function in vector of  old solutions
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
	-- setup time Disc for old solutions and timestep
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
	
	-- synchronize (for profiling)
	PclDebugBarrierAll()
	
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(linOp, u, b, linSolver) == false then
		print("Could not apply linear solver.")
		exit()
	end
	
	--[[
	-- log Vm and Ca
	if ProcRank() == 0 then
		vm_soma  = EvaluateAtClosestVertex(MakeVec(6.9e-07, 3.74e-06, -2.86e-06), u, "v", "soma", dom:subset_handler())
		vm_apic  = EvaluateAtClosestVertex(MakeVec(-4.05e-06, 6.736e-05, -1.341e-05), u, "v", "apical_dendrite", dom:subset_handler())
		vm_dend  = EvaluateAtClosestVertex(MakeVec(-4.631e-05, -0.0001252, 4.62e-06), u, "v", "dendrite", dom:subset_handler())
		measOutVm:write(time, "\t", vm_soma, "\t", vm_apic, "\t", vm_dend, "\t", -65, "\n")
	end
	--]]
	
	-- update to new time
	time = solTimeSeries:time(0) + curr_dt
	
	-- vtk output
	if generateVTKoutput then
		if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
			out:print(outPath.."vtk/solution", u, math.floor(time/pstep+0.5), time)
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
if generateVTKoutput then 
	out:write_time_pvd(outPath.."vtk/solution", u) 
end

--[[
-- close measure file
if ProcRank() == 0 then
	measOutVm:close()
end
--]]
	
