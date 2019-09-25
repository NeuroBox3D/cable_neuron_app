--------------------------------------------------------------
-- This script solves the cable equation with physiology    --
-- defined by the bAP model by Golding (2001)               --
-- The potential at a specified point is written to file.   --
--                                                          --
-- author: mbreit                                           --
-- date:   2019-06-13                                       --
--------------------------------------------------------------
ug_load_script("ug_util.lua")
AssertPluginsLoaded({"cable_neuron"})

-- init UG
InitUG(3, AlgebraType("CPU", 1));


---------------------------------
-- read command line arguments --
---------------------------------
-- choice of grid
gridName = util.GetParam("-grid", "cable_neuron_app/grids/vlachos_mouse_cell1.ugx")

-- cell model
cellModel = util.GetParamNumber("-cellModel", 0)

-- equilibration time
equil_time = util.GetParamNumber("-equilTime", 0.05)

-- parameters steering simulation
numRefs = util.GetParamNumber("-numRefs", 0)
dt = util.GetParamNumber("-dt", 1e-5) -- in s
simEndTime = util.GetParamNumber("-simEndTime", equil_time + 0.01)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")

-- file handling
outputPath = util.GetParam("-outName", ".")
outputPath = outputPath.."/"


-------------------------
-- biological settings --
-------------------------
-- settings are according to T. Branco

-- specific capacitance (in units of F/m^2)
spec_cap = 7.5e-3
spec_cap_myelin = 7.5e-4

-- resistivity (in units of Ohm m)
spec_res = 2.0

-- reversal potentials (in units of V)
e_k  = -0.09
e_na = 0.055
e_ca = 0.14
e_l = -0.066

-- membrane conductances (in units of S/m^2)
g_l = 0.25
g_l_node = 200.0

if model_cell == 0 or model_cell == 1 then
	g_k_a_p = 1e3  -- proximal A-type potassium starting density
	g_k_a_d = 1e3  -- distal A-type potassium starting density
else
	g_k_a_p = 1.3e3  -- proximal A-type potassium starting density
	g_k_a_d = 1.3e3  -- distal A-type potassium  starting density
end

g_k_dr = 400.0  -- delayed rectifier density

if cellModel == 0 then
	g_na = 420.0  -- sodium conductance 0.042 (S/m^2)
	naSlope = 0.0  -- slope of sodium channel density
elseif cellModel == 1 then
	g_na = 300.0  -- sodium conductance 0.030 (S/m^2)
	naSlope = 5e3  -- slope of sodium channel density (1/m)
end
g_na_node = 5e5  -- Na conductance at a node (S/m^2)

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



initSegFactor = 100
--initSegFrac = 0.8
naLimit = 2.75e-4  -- cut-off limit for increase of sodium conductance
function g_na_fct(x, y, z, t, si)
	if si == 2 then  -- axon_init
		-- has to be set individually for every geometry
		if y < -2.234e-5 then  -- lower fifth of initial segment
			return initSegFactor * g_na
		end
	elseif si == 4 then  -- ranvier
		return g_na_node
	elseif si == 6 or si == 7 then  -- apical
		local dist = math.sqrt(x*x + y*y + z*z)
		return (1.0 + math.min(dist, naLimit) * naSlope) * g_na
	end
	return g_na
end

function k_a_is_proximal(x, y, z, t, si)
	local dist = math.sqrt(x*x + y*y + z*z)
	if dist <= kProx then
		return 1
	end
	return 0
end

kLimit = 3e-4  -- cut-off for increase of K A-type density
kProx = 1e-4  -- distance to switch from proximal to distal type
kSlope = 1e4  -- slope of A-type density (1/m)
function g_k_a_fct(x, y, z, t, si)
	if si == 3 or si == 4 then  -- myelin/ranvier
		return 0.2 * g_k_a_p
	elseif si == 6 or si == 7 then  -- apical
		local dist = math.sqrt(x*x + y*y + z*z)
		if dist > kProx then
			return (1.0 + math.min(dist, kLimit) * kSlope) * g_k_a_d
		else
			return (1.0 + math.min(dist, kLimit) * kSlope) * g_k_a_p
		end
	end
	return g_k_a_p
end

spineLimit = 1e-4  -- distance beyond which to modify for spines
spineFactor = 2.0  -- factor by which to change passive properties
function g_l_fct(x, y, z, t, si)
	if si == 4 then  -- ranvier
		return g_l_node
	elseif si == 6 or si == 7 then  -- apical
		local dist = math.sqrt(x*x + y*y + z*z)
		if dist > spineLimit then
			return spineFactor * g_l
		end
	end
	return g_l
end

function cm_fct(x, y, z, t, si)
	if si == 3 then  -- myelin
		return spec_cap_myelin
	elseif si == 6 or si == 7 then  -- apical
		local dist = math.sqrt(x*x + y*y + z*z)
		if dist > spineLimit then
			return spineFactor * spec_cap
		end 
	end
	return spec_cap
end


------------------------------------
-- create domain and approx space --
------------------------------------
dom = Domain()
requiredSubsets = {"soma", "hillock", "axon_init", "myelin", "ranvier", "basal", "apic_rad", "apic_lm"}
dom = util.CreateDomain(gridName, numRefs, requiredSubsets)

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
allSubsetString = "soma, hillock, axon_init, myelin, ranvier, basal, apic_rad, apic_lm"
CE = CableEquation(allSubsetString, false)

CE:set_spec_cap("cm_fct")
CE:set_spec_res(spec_res)

CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)
CE:set_rev_pot_ca(e_ca)

CE:set_k_out(k_out)
CE:set_na_out(na_out)
CE:set_ca_out(ca_out)

CE:set_diff_coeffs({diff_k, diff_na, diff_ca})

CE:set_temperature_celsius(temp)



-- leakage
leak = ChannelLeak("v", allSubsetString)
leak:set_cond("g_l_fct")
leak:set_rev_pot(e_l)
CE:add(leak)

-- K A-type channels
k_a = KA_Golding01("v", allSubsetString)
k_a:set_conductance("g_k_a_fct")
k_a:set_proximality_fct("k_a_is_proximal")
CE:add(k_a)

k_dr = KDR_Golding01("v", allSubsetString)
k_dr:set_conductance(g_k_dr)
CE:add(k_dr)

n_ax = Nax_Golding01("v", allSubsetString)
n_ax:set_conductance("g_na_fct")
CE:add(n_ax)


-- electrode stimulation
-- current, x, y, z, begin, duration
-- the (x,y,z) coords need to specify an edge center!
CE:set_influx(5e-9, -5e-8, 0.0, 0.0, equil_time, 0.001)


domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)

-- some tuning for speed
assTuner = domainDisc:ass_tuner()

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)

-- create operator from discretization
linOp = AssembledLinearOperator(timeDisc)


------------------
-- solver setup	--
------------------
-- linear solver --
linConvCheck = CompositeConvCheck(approxSpace, 20, 2e-26, 1e-08)
linConvCheck:set_component_check("v", 1e-21, 1e-12)
linConvCheck:set_verbose(verbose)

ilu = ILU()
linSolver = LinearSolver()
linSolver:set_preconditioner(ilu)
linSolver:set_convergence_check(linConvCheck)

-------------------
-- time stepping --
-------------------
time = 0.0

-- init solution (equilibrium)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(v_eq)

--[[
-- prepare measurement point and write first measurement
measCoords = {0,0,0}
measPosVector = MakeVec(measCoords[1], measCoords[2], measCoords[3]) -- some arbitrary dendrite vertex pos

measFileVm = outputPath .. "meas/vm_" .. string.format("%.5f", time) .. ".dat"
measOutVm = assert(io.open(measFileVm, "w"))
vm_at_measPt = EvaluateAtClosestVertex(measPosVector, u, "v", "soma", dom:subset_handler())
-- VDCC_BG_VM2UG expects voltages in mV
measOutVm:write(measCoords[1], "\t", measCoords[2], "\t", measCoords[3], "\t", 1e3*vm_at_measPt, "\n")
measOutVm:close()
--]]

-- prepare vtk output
-- NOTE: subdirectory "vtk" needs to exist in output path
if generateVTKoutput then 
	out = VTKOutput()
end

-- store grid function in vector of old solutions
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

min_dt = 1e-10
curr_dt = dt
dtred = 2

lv = 0
cb_counter = {}
cb_counter[lv] = 0
while simEndTime-time > 0.001*curr_dt do
		-- setup time disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, curr_dt)
	
	-- reduce time step if cfl < curr_dt
	-- (this needs to be done AFTER prepare_step as channels are updated there)
	dtChanged = false
	cfl = CE:estimate_cfl_cond(solTimeSeries:latest())
	if cfl < min_dt then
		print("Required time step size is lower than admissible. Aborting.")
		break
	end
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
	
	-- apply linear solver
	ilu:set_disable_preprocessing(matrixIsConst)
	if ApplyLinearSolver(linOp, u, b, linSolver) == false then
		print("Could not apply linear solver.");
		if generateVTKoutput then 
			out:write_time_pvd(outputPath.."vtk/solution", u) 
		end
		exit()
	end
	
	-- update to new time
	time = solTimeSeries:time(0) + curr_dt
	
	--[[
	-- log vm and calcium at soma
	if math.abs(time/dt - math.floor(time/dt+0.5)) < 1e-5 then
		measFileVm = outputPath .. "meas/vm_" .. string.format("%.5f", time) .. ".dat"
		measOutVm = assert(io.open(measFileVm, "w"))
		vm_at_measPt = EvaluateAtClosestVertex(measPosVector, u, "v", "soma", dom:subset_handler())
		-- VDCC_BG_VM2UG expects voltages in mV
		measOutVm:write(measCoords[1], "\t", measCoords[2], "\t", measCoords[3], "\t", 1e3*vm_at_measPt, "\n")
		measOutVm:close()
	end
	--]]
	
	-- vtk output
	if generateVTKoutput then 
		if math.abs((time-equil_time)/pstep - math.floor((time-equil_time)/pstep+0.5)) < 1e-5 and time - equil_time > -pstep/2 then
			out:print(outputPath.."vtk/solution", u, math.floor((time-equil_time)/pstep+0.5), time)
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
	out:write_time_pvd(outputPath.."vtk/solution", u) 
end

	
