--------------------------------------------------------------------------------
-- This script solves the cable equation on a 1d version of Viet's            --
-- reconstructed spine (with extensions).                                     --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2019-01-17                                                         --
--------------------------------------------------------------------------------
ug_load_script("ug_util.lua")

-- init UG
InitUG(3, AlgebraType("CPU", 1));
AssertPluginsLoaded({"cable_neuron"})


--------------
-- settings --
--------------
-- grid
gridName = util.GetParam("-grid", "nernst_planck_app/grid/spines/vietSpine1d.ugx")

-- refinement
numRefs = util.GetParamNumber("-numRefs", 0)

-- time stepping
dt = util.GetParamNumber("-dt", 1e-5)  -- in s
endTime = util.GetParamNumber("-endTime", 0.05)  -- in s
nSteps = util.GetParamNumber("-nSteps", endTime/dt)

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt)

-- file handling
outDir = util.GetParam("-outName", "viet_1d")
outDir = outDir .. "/"


--------------------------
-- biological settings	--
--------------------------
-- membrane conductances (in units of S/m^2)
g_l_k = 0.5
g_l_na = 0.5

-- resistivity (in units of (V s m) / C)
spec_res = 0.407224494  -- deduced from PNP model

-- membrane specific capacitance (in units of C / (V m^2))
spec_cap = 3.585763867e-3  -- deduced from PNP model

-- reversal potentials (in units of V)
e_k  = -0.0940
e_na = 0.0640
e_ca = 0.14

-- equilibrium concentrations (in units of mM)
k_out  = 4.0
na_out = 145.0
ca_out = 2.0
cl_out = 123.0

k_in = 155.0
na_in = 12.0
ca_in = 5e-5
cl_in = 5.0

-- equilibrium potential (in units of V)
v_eq = -0.070

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.96e-9
diff_na	= 1.33e-9
diff_ca = 2.2e-10
diff_cl = 2.03e-9

-- temperature in units of K
temp = 298.15


--------------------------------
-- create approximation space --
--------------------------------
dom = util.CreateDomain(gridName, numRefs, {"dend"})

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true)

u = GridFunction(approxSpace)

-- cable equation
CE = CableEquation("dend, spine, syn", false)
CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)
CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)
CE:set_rev_pot_ca(e_ca)
CE:set_k_out(k_out)
CE:set_na_out(na_out)
CE:set_ca_out(ca_out)
CE:set_diff_coeffs({diff_k, diff_na, diff_ca})
CE:set_temperature(temp)

-- leakage
leakK = ChannelLeak("k", "dend")
leakK:set_cond(g_l_k, "dend")
leakK:set_rev_pot(v_eq, "dend")
CE:add(leakK)

leakNa = ChannelLeak("na", "dend")
leakNa:set_cond(g_l_na, "dend")
leakNa:set_rev_pot(v_eq, "dend")
CE:add(leakNa)

-- synapse
nGABAR = 50
local t_on = 0.0
local tau = 0.008  -- GABAR activation decay constant (in s)
local period = 0.005
function channelOpen(t)
	local t_p = t % period
	return math.exp(-t_p / tau)
end

FRT = 96485.0 / (8.31451 * 298.15)  -- in 1 / V
synArea = Integral(1.0, u, "syn")
areaFactor = 9.275e-14 / (math.pi * 0.5*(2.95945e-07 + 1.78e-07) * synArea)

equiv_channel_area = 2.257e-19  -- cross-sectional area of GABAR channel, in m^2,
                                -- deduced from GHK equation using conductance g = 24pS
                                -- and reference concentration c_r = 145mM (inside and outside)
channelAreaPerPSDArea = nGABAR*equiv_channel_area / (math.pi * 0.5*(2.95945e-07 + 1.78e-07) * synArea)

vmAtSyn = v_eq
function currentDensityFunction(x, y, z, t, si)
	local exp = math.exp(-FRT*vmAtSyn)
	local currentCl = 96485.0 * diff_cl / 1e-8 * FRT * vmAtSyn
		* (cl_out - cl_in * exp) / (1.0 - exp)
	return - channelOpen(t) * channelAreaPerPSDArea * currentCl * areaFactor
end

CE:set_influx_function("currentDensityFunction", "syn")

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
convCheck = CompositeConvCheck(approxSpace, 1, 2e-26, 1e-08)
convCheck:set_component_check("v", 1e-21, 1e-12)
convCheck:set_verbose(verbose)

ilu = ILU()
solver = LinearSolver()
solver:set_preconditioner(ilu)
solver:set_convergence_check(convCheck)
--solver:set_debug(dbgWriter)


----------------------
-- time stepping	--
----------------------
time = 0.0

-- init solution
b = GridFunction(approxSpace)
u:set(0.0)
Interpolate(v_eq, u, "v")

synArea = Integral(1.0, u, "syn")
areaFactor = 9.275e-14 / (math.pi * 0.5*(2.95945e-07 + 1.78e-07) * synArea)
currentOutFile = assert(io.open(outDir.."meas/current.dat", "a"))

-- write start solution
if generateVTKoutput then 
	out = VTKOutput()
	out:print(outDir.."vtk/solution", u, 0, time)
end
measFcts = "v"
take_measurement(u, time, "syn", measFcts, outDir.."meas/sol")

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
	
	vmAtSyn = Integral(u, "v", "syn") / synArea
	currentOutFile:write(time, "\t", -currentDensityFunction(0,0,0,time,0) * math.pi * 0.5*(2.95945e-07 + 1.78e-07) * synArea, "\n")
	
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
	if generateVTKoutput then
		if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
			out:print(outDir.."vtk/solution", u, math.floor(time/pstep+0.5), time)
		end
	end
	
	take_measurement(u, time, "syn", measFcts, outDir.."meas/sol")
	
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
	out:write_time_pvd(outDir.."vtk/solution", u) 
end

currentOutFile:close()

	
