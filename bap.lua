--------------------------------------------------------------
-- This script solves the cable equation with HH channels, 	--
-- activating synapses and transmission synapses.			--
--------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")


-- dimension
dim = 3

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"cable_neuron"})


--------------------------------------------------------------------------------
-- Settings
--------------------------------------------------------------------------------
cell = util.GetParam("-cellName", "12-L3pyr")
if not cell == "12-L3pyr" or not cell == "31o_pyr" then
	exit("Cell not specified correctly. Type '12-L3pyr' or '31o_pyr'.")
end

if cell == "12-L3pyr" then
	gridName = "grids/13-L3pyr-77.CNG.ugx"
	gridSyn  = "grids/13-L3pyr-77.CNG_syn.ugx"
	gridDeg  = "grids/13-L3pyr-77.CNG_syn_deg.ugx"
	distro   = {0.0, 0.0, 0.5, 0.5}
	neededSubsets = {"soma", "axon", "axon_myel", "basal", "apical_rad", "apical_lm", "axon_hillock", "ranvier_nodes"}
	dendSubsets = "dendrite, apical_dendrite"
else
	gridName = "grids/31o_pyramidal19aFI.CNG.ugx"
	gridSyn  = "grids/31o_pyramidal19aFI.CNG_syn.ugx"
	gridDeg  = "grids/31o_pyramidal19aFI.CNG_syn_deg.ugx"
	distro   = {0.0, 1.0, 0.0}
	neededSubsets = {"soma", "dendrite", "axon"}
	dendSubsets = "dendrite"
end

-- parameters steering simulation
numRefs		= util.GetParamNumber("-numRefs",		0)
dt			= util.GetParamNumber("-dt",			1e-5) -- in s
endTime		= util.GetParamNumber("-endTime",	  	1.0)  -- in s
nSteps 		= util.GetParamNumber("-nSteps",		endTime/dt)
pstep		= util.GetParamNumber("-pstep",			dt,		"plotting interval")

-- synapse activity parameters
avg_start = util.GetParamNumber("-avgStart"	,  0.03)
avg_dur = util.GetParamNumber(	"-avgDur"	,  2.4e-3)
dev_start = util.GetParamNumber("-devStart"	,  0.015)
dev_dur = util.GetParamNumber(	"-devDur"	,  0.0)
num_synapses = util.GetParamNumber("-nSyn", 140)

-- with simulation of single ion concentrations?
withIons = util.HasParamOption("-ions")

-- specify "-verbose" to output linear solver convergence
verbose	= util.HasParamOption("-verbose")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")

-- profiling?
doProfiling = util.HasParamOption("-profile")
SetOutputProfileStats(doProfiling)

-- file handling
filename = util.GetParam("-outName", "sol_new_clearance_1e-3")
filename = filename .. "/"

--------------------------------------------------------------------------------
-- Synapse distributions via plugin by Lukas Reinhardt
--------------------------------------------------------------------------------
--[[
synDistr = SynapseDistributor(gridName)
synDistr:clear() -- clear any synapses from grid
synDistr:place_synapses(distro, num_synapses, "AlphaPostSynapse")
export_succes = synDistr:export_grid(gridSyn)
print("SynapseDistributor grid export successful: " .. tostring(export_succes))
--]]

gridName = gridSyn

--------------------------------------------------------------------------------
-- Synapse degeneration
--------------------------------------------------------------------------------
---[[
deg_factor = util.GetParamNumber("-degFac", 0.5)
deg_factor = deg_factor + 0.5/num_synapses -- rounding instead of floor-ing

synDistr = SynapseDistributor(gridName)
--synDistr:print_status()
if cell == "12-L3pyr" then
	synDistr:degenerate_uniform(deg_factor, 2) -- first factor means: newNumber = (1-factor)*oldNumber
	synDistr:degenerate_uniform(deg_factor, 3) -- second param is the subset index
else
	synDistr:degenerate_uniform(deg_factor, 1)
end
synDistr:print_status()
synDistr:export_grid(gridDeg)

gridName = gridDeg
--]]

--------------------------------------------------------------
-- File i/o setup for sample calcium concentration measurement
-------------------------------------------------------------- 
measFileVm = filename.."meas/measVm.txt"
measFileCa = filename.."meas/measCa.txt"

if ProcRank() == 0 then
	measOutVm = assert(io.open(measFileVm, "a"))
	measOutCa = assert(io.open(measFileCa, "a"))
end

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


model_cell = 0 --?
if model_cell == 0 then
	na_slope = 0.0 
	g_na_bar = 0.042 --?units?
	g_na_ranvier = 1190.47619 * g_na_bar
else
	na_slope = 0.005
	g_na_bar = 0.03
	g_na_ranvier = 1666.67 * g_na_bar
end

if model_cell > 1 then
	ka_factor=1.3
	g_kap_bar = ka_factor*0.1
	g_kad_bar = g_kap_bar
else then 
	g_kap_bar = 0.1
	g_kad_bar = g_kap_bar
end

spinelimit=100

g_na_basal = g_na_bar
nalimit=275

dslope=0.01
dlimit=300
dprox=100

isegfactor=100
isegfrac=0.8

g_kdr=0.04

g_l_ax = 200.0
g_l_so = 1.0
g_l_de = 1.0

g_l_myelin = 2.5e-1 -- <-- new
g_l_ranvier = 2.0e2 -- <-- new

-- specific capacitance (in units of F/m^2)
spec_cap = 7.5e-3		 -- <-- new
spec_cap_myelin = 7.5e-4 -- <-- new


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
v_eq = -0.07	 -- <-- new

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.0e-9
diff_na	= 1.0e-9
diff_ca	= 2.2e-10

-- temperature in units of deg Celsius
temp = 35.0 	 -- <-- new


--------------------------------------------------------------------------------
-- Create, Load, Refine Domain
--------------------------------------------------------------------------------
dom = util.CreateDomain(gridName, numRefs, neededSubsets)

-- check domain is acyclic
isAcyclic = is_acyclic(dom)
if not isAcyclic then
	print("Domain is not acyclic!")
	exit()
end

--------------------------------------------------------------------------------
-- create Approximation Space
--------------------------------------------------------------------------------
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
if withIons then
	approxSpace:add_fct("k", "Lagrange", 1)
	approxSpace:add_fct("na", "Lagrange", 1)
	approxSpace:add_fct("ca", "Lagrange", 1)
end

approxSpace:init_levels();
approxSpace:init_surfaces();
approxSpace:init_top_surface();
approxSpace:print_layout_statistic()
approxSpace:print_statistic()
OrderCuthillMcKee(approxSpace, true);


-- cable equation
CE = CableEquation("soma, axon, " .. dendSubsets, withIons)

CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)

function blubb(x, y, z, si)
	if (math.sqrt(x*x + y*y + z*z)) < thresh
	then return 1.0
	end
	return 0.0
end
CE:set_spec_res("blubb")


CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)
CE:set_rev_pot_ca(e_ca)

CE:set_k_out(k_out)
CE:set_na_out(na_out)
CE:set_ca_out(ca_out)

CE:set_diff_coeffs({diff_k, diff_na, diff_ca})

CE:set_temperature_celsius(temp)


-- Hodgkin and Huxley channels
if withIons == true then
	HH = ChannelHHNernst("v", "axon, soma, " .. dendSubsets)
else
	HH = ChannelHH("v", "axon, soma, " .. dendSubsets)
end
HH:set_conductances(g_k_ax, g_na_ax, "axon")
HH:set_conductances(g_k_so, g_na_so, "soma")
HH:set_conductances(g_k_de, g_na_de, dendSubsets)

CE:add(HH)


-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leak = ChannelLeak("v", "axon, soma, " .. dendSubsets)
leak:set_cond(g_l_ax*tmp_fct, "axon")
leak:set_rev_pot(-0.066210342630746467, "axon") -- golding sets: -66mV in Cs+  -- <-- new
leak:set_cond(g_l_so*tmp_fct, "soma")
leak:set_rev_pot(-0.022074360525636, "soma")
leak:set_cond(g_l_de*tmp_fct, dendSubsets)
leak:set_rev_pot(-0.056314322586687, dendSubsets)

CE:add(leak)


-- Calcium dynamics
if withIons then
	vdcc = VDCC_BG_cable("ca", "soma, " .. dendSubsets)
	ncx = NCX_cable("ca", "soma, " .. dendSubsets)
	pmca = PMCA_cable("ca", "soma, " .. dendSubsets)
	caLeak = IonLeakage("ca", "soma, " .. dendSubsets)
	leakCaConst = -3.4836065573770491e-9 +	-- single pump PMCA flux density (mol/s/m^2)
				  -1.0135135135135137e-9 +	-- single pump NCX flux (mol/s/m^2)
				  3.3017662162505882e-11
	caLeak:set_perm(leakCaConst, ca_in, ca_out, v_eq, 2)
	
	CE:add(ncx)
	CE:add(pmca)
	CE:add(vdcc)
	CE:add(caLeak)
end


-- channels - by nmodl converter

-- for distance dependent values
function distance_dep_res(x, y, z, si)
	if (math.sqrt(x*x + y*y + z*z)) >= spinelimit
	then return spinefactor/spec_res
	else return spec_res
	end
end

function distance_dep_cap(x, y, z, si)
	if (math.sqrt(x*x + y*y + z*z)) >= spinelimit
	then return spinefactor*spec_cap
	else return spec_cap
	end
end

function distance_dep_nax(x, y, z, si)
	if (math.sqrt(x*x + y*y + z*z)) >= nalimit
	then return (1+limit*na_slope)*g_na_bar
	else return g_na_bar
	end
end

function iseg_nax(x, y, z, si)
	if (math.sqrt(x*x + y*y + z*z)) >= isegfrac
	then return isegfactor*g_na_bar
	else return g_na_bar
	end
end



function distance_dep_d(x, y, z, si)
	if (math.sqrt(x*x + y*y + z*z)) >= dlimit
	then xdist = dlimit
	end
	if xdist >= dprox then return g_kad_bar*(1+xdist*dslope)
	else then return g_kap_bar*(1+xdist*dslope)
	end
end


--passive
CE:set_spec_res("distance_dep_res", "apical_rad")
CE:set_spec_res("distance_dep_res", "apical_lm")
CE:set_spec_cap("distance_dep_cap", "apical_rad")
CE:set_spec_cap("distance_dep_cap", "apical_lm")


--active
--na
nax:set_conductances(g_na_bar, "axon")
nax:set_conductances(g_na_bar, "axon_hillock")
nax:set_conductances(g_na_bar, "soma")
nax:set_conductances(g_na_bar, "basal")
nax:set_conductances(g_na_ran, "ranvier_nodes")
nax:set_conductances("distance_dep_nax", "apical_rad")
nax:set_conductances("distance_dep_nax", "apical_lm")
--nax:set_conductances(g_na_bar, "iseg") intial segment??


--ka
kap:set_conductances(g_kap_bar, "axon")
kap:set_conductances(g_kap_bar, "axon_hillock")
kap:set_conductances(g_kap_bar, "soma")
kap:set_conductances(g_kap_bar, "basal")
kap:set_conductances(g_kap_bar*0.2, "ranvier_nodes")
kap:set_conductances("distance_dep_d", "apical_rad")
kap:set_conductances("distance_dep_d", "apical_lm")
--kap:set_conductances(g_kap_bar, "iseg")

--kdr
kdr:set_conductances(g_kdr, "axon")
kdr:set_conductances(g_kdr, "axon_hillock")
kdr:set_conductances(g_kdr, "soma")
kdr:set_conductances(g_kdr, "basal")
kdr:set_conductances(g_kdr, "ranvier_nodes")
kdr:set_conductances(g_kdr, "apical_rad")
kdr:set_conductances(g_kdr, "apical_lm")
--kdr:set_conductances(g_kdr, "iseg")

--[[
-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing_alpha(
	avg_start,	 -- average onset of synaptical activity in [s]
	avg_dur/6.0, -- average tau of activity function in [s]
	dev_start,   -- deviation of onset in [s]
	dev_dur/6.0, -- deviation of tau in [s]
	1.2e-9)		 -- peak conductivity in [S]
CE:set_synapse_handler(syn_handler)
--]]

--[[
-- electrode stimulation
-- 5nA seem to enervate the pyramidal cell with uniform diameters of 1um
-- (coords for 13-L3pyr-77.CNG.ugx, current given in C/ms)
CE:set_influx(5e-9, 6.54e-05, 2.665e-05, 3.985e-05, 0.0, 0.04)			-- near soma
CE:set_influx(5e-9, 3.955e-06, 1.095e-06, -3.365e-06, 0.001, 0.0025)		-- 1st edge soma to dend
CE:set_influx(0.3e-9, 3.955e-06, 1.095e-06, -3.365e-06, 0.0, 0.03)		-- 1st 1st edge soma to dend
CE:set_influx(0.095e-9, 0.0, 0.0, 0.0, 0.1, 0.1)							-- soma center vertex
CE:set_influx(0.2e-9, 0.0, 0.0, 0.0, 0.005, 0.0005)						-- soma center vertex
CE:set_influx(10.0e-9, 0.000139, 0.00020809, -2.037e-05, 0.005, 0.005)	-- distal apical dendrite vertex v1
CE:set_influx(10.0e-9, -3.96e-06, 0.0002173, -5.431e-05, 0.005, 0.005)	-- distal apical dendrite vertex v2
--]]


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
linConvCheck = CompositeConvCheck(approxSpace, 20, 2e-26, 1e-08)
linConvCheck:set_component_check("v", 1e-21, 1e-12)
linConvCheck:set_verbose(verbose)

ilu = ILU()
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(linConvCheck)
--cgSolver:set_debug(dbgWriter)

----------------------
-- time stepping	--
----------------------
time = 0.0

-- init solution
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
Interpolate(v_eq, u, "v")
if withIons then
	Interpolate(k_in, u, "k");
	Interpolate(na_in, u, "na");
	Interpolate(ca_in, u, "ca")
end


-- write start solution
if generateVTKoutput then 
	out = VTKOutput()
	out:print(filename.."vtk/solution", u, 0, time)
end

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
	if ApplyLinearSolver(linOp, u, b, cgSolver) == false then
		print("Could not apply linear solver.")
		exit()
	end
	
	-- log time and vm in Soma
	if ProcRank() == 0 then
		if cell == "12-L3pyr" then
			vm_soma  = EvaluateAtClosestVertex(MakeVec(0.0, 0.0, 0.0), 						u, "v", "soma", 		dom:subset_handler())
			vm_axon  = EvaluateAtClosestVertex(MakeVec(-3.828e-05, -0.00013166, -2.34e-05), u, "v", "axon", 		dom:subset_handler())
			vm_dend  = EvaluateAtClosestVertex(MakeVec(8.304e-05, -1.982e-05, -8.4e-06), 	u, "v", "dendrite", 		dom:subset_handler())
			vm_aDend = EvaluateAtClosestVertex(MakeVec(-3.84e-06, 0.00018561, -3.947e-05), 	u, "v", "apical_dendrite", 	dom:subset_handler())
			measOutVm:write(time, "\t", vm_soma, "\t", vm_axon, "\t", vm_dend, "\t", vm_aDend, "\n")
			if (withIons) then
				ca_soma  = EvaluateAtClosestVertex(MakeVec(0.0, 0.0, 0.0), 						u, "ca", "soma", 		dom:subset_handler())
				ca_axon  = EvaluateAtClosestVertex(MakeVec(-3.828e-05, -0.00013166, -2.34e-05), u, "ca", "axon", 		dom:subset_handler())
				ca_dend  = EvaluateAtClosestVertex(MakeVec(8.304e-05, -1.982e-05, -8.4e-06), 	u, "ca", "dendrite", 		dom:subset_handler())
				ca_aDend = EvaluateAtClosestVertex(MakeVec(-3.84e-06, 0.00018561, -3.947e-05), 	u, "ca", "apical_dendrite", 	dom:subset_handler())
				measOutCa:write(time, "\t", ca_soma, "\t", ca_axon, "\t", ca_dend, "\t", ca_aDend, "\n")
			end
		else	
			vm_soma  = EvaluateAtClosestVertex(MakeVec(6.9e-07, 3.74e-06, -2.86e-06), 		u, "v", "soma", 		dom:subset_handler())
			vm_axon  = EvaluateAtClosestVertex(MakeVec(-4.05e-06, 6.736e-05, -1.341e-05), 	u, "v", "axon", 		dom:subset_handler())
			vm_dend  = EvaluateAtClosestVertex(MakeVec(-4.631e-05, -0.0001252, 4.62e-06), 	u, "v", "dendrite", 	dom:subset_handler())
			measOutVm:write(time, "\t", vm_soma, "\t", vm_axon, "\t", vm_dend, "\t", -65, "\n")
			if (withIons) then
				ca_soma  = EvaluateAtClosestVertex(MakeVec(6.9e-07, 3.74e-06, -2.86e-06), 		u, "ca", "soma", 		dom:subset_handler())
				ca_axon  = EvaluateAtClosestVertex(MakeVec(-4.05e-06, 6.736e-05, -1.341e-05), 	u, "ca", "axon", 		dom:subset_handler())
				ca_dend  = EvaluateAtClosestVertex(MakeVec(-4.631e-05, -0.0001252, 4.62e-06), 	u, "ca", "dendrite", 	dom:subset_handler())
				measOutCa:write(time, "\t", ca_soma, "\t", ca_axon, "\t", ca_dend, "\t", -65, "\n")
			end
		end
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

-- close measure file
if ProcRank() == 0 then
	measOutVm:close()
	measOutCa:close()
end

if doProfiling then
	WriteProfileData(filename .."pd.pdxml")
end



	
