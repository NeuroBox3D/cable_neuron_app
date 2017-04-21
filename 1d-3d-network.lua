-------------------------------------------------------------------------
-- This script sets up a 1d/3d hybrid simulation with a 1d network and --
-- the 3d representation of one of the network cells.                  --
-- On the 1d domain, it solves the cable equation with HH channels,    --
-- activating specifically set synapses. On the 3d domain, it solves   --
-- a calcium problem (diffusion and buffering) with channels and pumps --
-- in the plasma membrane, where VDCCs are activated according to the  --
-- potential mapped from the 1d domain.                                --
--                                                                     --
-- author: mbreit                                                      --
-- date:   23-03-2017                                                  --
-------------------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"cable_neuron", "neuro_collection"})

dim = 3
InitUG(dim, AlgebraType("CPU", 1));


-- choice of grid and refinement level
gridName1d = util.GetParam("-grid1d", "../apps/cable_neuron_app/grids/NMO_05371_1d.ugx")
gridName3d = util.GetParam("-grid3d", "../apps/cable_neuron_app/grids/NMO_05371_3d_6.ugx")
numRefs = util.GetParamNumber("-numRefs", 0)

-- parameters for instationary simulation
dt1d = util.GetParamNumber("-dt1d", 1e-5) -- in s
dt3d = util.GetParamNumber("-dt3d", 1e-2) -- in s
dt3dStart = util.GetParamNumber("-dt3dstart", dt3d)
endTime = util.GetParamNumber("-endTime", 1.0)  -- in s

-- with simulation of single ion concentrations?
withIons = util.HasParamOption("-ions")

-- choice of solver setup
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0;
validSolverIDs["GS"] = 0;
validSolverIDs["ILU"] = 0;
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- specify "-verbose" to output linear solver convergence
verbose1d = util.HasParamOption("-verbose1d")
verbose3d = util.HasParamOption("-verbose3d")

-- vtk output?
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt3d, "plotting interval")

-- file handling
filename = util.GetParam("-outName", "hybrid_test")
filename = filename.."/"


-- choose length of time step at the beginning
-- if not timeStepStart = 2^(-n)*timeStep, take nearest lower number of that form
function log2(x)
	return math.log(x) / math.log(2)
end
startLv =  math.ceil(log2(dt3d / dt3dStart))
dt3dStartNew = dt3d / math.pow(2, startLv)
if (math.abs(dt3dStartNew - dt3dStart) / dt3dStart > 1e-5) then 
	print("dt3dStart argument ("..dt3dStart..") was not admissible; taking "..dt3dStartNew.." instead.")
end
dt3dStart = dt3dStartNew


print("Chosen parameters:")
print("    grid       = " .. gridName1d)
print("    numRefs    = " .. numRefs)
print("    dt1d       = " .. dt1d)
print("    dt3d       = " .. dt3d)
print("    dt3dStart  = " .. dt3dStart)
print("    endTime    = " .. endTime)
print("    pstep      = " .. pstep)
print("    ions       = " .. tostring(withIons))
print("    solver     = " .. solverID)
print("    verbose1d  = " .. tostring(verbose1d))
print("    verbose3d  = " .. tostring(verbose3d))
print("    vtk        = " .. tostring(generateVTKoutput))
print("    outname    = " .. filename)


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
v_eq = -0.065

-- diffusion coefficients (in units of m^2/s)
diff_k 	= 1.0e-9
diff_na	= 1.0e-9
diff_ca	= 2.2e-10

-- temperature in units of deg Celsius
temp = 37.0


------------------------------------
-- create 1d domain and approx space --
------------------------------------
neededSubsets1d = {}
dom1d = util.CreateDomain(gridName1d, 0, neededSubsets1d)

approxSpace1d = ApproximationSpace(dom1d)
approxSpace1d:add_fct("v", "Lagrange", 1)
if withIons == true then
	approxSpace1d:add_fct("k", "Lagrange", 1)
	approxSpace1d:add_fct("na", "Lagrange", 1)
	approxSpace1d:add_fct("ca", "Lagrange", 1)
end

approxSpace1d:init_levels();
approxSpace1d:init_surfaces();
approxSpace1d:init_top_surface();
approxSpace1d:print_layout_statistic()
approxSpace1d:print_statistic()

------------------------------
-- create 1d discretization --
------------------------------
ss_axon = "AXON__L4_STELLATE, AXON__L23_PYRAMIDAL, AXON__L5A_PYRAMIDAL, AXON__L5B_PYRAMIDAL"
ss_dend = "DEND__L4_STELLATE, DEND__L23_PYRAMIDAL, DEND__L5A_PYRAMIDAL, DEND__L5B_PYRAMIDAL"
ss_soma = "SOMA__L4_STELLATE, SOMA__L23_PYRAMIDAL, SOMA__L5A_PYRAMIDAL, SOMA__L5B_PYRAMIDAL"
	
-- cable equation
CE = CableEquation(ss_axon..", "..ss_dend..", "..ss_soma, withIons)
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
if withIons == true then
	HH = ChannelHHNernst("v", ss_axon..", "..ss_dend..", "..ss_soma)
else
	HH = ChannelHH("v", ss_axon..", "..ss_dend..", "..ss_soma)
end
HH:set_conductances(g_k_ax, g_na_ax, ss_axon)
HH:set_conductances(g_k_so, g_na_so, ss_soma)
HH:set_conductances(g_k_de, g_na_de, ss_dend)

CE:add(HH)


-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

-- calculate leakage reversal potential for equilibrium at given potential
function leakageReversal(gNa, gK, gL, ENa, EK, Veq, m, h, n)
	return Veq + (gNa*m*m*m*h*(Veq - ENa) + gK*n*n*n*n*(Veq - EK)) / gL
end

m_inf = 0.0529325
h_inf = 0.596121
n_inf = 0.317677
leakRev_ax = leakageReversal(g_na_ax, g_k_ax, g_l_ax, e_na, e_k, v_eq, m_inf, h_inf, n_inf)
leakRev_de = leakageReversal(g_na_de, g_k_de, g_l_de, e_na, e_k, v_eq, m_inf, h_inf, n_inf)
leakRev_so = leakageReversal(g_na_so, g_k_so, g_l_so, e_na, e_k, v_eq, m_inf, h_inf, n_inf)

leak = ChannelLeak("v", ss_axon..", "..ss_dend..", "..ss_soma)
leak:set_cond(g_l_ax*tmp_fct, ss_axon)
leak:set_rev_pot(leakRev_ax, ss_axon)
leak:set_cond(g_l_so*tmp_fct, ss_soma)
leak:set_rev_pot(leakRev_so, ss_soma)
leak:set_cond(g_l_de*tmp_fct, ss_dend)
leak:set_rev_pot(leakRev_de, ss_dend)

CE:add(leak)


-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
CE:set_synapse_handler(syn_handler)


-- domain discretization
domDisc1d = DomainDiscretization(approxSpace1d)
domDisc1d:add(CE)

-------------------------
-- 1d domain distribution --
-------------------------
-- Domain distribution needs to be performed AFTER addition
-- of the synapse handler to the CE object and addition of the
-- CE object to the domain disc (i.e.: when the synapse handler
-- has got access to the grid).
-- The reason is that the synapse handler needs access to the grid
-- to correctly distribute the synapse* attachments.

balancer.partitioner = "parmetis"
balancer.firstDistLvl = -1
balancer.redistSteps = 0
balancer.imbalanceFactor = 1.05
balancer.staticProcHierarchy = true

balancer.ParseParameters()
balancer.PrintParameters()

loadBalancer1d = balancer.CreateLoadBalancer(dom1d)
if loadBalancer1d ~= nil then
	loadBalancer1d:enable_vertical_interface_creation(false)
	balancer.Rebalance(dom1d, loadBalancer1d)
	
	edgeCut = balancer.defaultPartitioner:edge_cut_on_lvl(0)
	if edgeCut ~= 0 then
		print("Network is not partitioned into whole cells.")
		print("Edge cut on base level: " .. edgeCut)
		exit(1)
	end
	loadBalancer1d:estimate_distribution_quality()
	loadBalancer1d:print_quality_records()
end

--SaveParallelGridLayout(dom1d:grid(), filename .. "grid/parallel_grid1d_layout_p"..ProcRank()..".ugx", 0)

-- ordering; needs to be done after distribution!
order_cuthillmckee(approxSpace1d)


-- find neuron IDs for 3d simulation
nid = innermost_neuron_id_in_subset("SOMA__L5A_PYRAMIDAL", dom1d:subset_handler())
print("Innermost neuron ID: "..nid)



----------------------------------
-- constants for the 3d problem --
----------------------------------
-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = 4*2.5e-6

-- diffusion coefficients
D_cac = 220.0
D_clb = 20.0

-- calbindin binding rates
k_bind_clb = 	27.0e06
k_unbind_clb = 	19

-- initial concentrations
ca_cyt_init = 5.0e-08
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)

pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 1.0

leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				+ vdccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end

----------------------------------
-- setup 3d approximation space --
----------------------------------
-- create, load, refine and distribute domain
reqSubsets = {"cyt", "pm"}
dom3d = util.CreateDomain(gridName3d, 0)
balancer.partitioner = "parmetis"

balancer.staticProcHierarchy = true
balancer.firstDistLvl = -1
balancer.redistSteps = 0

balancer.ParseParameters()
balancer.PrintParameters()

-- in parallel environments: use a load balancer to distribute the grid
-- actual refinement and load balancing after setup of disc.
loadBalancer3d = balancer.CreateLoadBalancer(dom3d)

-- refining and distributing
-- manual refinement (need to update interface node location in each step)
if loadBalancer3d ~= nil then
	loadBalancer3d:enable_vertical_interface_creation(false)
	balancer.Rebalance(dom3d, loadBalancer3d)
end

if numRefs > 0 then	
	refiner = GlobalDomainRefiner(dom3d)	
	for i = 1, numRefs do
		refiner:refine()
	end
end

if loadBalancer3d ~= nil then
	print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
	loadBalancer3d:estimate_distribution_quality()
	loadBalancer3d:print_quality_records()
end
print(dom3d:domain_info():to_string())


--[[
--print("Saving domain grid and hierarchy.")
--SaveDomain(dom3d, "refined_grid_p" .. ProcRank() .. ".ugx")
--SaveGridHierarchyTransformed(dom3d:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
--print("Saving parallel grid layout")
SaveParallelGridLayout(dom3d:grid(), filename.."parallel_grid3d_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
approxSpace3d = ApproximationSpace(dom3d)

cytVol = "cyt"
plMem = "pm"
plMem_vec = {"pm"}

outerDomain = cytVol .. ", " .. plMem

approxSpace3d:add_fct("ca_cyt", "Lagrange", 1)
approxSpace3d:add_fct("clb", "Lagrange", 1)

approxSpace3d:init_levels();
approxSpace3d:init_surfaces();
approxSpace3d:init_top_surface();
approxSpace3d:print_layout_statistic()
approxSpace3d:print_statistic()

OrderCuthillMcKee(approxSpace3d, true)

--------------------------
-- setup discretization --
--------------------------
-- diffusion --
diffCa = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
diffCa:set_diffusion(D_cac)

diffClb = ConvectionDiffusion("clb", cytVol, "fv1")
diffClb:set_diffusion(D_clb)

-- buffering --
elemDiscBuffering = BufferFV1(cytVol) -- where buffering occurs
elemDiscBuffering:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant

-- plasma membrane transport systems
pmca = PMCA({"ca_cyt", ""})
pmca:set_constant(1, 1.0)
pmca:set_scale_inputs({1e3,1.0})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ncx = NCX({"ca_cyt", ""})
ncx:set_constant(1, 1.0)
ncx:set_scale_inputs({1e3,1.0})
ncx:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakPM = Leak({"", "ca_cyt"})
leakPM:set_constant(0, 1.0)
leakPM:set_scale_inputs({1.0,1e3})
leakPM:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)

vdcc = VDCC_BG_CN({"ca_cyt", ""}, plMem_vec, approxSpace1d, approxSpace3d, "v")
vdcc:set_domain_disc_1d(domDisc1d)
vdcc:set_cable_disc(CE)
vdcc:set_3d_neuron_ids({nid})
vdcc:set_coordinate_scale_factor_3d_to_1d(1e-6)
if withIons then
	vdcc:set_initial_values({v_eq, k_in, na_in, ca_in})
else
	vdcc:set_initial_values({v_eq})
end
vdcc:set_time_steps_for_simulation_and_potential_update(dt1d, dt1d)
vdcc:set_solver_output_verbose(verbose1d)
vdcc:set_vtk_output(filename.."vtk/solution1d", pstep)
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L() -- default, but to be sure


discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)

discNCX = MembraneTransportFV1(plMem, ncx)
discNCX:set_density_function(ncxDensity)

discLeak = MembraneTransportFV1(plMem, leakPM)
discLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))

discVDCC = MembraneTransportFV1(plMem, vdcc)
discVDCC:set_density_function(vdccDensity)

-- synaptic activity --
synapseInflux = HybridSynapseCurrentAssembler(approxSpace3d, approxSpace1d, syn_handler, {"pm"}, "ca_cyt")
synapseInflux:set_current_percentage(0.01)
synapseInflux:set_3d_neuron_ids({nid})
synapseInflux:set_scaling_factors(1e-15, 1e-6, 1.0)
synapseInflux:set_valency(2)


-- domain discretization --
domDisc3d = DomainDiscretization(approxSpace3d)

domDisc3d:add(diffCa)
domDisc3d:add(diffClb)

domDisc3d:add(elemDiscBuffering)

domDisc3d:add(discPMCA)
domDisc3d:add(discNCX)
domDisc3d:add(discLeak)
domDisc3d:add(discVDCC)

domDisc3d:add(synapseInflux)


-- setup time discretization --
timeDisc = ThetaTimeStep(domDisc3d)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()


------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace3d)
dbgWriter:set_base_dir(filename)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-8)
convCheck:set_verbose(verbose3d)

if (solverID == "ILU") then
    bcgs_steps = 10000
    bcgs_precond = ILU()
elseif (solverID == "GS") then
    bcgs_steps = 10000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace3d)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(0)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	gmg:set_base_solver(SuperLU())
	
	ilu_gmg = ILU()
	ilu_gmg:set_sort(true)		-- <-- SUPER-important!
	gmg:set_smoother(ilu_gmg)
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_rap(true) -- causes error in base solver!!
	--gmg:set_debug(GridFunctionDebugWriter(approxSpace))
	
    bcgs_steps = 1000
	bcgs_precond = gmg
end

convCheck:set_maximum_steps(bcgs_steps)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace3d, 10, 1e-19, 1e-12)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
--newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)


-------------
-- solving --
-------------
-- get grid function
u = GridFunction(approxSpace3d)

-- set initial value
InterpolateInner(ca_cyt_init, u, "ca_cyt", 0.0)
InterpolateInner(clb_init, u, "clb", 0.0)

-- timestep in seconds
dt = dt3dStart
time = 0.0
step = 0

-- initial vtk output
if (generateVTKoutput) then
	out = VTKOutput()
	out:print(filename .. "vtk/solution3d", u, step, time)
end

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

min_dt = dt3d / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = 0.01
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- apply newton solver
	if newtonSolver:apply(u) == false
	then
		-- in case of failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
		
		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Trying with half the time step...")
			cb_counter[lv] = 0
		end
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- update check-back counter and, if applicable, reset dt
		cb_counter[lv] = cb_counter[lv] + 1
		while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 and (time >= levelUpDelay or lv > startLv) do
			print ("Doubling time due to continuing convergence; now: " .. 2*dt)
			dt = 2*dt;
			lv = lv - 1
			cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
			cb_counter[lv+1] = 0
		end
		
		-- plot solution every pstep seconds
		if (generateVTKoutput) then
			if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then
				out:print(filename .. "vtk/solution3d", u, math.floor(time/pstep+0.5), time)
			end
		end
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end

end

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(filename .. "vtk/solution3d", u) end

