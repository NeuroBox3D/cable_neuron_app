-----------------------------------------------
-- This script sets up a network simulation. --
-----------------------------------------------

-- load utility functions
ug_load_script("cable_util.lua")

-- choice of morphology
gridName = "testNetwork.ugx"

-- the following geometric subsets need to be defined in the network file
neededSubsets = {"Axon", "Dendrite", "Soma", "PreSynapseEdges", "PostSynapseEdges"}

-- init UG (for 3-dim coordinate space)
InitUG(3, AlgebraType("CPU", 1));

-----------------------------
-- biophysical model setup --
-----------------------------
-- equilibrium potential (in units of V)
v_eq = -0.065

-- temperature in units of deg Celsius
temp = 37.0

-- cable equation (here: only potential, no ions)
CE = CableEquation("Axon, Dendrite, Soma, PreSynapseEdges, PostSynapseEdges", false)
CE:set_spec_cap(1e-2)		-- specific capacitance (units of F/m^2)
CE:set_spec_res(1.5)		-- specific resistance (units of Ohm m)
CE:set_rev_pot_k(-0.09)			-- reversal potential K (V)
CE:set_rev_pot_na(0.06)			-- reversal potential Na (V)
CE:set_temperature_celsius(temp)

-- Hodgkin and Huxley channels --
-- specific membrane conductances (in units of S/m^2)
g_k_ax = 400.0	-- axon
g_k_so = 200.0	-- soma
g_k_de = 30.0	-- dendrite

g_na_ax = 3.0e4
g_na_so = 1.5e3
g_na_de = 40.0

HH = ChannelHH("", "Axon, Dendrite, Soma, PreSynapseEdges, PostSynapseEdges")
HH:set_conductances(g_k_ax, g_na_ax, "Axon, PreSynapseEdges")
HH:set_conductances(g_k_so, g_na_so, "Soma")
HH:set_conductances(g_k_de, g_na_de, "Dendrite, PostSynapseEdges")
CE:add(HH)

-- leakage (needs to be precisely calibrated to attain the aspired
-- resting potential in each part of the network! -- in our case: -65mV everywhere)
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)
leak = ChannelLeak("v", "Axon, Dendrite, Soma, PreSynapseEdges, PostSynapseEdges")
leak:set_cond(200.0*tmp_fct, "Axon, PreSynapseEdges")
leak:set_rev_pot(-0.066148458, "Axon, PreSynapseEdges")
leak:set_cond(1.0*tmp_fct, "Soma")
leak:set_rev_pot(-0.030654022, "Soma")
leak:set_cond(1.0*tmp_fct, "Dendrite, PostSynapseEdges")
leak:set_rev_pot(-0.057803624, "Dendrite, PostSynapseEdges")
CE:add(leak)

-- synapses
syn_handler = NETISynapseHandler()
syn_handler:set_cable_equation(CE)

-- tell synapse handler where presynapses are located
syn_handler:set_presyn_subset("PreSynapse")

-- define an excitatory input pattern for primary synapses
-- activating the network
syn_handler:set_activation_timing(
	0.0,	-- average start time of synaptical activity in ms
	2.4,	-- average duration of activity in ms (10)
	0.0,	-- deviation of start time in ms
	0.0,	-- deviation of duration in ms
	1.2e-3)	-- peak conductivity in units of uS
CE:set_synapse_handler(syn_handler)

----------------------------------------------------
-- handing the problem over to UG4 and running it --
----------------------------------------------------
-- reads the network file and distributes network to available processors
-- sets up approximation space for unknown membrane potential
-- sets up space and time discretizations for the given model
-- sets up a solver
ugEnv = setup_problem_ug4(CE, gridName, neededSubsets)

-- parameters steering the simulation
-- can be provided on the command line,
-- default value defined here is used otherwise
simParams = {}
simParams.v_init = v_eq
simParams.startTime = 0.0
simParams.dt = util.GetParamNumber("-dt", 2e-5, "integration time step (s)")
simParams.endTime = util.GetParamNumber("-endTime",	0.1,	"end time of simulation (s)")
simParams.generateVTKOutput = util.HasParamOption("-vtk")
simParams.vtkFolder = util.GetParam("-outName", "simulation_output")
simParams.plotStep = util.GetParamNumber("-plotStep",	dt,		"plotting interval")

-- runs the specified simulation
run_simulation(ugEnv, simParams)


