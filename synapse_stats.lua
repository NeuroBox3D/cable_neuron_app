------------------------------------------------------
-- This script is intended for testing purposes.	--
-- It solves the cable equation with HH channels,	--
-- activating synapses and transmission synapses.	--
------------------------------------------------------

---------------------------------------------------
-- THIS SCRIPT IS NOT FUNCTIONAL AT THE MOMENT ! --
-- Repairing this would require re-implementing  --
-- the print_synapse_statistics() method for the --
-- current SynapseHandler.                       --
---------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- choice of grid
gridName = util.GetParam("-grid", "testNetwork.ugx")

-- dimension
dim = 3

-- init UG
InitUG(dim, AlgebraType("CPU", 1));
AssertPluginsLoaded({"cable_neuron"})


-- create, load, refine and distribute domain
neededSubsets = {"Axon", "Dendrite", "Soma"}
dom = util.CreateDomain(gridName, 0, neededSubsets)

balancer.partitioner = "parmetis"
balancer.firstDistLvl = -1
balancer.redistSteps  = 0
balancer.staticProcHierarchy = true
balancer.ParseParameters()

-- in parallel environments: use a load balancer to distribute the grid
loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	loadBalancer:enable_vertical_interface_creation(false)
	balancer.RefineAndRebalanceDomain(dom, 0, loadBalancer)
end
print(dom:domain_info():to_string())

-- create Approximation Space
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:init_levels();
approxSpace:init_surfaces();
approxSpace:init_top_surface();


-- cable equation
CE = CableEquation("Axon, Dendrite, Soma, PreSynapseEdges, PostSynapseEdges", false)

-- synapses
syn_handler = SynapseHandler()
syn_handler:set_presyn_subset("PreSynapse")
syn_handler:set_ce_object(CE)
syn_handler:set_activation_timing(
	5.0,	-- average start time of synaptical activity in ms
	2.5,	-- average duration of activity in ms (10)
	2.5,	-- deviation of start time in ms
	0.1,	-- deviation of duration in ms
	1.2e-3)	-- peak conductivity in units of uS (6e-4)
CE:set_synapse_handler(syn_handler)

-- domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)


syn_handler:print_synapse_statistics(2)

