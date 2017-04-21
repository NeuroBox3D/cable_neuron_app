-------------------------------------------------------------------------
-- This script gives some information about neuronal networks.         --
--                                                                     --
-- author: mbreit                                                      --
-- date:   15-04-2017                                                  --
-------------------------------------------------------------------------

ug_load_script("ug_util.lua")

AssertPluginsLoaded({"cable_neuron"})

dim = 3
InitUG(dim, AlgebraType("CPU", 1));

gridName = util.GetParam("-grid", "../apps/cable_neuron_app/grids/NMO_05371_1d.ugx")

neededSubsets = {}
dom = util.CreateDomain(gridName, 0, neededSubsets)

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:init_top_surface();

ss_axon = "AXON__L4_STELLATE, AXON__L23_PYRAMIDAL, AXON__L5A_PYRAMIDAL, AXON__L5B_PYRAMIDAL"
ss_dend = "DEND__L4_STELLATE, DEND__L23_PYRAMIDAL, DEND__L5A_PYRAMIDAL, DEND__L5B_PYRAMIDAL"
ss_soma = "SOMA__L4_STELLATE, SOMA__L23_PYRAMIDAL, SOMA__L5A_PYRAMIDAL, SOMA__L5B_PYRAMIDAL"
	
CE = CableEquation(ss_axon..", "..ss_dend..", "..ss_soma, false)
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
CE:set_synapse_handler(syn_handler)

domDisc = DomainDiscretization(approxSpace)
domDisc:add(CE)

prim_iter = syn_handler:begin_OnsetPreSynapse()
prim_iter_end = syn_handler:end_OnsetPreSynapse()
numPrim = 0
while prim_iter:inequal(prim_iter_end) do
	numPrim = numPrim + 1
	prim_iter:next()
end
print("Number of OnsetPreSynapses: " .. numPrim)

sec_iter = syn_handler:begin_ThresholdPreSynapse()
sec_iter_end = syn_handler:end_ThresholdPreSynapse()
numSec = 0
while sec_iter:inequal(sec_iter_end) do
	numSec = numSec + 1
	sec_iter:next()
end
print("Number of ThresholdPreSynapses: " .. numSec)
