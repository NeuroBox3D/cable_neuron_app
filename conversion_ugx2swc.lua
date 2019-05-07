--------------------------------------------------------------
-- This script will create a .ugx file from NeuGen output	--
--------------------------------------------------------------

ug_load_script("ug_util.lua")

ugxName = util.GetParam("-grid", "testNetwork_splitSyn.ugx")

InitUG(3, AlgebraType("CPU", 1));
dom1d = util.CreateDomain(ugxName, 0)

nid = innermost_neuron_id_in_subset("SOMA__L5A_PYRAMIDAL", dom1d:subset_handler())
ind = util.GetParamNumber("-ni", nid)
swcName = util.GetParam("-out", "testNeuron.swc")
scale = util.GetParamNumber("-scale", 1e6)


print("Extracting neuron with index " .. ind .. " from network.")
save_neuron_to_swc(ugxName, ind, swcName, scale)
