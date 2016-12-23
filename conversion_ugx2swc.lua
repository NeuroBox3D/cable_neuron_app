--------------------------------------------------------------
-- This script will create a .ugx file from NeuGen output	--
--------------------------------------------------------------

ug_load_script("ug_util.lua")

ugxName = util.GetParam("-grid", "testNetwork_splitSyn.ugx")
ind = util.GetParamNumber("-ni", 0)
swcName = util.GetParam("-out", "testNeuron.swc")
scale = util.GetParamNumber("-scale", 1e6)

save_neuron_to_swc(ugxName, ind, swcName, scale)
