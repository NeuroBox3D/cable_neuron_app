--------------------------------------------------------------------------------
-- This script will extract one cell from a NeuGen-created network.	          --
-- It is chosen to be the innermost cell of a specific subset.                --
-- The cell is then exported as separate SWC file, retaining all diameter,    --
-- but losing all synapse information.                                        --
--                                                                            --
-- Author:  Markus Breit                                                      --
-- Date:    2017-03-25                                                        --
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"cable_neuron"})

ugxName = util.GetParam("-grid", "cable_neuron_app/grids/testNetwork.ugx")

InitUG(3, AlgebraType("CPU", 1))
dom1d = util.CreateDomain(ugxName, 0)

nid = innermost_neuron_id_in_subset("SOMA__L5A_PYRAMIDAL", dom1d:subset_handler())
ind = util.GetParamNumber("-ni", nid)
swcName = util.GetParam("-out", "testNeuron.swc")
scale = util.GetParamNumber("-scale", 1e6)


print("Extracting neuron with index " .. ind .. " from network.")
save_neuron_to_swc(ugxName, ind, swcName, scale)
