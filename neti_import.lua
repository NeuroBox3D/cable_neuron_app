--------------------------------------------------------------
-- This script will create a .ugx file from NeuGen output	--
--------------------------------------------------------------

ug_load_script("ug_util.lua")

-- set debug level
GetLogAssistant():set_debug_level("NETI_DID.import_geometry", 1)
GetLogAssistant():set_debug_level("NETI_DID.generate_grid", 1)

-- get a geometry importer provider
gip = NeuronalTopologyImporterProvider()

-- get the default geometry importer
gi = gip:getDefaultNeuronalTopologyImporter()

-- base name
baseName = util.GetParam("-name", "testNetwork")
-- import
method = util.GetParam("-method", "txt")

if method == "txt" then
gi:import_txt(baseName.."_secs.txt", baseName.."_connex.txt", baseName.."_synapses.txt", baseName.."_identifier.txt")
elseif method == "swc" then
gi:import_geometry_and_generate_grid(baseName, "swc")
else 
print("Unknown Method:"..method)
end



