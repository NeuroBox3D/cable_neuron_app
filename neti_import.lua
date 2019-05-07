--------------------------------------------------------------
-- This script will create a .ugx file from NeuGen output	--
--------------------------------------------------------------

ug_load_script("ug_util.lua")

-- set debug level
GetLogAssistant():set_debug_level("NETI_DID.import_geometry", 6)
GetLogAssistant():set_debug_level("NETI_DID.generate_grid", 6)

-- get the default geometry importer
neti = NeuronalTopologyImporter()

-- base name
baseName = util.GetParam("-name", "testNetwork")
-- import
method = util.GetParam("-method", "txt")

if method == "hoc" or method == "ngx" then
	neti:add_joining_criterion("soma")
	neti:add_joining_criterion("dend")
	neti:add_joining_criterion("apic")
	neti:add_joining_criterion("axon")
end

neti:import_geometry_and_generate_grid(baseName, method)

