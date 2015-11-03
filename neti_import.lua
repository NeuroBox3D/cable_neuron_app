--------------------------------------------------------------
-- This script will create a .ugx file from NeuGen output	--
--------------------------------------------------------------

ug_load_script("ug_util.lua")

-- get a geometry importer provider
gip = NeuronalTopologyImporterProvider()

-- get the default geometry importer
gi = gip:getDefaultNeuronalTopologyImporter()

-- base name
baseName = util.GetParam("-name", "testNetwork")
-- import
gi:import_txt(baseName.."_secs.txt", baseName.."_connex.txt", baseName.."_synapses.txt", baseName.."_identifier.txt")