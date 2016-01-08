--[[
--
-- ******************************************
-- *** Synapse Provider example script
-- ******************************************
--
-- provides synapses for grids
--]]

-- build configuration, utility load and assert plugin loaded
PrintBuildConfiguration()

-- splits a path and takes the last field (filename)
function string:splitPath(sep, take)
   local sep, fields = sep or "/", {}
   local pat = string.format("([^%s]+)", sep)
   self:gsub(pat, function(z) fields[#fields+1] = z end)
   return fields[take or #fields]
end

-- load ug script util
ug_load_script("../scripts/ug_util.lua")
ug_load_script("../apps/synapse_provider/util.lua")
AssertPluginsLoaded(
                     {
                        "SynapseProvider",
                        "NeuronalTopologyImporter"
                        -- extend required plugins list here 
                        -- SynapseDistributor,
                        -- NeuronalTopologyImporter,
                        -- etc ...
                     }
                   )

-- cli arguments
file = util.GetParam("-file", "synapses.hoc")
grid = file:splitPath(".", 1) .. ".ugx"

-- get the NeuronalTopologyImporter Plugin Geometry Importer Provider
gip = NeuronalTopologyImporterProvider()

-- get the default Geometry Importer
gi = gip:getDefaultNeuronalTopologyImporter()

-- setup joining criteria
joiningCriteria = {
   "apic",
   "axon",
   "soma",
   "dend",
   "Exp2Syn",
   "AlphaSyn"
}

for _, v in pairs(joiningCriteria) do
   gi:add_joining_criteria(v)
end

-- import the geometry and generate the grid with the importer for hoc/swc files
import_status = gi:import_geometry(file, ""), "GeometryImport"
common:printfn("Status: %s", tostring(import_status))

-- create the Synapse Provider Factory for 3d 
SPF = SynapseProviderFactory3d()
NETIProvider = SPF:get_neti_provider()
common:printfn("%s", " done!")
x = 0
y = 0
z = 0
time = 0 
current = 0
status = NETIProvider:synapse_at_location(x, y, z, time)
common:printfn("Synapse at location (%d, %d, %d) at time (%d)? %s", x, y, z, time, tostring(status))

--- create domain from grid in memory
dim = 3
ug_load_script("../scripts/ug_util.lua")
InitUG(dim, AlgebraType("CPU", 1));
dom = Domain()
LoadDomain(dom, grid)
--grid = gi:get_grid_p()
--sh = gi:get_sh_p()
numRefs = 1
numPreRefs = 0
neededSubsets = {"axon", "soma", "dend"}
--LoadDomainFromGridInMemory(dom, grid, sh)
dom = util.CreateAndDistributeDomain(grid, numRefs, numPreRefs, neededSubsets, "metis")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("ca", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()
VMD = VMDisc("axon", approxSpace)
--Diameter = 10
--common:printf("Setting diameter and NETIProvider...")
--VMD:set_diameter(Diameter)
--NETIProvider:set_diameter_attachment(VMD)
VMD:set_synapse_provider(NETIProvider)
NETIProvider:set_vmdisc(VMD)
common:printf("... done!")
status = NETIProvider:synapse_at_location(x, y, z, time)
common:printfn("Synapse at location (%d, %d, %d) at time (%d)? %s", x, y, z, time, tostring(status))

