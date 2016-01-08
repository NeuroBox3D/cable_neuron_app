ug_load_script("ug_util.lua")

-- Gathering standard values
infile			=	util.GetParam("-infile", "grids/31o_pyramidal19aFI.CNG_with_subsets_and_diams.ugx")
--infile			=	util.GetParam("-infile", "grids/test_cell_small_ref_1.ugx")
outfile			=	util.GetParam("-outfile", "grids/test_out.ugx")

-- Init UG Algebra for use of domain
dim = util.GetParamNumber("-dim", 3)
InitUG(dim, AlgebraType("CPU", 1));
numRefs = 0
neededSubsets = {"axon", "dend", "soma"}
dom = util.CreateDomain(infile, numRefs, neededSubsets)
--util.DistributeDomain(dom, "metis")

-- Instantiate SynapseDistributor object
--synDistr = SynapseDistributor(dom, "grids/grid_out.ugx")
synDistr = SynapseDistributor(infile, outfile)

--Status of grid before placing of synapses
synDistr:print_status()
print("Place synapses uniform.")
synDistr:place_synapses_uniform(1, 1000)
synDistr:print_status()
synDistr:degenerate_uniform(0.5, 1)
synDistr:print_status()

synDistr:export_grid()

-- PARAMETER SET FOR ACTIVITY TIMING
--[[
15,					--average start time of synaptical activity in ms
1,					--average duration of activity in ms
13, 				--deviation of start time in ms
0.5 				--deviation of duration in ms
--]]

--synDistr:set_activation_timing(15, 1, 13, 0.5)
--[[
for t = 0, 35, 1 do
	print(synDistr:num_active_synapses(t))		
end
--]]

--synDistr:activity_info()