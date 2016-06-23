ug_load_script("ug_util.lua")

-- Gathering standard values
infile			=	util.GetParam("-infile", "grids/31o_pyramidal19aFI.CNG.ugx")
--infile			=	util.GetParam("-infile", "grids/test_cell_small_ref_1.ugx")
outfile			=	util.GetParam("-outfile", "grids/test_out.ugx")

-- Init UG Algebra for use of domain
dim = util.GetParamNumber("-dim", 3)
InitUG(dim, AlgebraType("CPU", 1));
numRefs = 0
neededSubsets = {"axon", "dend", "soma"}
dom = util.CreateDomain(infile, numRefs, neededSubsets)
--util.DistributeDomain(dom, "metis")

alphasyns = AlphaSynapses(0,15)
alphasyns:set_mean_gMax(3)
alphasyns:set_dev_gMax(1)
alphasyns:set_mean_onset(2)
alphasyns:set_dev_onset(0.5)
alphasyns:set_mean_tau(1)
alphasyns:set_dev_tau(0.3)
alphasyns:set_mean_e(0.5)
alphasyns:set_dev_e(0.1)

exp2syns = Exp2Synapses(14,30)
exp2syns:set_mean_tau1(3)
exp2syns:set_dev_tau1(1)
exp2syns:set_mean_tau2(4)
exp2syns:set_dev_tau2(1)
exp2syns:set_mean_e(5)
exp2syns:set_dev_e(0.5)
exp2syns:set_mean_w(1)
exp2syns:set_dev_w(0.2)



-- Instantiate SynapseDistributor object
--synDistr = SynapseDistributor(dom, "grids/grid_out.ugx")
synDistr = SplitSynapseDistributor(infile, outfile, false)

--Status of grid before placing of synapses
print("Place synapses uniform.")
synDistr:place_synapses_uniform(alphasyns:get_synapses())
synDistr:place_synapses_uniform(exp2syns:get_synapses())

--synDistr:print_status()
--synDistr:degenerate_uniform(0.5, 1)
synDistr:print_status()

print(synDistr:export_grid())

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
