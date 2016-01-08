function util.CreateAndDistributeDomainFromGrid(grid, sh, numRefs, numPreRefs,
										neededSubsets, distributionMethod,
										verticalInterfaces, numTargetProcs,
										distributionLevel, wFct)

	-- create Instance of a Domain
	local dom = Domain()
	
	-- load domain
	write("Loading Domain  ... ") 
LoadDomainFromGridInMemory(dom, grid, sh)
	write("done. ")
	
	-- create Refiner
	ug_assert(numPreRefs <= numRefs, "numPreRefs must be smaller than numRefs. Aborting.");
	
	if numPreRefs > numRefs then
		numPreRefs = numRefs
	end
	
	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner = nil
	if numRefs > 0 then
		refiner = GlobalDomainRefiner(dom)
	end
	
	write("Pre-Refining("..numPreRefs.."): ")
	-- Performing pre-refines
	for i=1,numPreRefs do
		TerminateAbortedRun()
		write(i .. " ")
		refiner:refine()
	end
	write("done. Distributing...")
	-- Distribute the domain to all involved processes
	if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
		ug_error("Error while Distributing Grid. Aborting.")
	end
	write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
	
	if numRefs > 0 then
		-- Perform post-refine
		for i=numPreRefs+1,numRefs do
			TerminateAbortedRun()
			refiner:refine()
			write(i-numPreRefs .. " ")
		end
end
write("done.\n")
	
	-- Now we loop all subsets an search for it in the SubsetHandler of the domain
	if neededSubsets ~= nil then
		if util.CheckSubsets(dom, neededSubsets) == false then 
			ug_error("Something wrong with required subsets. Aborting.");
		end
	end
	
	
	--clean up
	if refiner ~= nil then
		delete(refiner)
	end
	
	-- return the created domain
	return dom
end
