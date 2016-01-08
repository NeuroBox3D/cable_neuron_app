-- for profiler output
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")

--set_debug(debugID.MAIN, 2)

--------------------------------------------------------------------------------
-- UG4-Standard-Settings
--------------------------------------------------------------------------------

-- choose algebra
dim = util.GetParamNumber("-dim", 3)
InitUG(dim, AlgebraType("CPU", 1));


numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numRefs    = util.GetParamNumber("-numRefs",    3)

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
                        "SynapseHandler",
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
grid = "foo_clean.ugx"
grid = "test3.ugx"
--grid = "small.ugx"

SPF = SynapseHandlerFactory3d()
NETIHandler = SPF:get_neti_handler()
NETIHandler:reg_alpha_syn()
NETIHandler:reg_exp2_syn()
common:printfn("%s", " done!")
x = 0
y = 0
z = 0
time = 0 
current = 0
--status = NETIHandler:synapse_at_location(x, y, z, time)
--common:printfn("Synapse at location (%d, %d, %d) at time (%d)? %s", x, y, z, time, tostring(status))

-- get the NeuronalTopologyImporter Plugin Geometry Importer Provider
gip = NeuronalTopologyImporterProvider()

-- get the default Geometry Importer
gi = gip:getDefaultNeuronalTopologyImporter()

-- setup joining criteria
joiningCriteria = {
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

dt = util.GetParamNumber("-dt", 0.01)
time = 10 -- in ms
timefaktor = 1/dt
-- 100 steps = 1ms
numbersteps = util.GetParamNumber("-nSteps", timefaktor*time)
StartValue = -65.0

-- choose number of time steps
--NumPreTimeSteps = util.GetParamNumber("-numPreTimeSteps", 1)
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", numbersteps)


print(" chosen parameters:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. grid)


--------------------------------------------------------------------------------
-- Synapse distributions via plugin by Lukas Reinhardt
--------------------------------------------------------------------------------

Flux_ac = 1e-5


--------------------------------------------------------------------------------
-- Settings for HH-Fluxes different fluxes are seperated by ","
--------------------------------------------------------------------------------
-- constants
InfluxValue = {8e-2} -- in C/m^2/ms
--0.125, 0.25, 0
InfluxPlacex = {0} --{-2.0e-5}-- in m 
InfluxPlacey = {0} --{0.0}-- in m
InfluxPlacez = {0} --{0.0}-- in m
InfluxStart = {0*timefaktor}--{10*timefaktor}		-- after 10 ms
InfluxDuration = {1*timefaktor}	-- duration of 2 ms
Flux_ac = 1e-9

--------------------------------------------------------------------------------
-- testing if all influx values are right
--------------------------------------------------------------------------------
Number_of_Influxes = table.getn(InfluxPlacex)

if 	(table.getn(InfluxValue) == Number_of_Influxes and 
	table.getn(InfluxPlacex) == Number_of_Influxes and
	table.getn(InfluxPlacey) == Number_of_Influxes and
	table.getn(InfluxPlacez) == Number_of_Influxes and
	table.getn(InfluxStart) == Number_of_Influxes and
	table.getn(InfluxDuration) == Number_of_Influxes)
then print("Influxes are okay!")
else print("Mistake in influx settings"); exit()
end
--------------------------------------------------------------------------------
-- Settings for HH-Constants
--------------------------------------------------------------------------------
-- in c/m^2/mV/ms
g_Na = 1.2e-3
g_K  = 0.36e-3
g_L  = 0.003e-3
-- accuracy because of gating params
ac = 1e-6
--------------------------------------------------------------------------------
-- Settings for Dendrit
--------------------------------------------------------------------------------
Diameter = 1.0e-6 -- in m
spec_cap = 1.0e-5 -- in C/mV/m^2
spec_res = 1.0e6 -- in mV ms m / Q


--------------------------------------------------------------------------------
-- User-Inputs for Data output
--------------------------------------------------------------------------------
-- Help-Vars for Output
XWerte = {}
YWerte = {}
ZWerte = {}
Value = {}
Erg = {}
it = 0
-- Point of Output
xStelle = 0.0
yStelle = 0.0
zStelle = 0.0 


-- Create, Load, Refine and Distribute Domain
neededSubsets = {"axon", "dend", "soma", "Exp2Syn", "AlphaSyn"}
dom = util.CreateAndDistributeDomain(grid, numRefs, numPreRefs, neededSubsets, "metis")

--print("Saving parallel grid layout")
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 1e-5)

-- create Approximation Space
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("ca", "Lagrange", 1)

-- gating functions for HH-Fluxes
approxSpace:print_layout_statistic()
approxSpace:print_statistic()




---------------------------------------------------
--Setup for injection flux (all dimensions)
---------------------------------------------------


-- injektion fluxes for all dimensions
function injection3d(t, x, y, z)
	--print("x: "..x.."y: "..y.."z: "..z.." time: "..time)
	for i = 1, Number_of_Influxes, 1 do
		if  t >= InfluxStart[i]*dt
			and t <= (InfluxStart[i]*dt + InfluxDuration[i]*dt)
			and math.abs(x - InfluxPlacex[i]) < Flux_ac
			and math.abs(y - InfluxPlacey[i]) < Flux_ac
			and math.abs(z - InfluxPlacez[i]) < Flux_ac
		then
		--print("flux is done")
		return InfluxValue[i]
		else
		return 0.0
		end
	end
end

function injection2d(t, x, y)
	--print("x: "..x.." time: "..time)
	for i = 1, Number_of_Influxes, 1 do
		if  t >= InfluxStart[i]*dt
			and t <= (InfluxStart[i]*dt + InfluxDuration[i]*dt)
			and math.abs(x - InfluxPlacex[i]) < Flux_ac
			and math.abs(y - InfluxPlacey[i]) < Flux_ac
		then
		--print("flux is done")
		return InfluxValue[i]
		else
		return 0.0
		end
	end
end


function injection1d(t, x)
	--print("x: "..x.." time: "..time)
	for i = 1, Number_of_Influxes, 1 do
		if  t >= InfluxStart[i]*dt
			and t <= (InfluxStart[i]*dt + InfluxDuration[i]*dt)
			and math.abs(x - InfluxPlacex[i]) < Flux_ac
		then
		--print("flux is done")
		return InfluxValue[i]
		else
		return 0.0
		end
	end
end

function Injektion_Bound(x, y, z, t)
	--print("flux function is used")
	for i = 1, Number_of_Influxes, 1 do
		if  t >= InfluxStart[i]*dt
			and t <= (InfluxStart[i]*dt + InfluxDuration[i]*dt)
			and math.abs(x - InfluxPlacex[i]) < Flux_ac
			and math.abs(y - InfluxPlacey[i]) < Flux_ac
			and math.abs(z - InfluxPlacez[i]) < Flux_ac
		then
		--print("flux is done")
		return true, InfluxValue[i]
		else
		return false, 0.0
		end
	end
end


-- Flux with lua callback
--injection_flux = LuaFunctionNumber()
--injection_flux:set_lua_callback("injection"..dim.."d", (1+dim))

ena = 63.5129
ek = -74.1266

-- Adding hodgkin and huxley fluxes
HH = ChannelHHNernst("v, k, na", "axon, dend, soma")

--constructor creates every needed concentration out of added Channels from Channel list
VMD = VMDisc("axon, dend, soma, Exp2Syn", approxSpace)
VMD:set_diameter(Diameter)
--VMD:set_diff_coeffs({1.0e-12, 1.0e-12, 2.2e-13}) --m^2/ms
VMD:add_channel(HH)
--VMD:set_synapse_provider(NETIHandler)
print("Setup Presyn Subset...")
NETIHandler:set_presyn_subset("PreSynapse")
print("Setup VMDisc...")
NETIHandler:set_vmdisc(VMD) -- segfaults, investigate why. TODO

--[[
print("Setting influxes")
for i = 1, Number_of_Influxes, 1 do
	VMD:set_influx(InfluxValue[i], InfluxPlacex[i], InfluxPlacey[i], InfluxPlacez[i], InfluxStart[i]*dt, (InfluxStart[i]*dt + InfluxDuration[i]*dt))
end
--]]

-- Set dendritic geo Values


--VMD:set_spec_cap(spec_cap)
--VMD:set_spec_res(spec_res)
--VMD:set_influx_ac(Flux_ac)
--building of new Gridfunction is needed for the added channels so all channels have to be added before this




approxSpace:init_levels();
approxSpace:init_surfaces();
approxSpace:init_top_surface();


------------------------------------------------------------
-- creating grid function
------------------------------------------------------------
-- get grid function
u = GridFunction(approxSpace)
time = 0.0
step = 0

-- Startvalues can now beiing iterpoplate VMD construktor build unknowns in dof
-- Getting new ApproxSpace and writting new Gridfunction



-------------------------------------
--Setting all Startvalues
-------------------------------------

-- set initial value
Interpolate(-65.0, u, "v", time)
Interpolate(54.4, u, "k", time);
Interpolate(10.0, u, "na", time);
Interpolate(5e-5, u, "ca", time)

-- Init interpolates all needed concentrations
-- has to be written in IChannel:init()
--HH:init(dt, u)
--HHa:init(dt, u)



-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(VMD)
--domainDisc:add(HH)
--domainDisc:add(injection_bound)

-------------------------------------------
--  Setup Time Discretization
-------------------------------------------

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)

-------------------------------------------
--  Algebra
-------------------------------------------
-- create operator from discretization
op = AssembledOperator(timeDisc)
op:init()
-------------------------------------------
-- linear solver
-------------------------------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(1.0)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
--bgs = BackwardGaussSeidel()
ilu = ILU()

lu = LU()
--slu = SuperLU()
--ilu:set_debug(dbgWriter)
--ilut = ILUT()

-- create GMG ---
-----------------

	-- Base Solver
	--baseConvCheck = ConvCheck()
	--baseConvCheck:set_maximum_steps(500)
	--baseConvCheck:set_minimum_defect(1e-8)
	--baseConvCheck:set_reduction(1e-30)
	--baseConvCheck:set_verbose(false)
	base = LU()
	--base = LinearSolver()
	--base:set_convergence_check(baseConvCheck)
	--base:set_preconditioner(jac)
	--base = sgs
	-- Gemoetric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	gmg:set_base_solver(base)
	gmg:set_smoother(ilu)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_debug(dbgWriter)


-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(2000)
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-08)
convCheck:set_verbose(false)

-- create Linear Solver
--linSolver = LinearSolver()
--linSolver:set_preconditioner(ILU())
--linSolver:set_convergence_check(convCheck)

-- create ILU Solver
iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilu)
iluSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(ilu)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

-- linear iterative solver
linSolver = LinearSolver()
--linSolver:set_preconditioner(slu)
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)


-- convergence check
--[[
newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(30)
newtonConvCheck:set_minimum_defect(1e-21)
newtonConvCheck:set_reduction(1e-10)
]]--


newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 20, 2e-26, 1e-08)
newtonConvCheck:set_component_check("v", 1e-21, 1e-12)
--newtonConvCheck:set_component_check("Na", 1e-15, 1e-10)
--newtonConvCheck:set_component_check("K", 1e-15, 1e-10)
newtonConvCheck:set_verbose(true)
newtonLineSearch = StandardLineSearch()


-- create Newton Solver
newtonSolver = NewtonSolver()
--newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)
--newtonSolver:set_line_search(newtonLineSearch)

------------------------------------------------
--  Solve stationary problem for start solution
------------------------------------------------
--[[ this is not meaningful unless a Dirichlet bnd is set!
print("Solving stationary problem")
newtonSolver:set_linear_solver(lu)

op_stat = AssembledOperator(domainDisc)
newtonSolver:init(op_stat)
newtonSolver:prepare(u)

if newtonSolver:apply(u) == false then
	print("Newton Solver failed."); exit();
end

out = VTKOutput()
out:select_all(true)
out:print("HHKabel_start", u)

print("stationary Problem solved")
--]]
-----------------------
-- Start time stepping
-----------------------
--newtonSolver:set_linear_solver(slu)
newtonSolver:init(op)

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")


-- filename
filename = "Solvung"

-- write start solution
if (generateVTKoutput) then 
	out = VTKOutput()
	out:print(filename, u, step, time)
end

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   NumTimeSteps is " .. NumTimeSteps   )--.. ",  NumPreTimeSteps is " .. NumPreTimeSteps )
step = 0

-- store grid function in vector of  old solutions
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)
i=0

-- init of all channels
--HH:init(dt, u)



for step = 1, NumTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	--[[
	----------------------------------
	-- Ausgabe vorbereiten
	----------------------------------
	SaveVectorCSV(u, "testvector"..i.."pr"..ProcRank()..".csv")
	
	zeilen = {}
	j=0
	for line in io.lines("testvector"..i.."pr"..ProcRank()..".csv", "r") do
		zeilen[j] = line
		j=j+1
	end
	---------------------------------------------------------------------------
	--Ausgabe parsen
	---------------------------------------------------------------------------
	k=1    
k=1    

	while zeilen[k]~=nil do
		nBeginn, nEnde = string.find(zeilen[k], ",")
		-- X-Wert
		--print(nBeginn.." - "..nEnde)
		XWerte[k] = string.sub(zeilen[k], 0, nEnde-1)

		--print(XWert)
		if nEnde~=nil
		then
			nBeginn1, nEnde1 = string.find(zeilen[k], ",", (nEnde+1))
		--print(nBeginn1.." - "..nEnde1)
		YWerte[k] = string.sub(zeilen[k], nEnde+2, nEnde1-1)
		--print(YWert)
		end
		if nEnde1~=nil
		then
			nBeginn2, nEnde2 = string.find(zeilen[k], ",", (nEnde1+1))
		ZWerte[k] = string.sub(zeilen[k], nEnde1+2, nEnde2-1)
		end
		if nEnde2~=nil
		then
			Value[k] = string.sub(zeilen[k], nEnde2+2)
		end
		--print(Werte[k])
		--print(Value[k])
		if math.abs(XWerte[k] - xStelle) < 1e-10 and math.abs(YWerte[k]-yStelle) < 1e-10 and math.abs(ZWerte[k] - zStelle) < 1e-10
		then
			--print("Ausgabe")
			Erg[it] = Value[k]
			it=it+1
		end
		k=k+1
	end
	i=i+1
	--]]
	
	-- choose time step
	do_dt = dt
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step_elem(solTimeSeries, do_dt)


	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver prepare failed at step "..step.."."); 
		if (generateVTKoutput) then 
			out:write_time_pvd(filename, u);
		end
		exit(); 
	end 
	
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then print ("Newton solver apply failed at step "..step.."."); out:write_time_pvd(filename, u); exit(); end 
	
	-- update new time
	time = solTimeSeries:time(0) + do_dt
	
	if (generateVTKoutput) then 
		out:print(filename, u, step, time)
	end
	
	oldestSol = solTimeSeries:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)

	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

--[[
test = io.open("Ausgabe.txt","w")

for i = 0, table.getn(Erg)-1 ,1 do
	-- for Ausgabe of m, n, h write +1,+2,+3
	test:write(Erg[i])
	test:write("\n")
end
test:close("Ausgabe.txt")
--]]

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(filename, u) end
