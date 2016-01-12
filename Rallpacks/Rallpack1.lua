-- Haut nach oben ab beim Anfang


-- for profiler output
SetOutputProfileStats(true)

ug_load_script("ug_util.lua")

--set_debug(debugID.MAIN, 2)

--------------------------------------------------------------------------------
-- UG4-Standard-Settings
--------------------------------------------------------------------------------
dim = util.GetParamNumber("-dim", 3)
dt = util.GetParamNumber("-dt", 1e-5)
time = 0.01
timefaktor = 1/dt
-- 100 steps = 1ms
numbersteps = util.GetParamNumber("-nSteps", timefaktor*time)
StartValue = -0.065
--------------------------------------------------------------------------------
-- Settings for HH-Fluxes different fluxes are seperated by ","
--------------------------------------------------------------------------------
-- constants
InfluxValue = {1e-10} -- in A
--(4.09897e-06, -7.83906e-08, 1.11741e-06)
--influx for 4 refines
InfluxPlacex = {0} --{-2.0e-5}-- in m 
InfluxPlacey = {0} --{0.0}-- in m
InfluxPlacez = {0} --{0.0}-- in m
InfluxStart = {0*timefaktor}--{10*timefaktor}		-- after 10 ms
InfluxDuration = {1e-3*timefaktor}	-- duration of 1 ms
Flux_ac = 1e-5

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
-- in S/m^2
g_Na = 1.2e3
g_K  = 360.0
g_L  = 3.0
-- accuracy because of gating params
ac = 1e-9
--------------------------------------------------------------------------------
-- Settings for Dendrit
--------------------------------------------------------------------------------
Diameter = 1.0e-6 -- in m
spec_cap = 1.0e-2 -- in F/m^2
spec_res = util.GetParamNumber("-spec_res", 1.0) -- in Ohm m


--------------------------------------------------------------------------------
-- User-Inputs for Data output
--------------------------------------------------------------------------------
-- true output is written false not 
ausgabe = util.GetParamBool("-output", true)

-- Help-Vars for Output
XWerte = {}
YWerte = {}
ZWerte = {}
Value = {}
Erg = {}
it = 0
-- Point of Output
--0.00035087, -0.000130453, 7.6404e-05)
xStelle = 0--0.00035087
yStelle = 0---0.000130453
zStelle = 0--7.6404e-05


-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

if 	dim == 3 then gridName = util.GetParam("-grid", "../../calcium_channels_mitochondria/grids/Ralls/1mm1000comps.ugx")
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    0)

-- choose number of time steps
--NumPreTimeSteps = util.GetParamNumber("-numPreTimeSteps", 1)
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", numbersteps)


print(" chosen parameters:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)


-- Create, Load, Refine and Distribute Domain
neededSubsets = {"subset"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets, "metis")

print("Saving parallel grid layout")
--SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 1e-5)

-- create Approximation Space
--print("Create ApproximationSpace needs to be somewhere else")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("k", "Lagrange", 1)
approxSpace:add_fct("na", "Lagrange", 1)
approxSpace:add_fct("ca", "Lagrange", 1)



--constructor creates every needed concentration out of added Channels from Channel list
CE = CableEquation("subset")
CE:set_diameter(Diameter)
CE:set_diff_coeffs({1.0e-9, 1.0e-9, 2.2e-10}) --m^2/s

---[[
print("Setting influxes")
for i = 1, Number_of_Influxes, 1 do
	CE:set_influx(InfluxValue[i], InfluxPlacex[i], InfluxPlacey[i], InfluxPlacez[i], InfluxStart[i]*dt, (InfluxDuration[i]*dt))
end
--]]

-- Set dendritic geo Values


CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)
--CE:set_influx_ac(Flux_ac)
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

-- Startvalues can now beiing iterpoplate CE construktor build unknowns in dof
-- Getting new ApproxSpace and writting new Gridfunction



-------------------------------------
--Setting all Startvalues
-------------------------------------

-- set initial value
Interpolate(-0.065, u, "v", time)
Interpolate(54.4, u, "k", time);
Interpolate(10.0, u, "na", time);
Interpolate(5e-5, u, "ca", time)

-- Init interpolates all needed concentrations
-- has to be written in ICableMembraneTransport:init()
--HH:init(dt, u)
--HHa:init(dt, u)



-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(CE)


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
jac:set_damp(0.8)
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
convCheck:set_reduction(1e-04)
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


newtonConvCheck = CompositeConvCheck3dCPU1(approxSpace, 100, 2e-23, 1e-14)
newtonConvCheck:set_component_check("v", 1e-20, 1e-12)
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
filename = ("Solvung_ddddd"..dt)

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
	

	--if ausgabe==true then
		----------------------------------
		-- Preparing Output
		----------------------------------
		SaveVectorCSV(u, "testvector"..i.."pr"..ProcRank()..".csv")
	
		zeilen = {}
		j=0
		for line in io.lines("testvector"..i.."pr"..ProcRank()..".csv", "r") do
			zeilen[j] = line
			j=j+1
		end
		---------------------------------------------------------------------------
		-- Parsing Output
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
	--end
	
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

print("before output")
--if ausgabe==true then
test = io.open(("Ausgabeconv_exp_".. dt ..".txt"),"w")

	for i = 0, table.getn(Erg)-1 ,4 do
		-- for Ausgabe of m, n, h write +1,+2,+3
		print("output generating")
		test:write(Erg[i])
		test:write("\n")
	end
	test:close("Ausgabeconv".. dt ..".txt")
--end

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(filename, u) end

