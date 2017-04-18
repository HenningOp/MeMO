# MeMO
A sensitivity driven optimisation tool for hydrological models
MeMO - A sensitivity driven meta-model optimisation tool

MeMO is designed to optimise hydrological models comprising a high number of parameters. 
The limitation to hydrological models is certainly not strict. Basic concept is a model that
calculates series of data and an observed series that the simulation can be compared to.
MeMO will try to identify which parameter of the model affects which component of the simulated
data series. The components are: Deviation of Average (BIAS), Deviation of Standard-deviation (Reactivity),
and linear correlation coefficient (Timing). Since MeMO is designed for hydrologic models the
components can be (straightforwardly) described as flood volume, relation of flood peaks and
the timing of simulated floods.



**System requirements**:

- OS:
	Linux
	Windows 7/8/10	64Bit
- Script Language:
	Python	Version	2.7.10
- Required Python Addins (Minimum):
	numpy	Version	1.11.0
	


**Installation**:

At the current state (BETA) MeMO is not availiable in the python package index so I evaded the
compilation of a setup file. So please install it the (very) plain way:

1. Unzip/Copy script to destination folder
2. Script can be started without further installation




**Usage**:

To use MeMO to optimise your model you have built three instances:

1. 	An array/list of a first parameter guess (can be random) par0 comprising all parameters that shall
	be used in calibration. All elements are floating values.
2. 	Two arrays containing the parameter contraints (upper and lower). All elements are floating values
3. 	A function that gives estimated parameters to the model, runs the model and returns a list
	with calculated components of the Kling-Gupta Efficiency. 

An example for such a function will be given in the next section. To apply a meta-optimisation
first import MeMO to your script:

	import MeMO

or 

	import MeMO as memo

I will refer to the latter import command in the following. 
The optimisation run (opt) is instanced with the number of parameters MeMO has to hande (numPars):

	opt = memo.memo(numPars)

Next some definitions have to be made. First Parameter constraints, please note that we refer them as "right" and "left" boundaries, not as
usual with "upper" and "lower" (due to the fact that boundary values can be switched in course
of calibration):

	opt.set_rightBoundary(numpyArray_1)
	opt.set_leftBoundary(numpyArray_2)

Note that numpyArray_1 and 2 have to comprise exactly numPars elements. Defined right/left 
parameter constraints can be returned with a get statement. This might be interesting in 
case of application on unknown parameters (see next steps):

	opt.get_rightBoundary()
	opt.get_leftBoundary()

Next step is to define the function that MeMO has to optimise (please see the following section
for a detailed description of the requirements for these functions):

	opt.set_objective_function(function)

After these basic settings there are two ways to proceed the calibration. Either you start at the 
very bottom and try to identify parameter groups with the sensitivity analysis of the parameters or
you already know the model, its sensitivity and the parameter groups. If the letter case is true you
can easily give the (meta)-groups to MeMO:

	opt.g_a = listA
	opt.g_b = listB
	opt.g_r = listR

g_a is the group of parameters assigned to the deviation of standard deviations, g_b is assigned to 
the deviation of averages and g_r to the linear correlation. Note that the number of elements in 
all three groups must sum up to the defined number of parameters. All elements have to be integers, 
indicating the position in the initial parameter guess array par0.
If you have not conducted a sensitivity analysis yet and the parameter groups are unknown, call 
the MeMO implemented sensitivity analysis for single parameters:

	opt.singleSensitivity(par0, r_except = True, printSummary = True, parNames = None, plotSummary = False)
	
	where:
		> par0			List containing initial parameter guesses. (1)
		> r_except		Boolean indicating if any parameter affecting linear correlation should be assigned to g_r (Most cases affection on r is rare) ignoring significance level
		> printSummary		Boolean indicating if a .txt file should be created giving performance criteria for each parameter iteration. Filename is drawn from parNames
		> plotSummary		Boolean indicating if a .png plot of parameter sensitivity should be created. Filenames are drawn from parNames
		> parNames		List/Array of strings, used for plot labels and filenames. If parNames is None parNames are replaced by integers representing parameter position in par0 list.

This command starts a one-factor-at-a-time sensitivity analysis. Each parameter in the par0 list
will be varied in its (beforehand) given constraints. All other parameters remain on the initial
guess par0, hence, (1) defines the "base performance" of the model on which the sensitivity is analysed.
Each parameter will be assigned to the meta-group (g_a, g_b or g_r) on which is has the highest influence.
If necessary the parameter is rearranged in order to direct the impact of all parameters in the
same "direction", i.e. all parameters rise the value of standard deviation, bias or correlation.
Please see the section below for further notes on the singleSensitivity analysis!

Please note that the rearrangement has to be made manually if the meta-groups are manually defined.

After the definition of the meta-groups the calibration can be started with the following command:

	opt.optimise(par0)

MeMO applies a linear adjustment of the meta-groups and its assinged parameters to the model performance.
The adjustment is made recursively nIter-times, the number of iterations is by default set to 2, but
can be altered, see following the following section (Additional settings).
The optimise command returns an array/list of optimal parameters. Order and number of elements is
identical to par0.
Performance criteria of the returned parameter set can be drawn by:

	opt.last_opt_val()
	opt.last_opt_components()




**Additional settings**:

On initialisation of MeMO some default settings are made that can be altered, but are recommended to
retained on default values.

	opt._maxIter = 2		#Number of optimisation iterations
	opt._significance = 0.1		#Significance level for single parameter sensitivity
	opt._pp				#Array of floating values in range (0,1] comprising by default 100 elements
					#	defines the steps parameter(-groups) are altered in
					#	in sensitivity analyses
	opt._overwrite = True		#Boolean if existing summaries and plots should be deleted.




**Further Notes**:

- Single sensitivity analysis
	1.	If a parameter has no impact on model performance, MeMO throws a warning that (if parNames are given)
		the respective parameter, or number of parameter, has no impact.
		The parameter has to be omitted from par0 or assigned to a meta-group manually!
		
	2.	If a parameter has an impact that on all critera that is greater than opt._significance
		MeMO throws an according warning.
		The parameter should be omitted from par0 or can be assigned to a meta-group manually!
		Alternatively the significance level can be increased.
		
	3.	At this point the outcome of the singleSensitivity routine is bound to par0. Hence, 
		results are sufficient but biased. Parameter impacts are sometimes over- or underestimated.
		Future work will be made to employ a global sensitivity analysis. Note that you can apply
		your own sensitivity analysis and give its result to model.

- Objective function
	The schematic structure of the objective function should look something like this:
		
		def objfunction(pars):
			#Instance model
			mod = initializeModell()
			
			#set given parameters into model
			mod.setControlPars(pars)
			
			#run simulation and read observed data
			sim = mod.run()
			obs = mod.getObservedData()
			
			#calculate KGE-Components
			alpha = numpy.std(sim)/numpy.std(obs)
			beta = numpy.average(sim)/numpy.average(obs)
			r = numpy.corrcoef(sim,obs)[0][0]
			
			return [alpha, beta, r]

If you're having trouble appliying MeMO to your model, or experience any bugs (I am sure thre plenty of 'em) 
feel free to contact me (henning.oppel@rub.de)!
