import numpy as np
import copy, os
import matplotlib.pyplot as plt

class memo(object):
    """Meta-model optimisation tool, previously named cal 'em all"""

    def __init__(s, nPars, incrementals = 0.01, nIter = 2, SingleSig = 0.1):
        """Instance calibration object"""
        #Parameter boundaries and numbers
        s._left = np.zeros(nPars)
        s._right = np.zeros(nPars)
        s._nP = nPars

        #Meta-Parameters
        s._pp = np.arange(0., 1., incrementals)

        #Remaining stuff
        s._maxIter = nIter
        s._significance = SingleSig
        s._overwrite = True

        #Objective function
        s._function = None

        #Meta-parameter groups
        s.g_a = []
        s.g_b = []
        s.g_r = []

        #Storage for performance criteria
        s.__alpha = 0.
        s.__beta = 0.
        s.__r2 = 0.
        s.__kge = 0.

    def set_objective_function(s, function):
        """Define function to optimise"""
        s._function = function
        
    def set_leftBoundary(s, lowerBoundary):
        """Define (left/lower) parameter constraints"""
        if len(lowerBoundary) == s._nP:
            s._left = lowerBoundary
        else:
            raise ValueError('Dimension Mismatch of given Number of parameters and given Boundaries')
        
    def set_rightBoundary(s, upperBoundary):
        """Define (right/upper) parameter constraints"""
        if len(upperBoundary) == s._nP:
            s._right = upperBoundary
        else:
            raise ValueError('Dimension Mismatch of given Number of parameters and given Boundaries')

    def get_leftBoundary(s):
        """Return (left/lower) parameter constraints"""
        return s._left

    def get_rightBoundary(s):
        """Return (right/upper) parameter constraints"""
        return s._right
                          

    def singleSensitivity(s, r_except = True, printSummary = True, parNames = None, plotSummary = True):
        """Perform a sensitivity analysis for given parameters"""
        
        #Initial parameter setting to mid of parameter range
        cPars = np.add(np.multiply(np.subtract(s._right, s._left), 0.5), s._left)
        groups = [[],[],[]]

        for i, par in enumerate(cPars):
            a, b, r = [],[],[]
            #Iterate meta-par
            for ppx in s._pp:
                #Manipulate parameter array
                cPars[i] = (s._right[i]-s._left[i]) * ppx + s._left[i]

                #Calculate performance
                alpha, beta, r2 = s._function(cPars)

                #Save performance
                a.append(alpha)
                b.append(beta)
                r.append(r2)
            
            #Analysis and asignment
            ranges = np.array([max(a) - min(a), max(b) - min(b), max(r) - min(r)])

            #Test if parameter has any impact on model
            if sum(ranges) == 0.0:
                if parNames == None:
                    print 'WARNING: Parameter No.' + str(i) + ' has no impact on model performance...'
                else:
                    print 'WARNING: Parameter ' + parNames[i] + ' has no impact on model performance...'
            else:
                if len(ranges[np.greater(ranges, s._significance)]) == 3:
                    #Exception 1: "Superparameter" remain unattended
                    if parNames == None:
                        print 'WARNING: Parameter No.' + str(i) + ' has significant impact on all performance criteria...'
                    else:
                        print 'WARNING: Parameter ' + parNames[i] + ' has significant impact on all performance criteria...'
                    superpar = True
                elif ranges[2] > 1.5 * s._significance and r_except == True:
                    #Exception 2: Any parameter affecting r2 is assigned to g_r2
                    groups[2].append(i)
                    superpar = False
                    #Alignement parameters
                    mini = np.argmin(r)
                    maxi = np.argmax(r)
                    array = r
                else:
                    #Regular case: Maximum impact defines group assignement
                    groups[np.argmax(ranges)].append(i)
                    superpar = False
                    #Alignment parameters
                    if np.argmax(ranges) == 0:
                        mini = np.argmin(a)
                        maxi = np.argmax(a)
                        array = a
                    elif np.argmax(ranges) == 1:
                        mini = np.argmin(b)
                        maxi = np.argmax(b)
                        array = b
                    elif np.argmax(ranges) == 2:
                        mini = np.argmin(r)
                        maxi = np.argmax(r)
                        array = r

                #Alignement
                if superpar == False:
                    right = copy.deepcopy(s._right[i])
                    left = copy.deepcopy(s._left[i])
                    #Case 1: Minimum is not at parameter boundaries
                    if a[mini] != a[0] and a[mini] != a[-1]:
                        s._left[i] =  (right-left) * s._pp[mini] + left
                    #Case 2: Maximum is not at parameter boundaries
                    elif a[maxi] != a[0] and a[maxi] != a[-1]:
                        s._right[i] = (right-left) * s._pp[maxi] + left
                    #Case 3: Par-Movement left to right decreases performance
                    if mini > maxi:
                        right = copy.deepcopy(s._right[i])
                        left = copy.deepcopy(s._left[i])
                        s._right[i] = copy.deepcopy(left)
                        s._left[i] = copy.deepcopy(right)
                    

                #Write summary to textfile
                if printSummary == True:
                    if parNames == None:
                        j = 0
                        out = 'sensitivity_par-' + str(i) + '.txt'
                        if s._overwrite == False:
                            while os.path.isfile(out) == True:
                                j += 1
                                out = 'sensitivity_par-' + str(i) + '_' + str(j) + '.txt'
                        f = open(out, 'w')
                    else:
                        j = 0
                        out = 'sensitivity_' + parNames[i] + '.txt'
                        if s._overwrite == False:
                            while os.path.isfile(out) == True:
                                j += 1
                                out = 'sensitivity_' + parNames[i] + '_' + str(j) + '.txt'
                        f = open(out, 'w')
                        
                    f.write('PP\talpha\tbeta\tr2\n')
                    for ppx, ia, ib, ir in zip(s._pp, a, b, r):
                        f.write(str(ppx) + '\t' + str(round(ia, 4)) + '\t' + str(round(ib, 4)) + '\t' + str(round(ir, 4)) + '\n')
                    f.close()

                #Create 3-panel plot
                if plotSummary == True:
                    fig = plt.figure(figsize = (7,10))
                    #Alpha
                    ax1 = fig.add_subplot(311)
                    ax1.plot(s._pp, a, color = 'b')
                    #Beta
                    ax2 = fig.add_subplot(312)
                    ax2.plot(s._pp, b, color = 'b')
                    #r2
                    ax3 = fig.add_subplot(313)
                    ax3.plot(s._pp, r, color = 'b')
                 
                    #Annotations
                    if parNames != None:
                        ax3.set_xlabel(r'rel. Parameter ' + parNames[i] + ' [-]', fontsize = 'medium')
                    else:
                        ax3.set_xlabel(r'rel. Parameter [-]', fontsize = 'medium')

                    ax1.set_ylabel(r'$\alpha$ [-]', fontsize = 'medium')
                    ax2.set_ylabel(r'$\beta$ [-]', fontsize = 'medium')
                    ax3.set_ylabel(r'$r^2$ [-]', fontsize = 'medium')

                    if ranges[0] < 0.02:
                        ax1.set_ylim(round(min(a)-0.05, 1), round(max(a)+0.05, 1))
                    if ranges[1] < 0.02:
                        ax2.set_ylim(round(min(b)-0.05, 1), round(max(b)+0.05, 1))
                    if ranges[2] < 0.02:
                        ax3.set_ylim(round(min(r)-0.05, 1), round(max(r)+0.05, 1))
                        

                    #Save figure
                    if parNames == None:
                        out = 'par_' + str(i)  + '.png'
                        if s._overwrite == False:
                            j = 0
                            while os.path.isfile(out) == True:
                                j += 1
                                out = 'par_' + str(i)  + '_' + str(j) + '.png'
                    else:
                        out = 'par_' + parNames[i]  + '.png'
                        if s._overwrite == False:
                            j = 0
                            while os.path.isfile(out) == True:
                                j += 1
                                out = 'par_' + parNames[i]  +  + '_' + str(j) + '.png'
                    plt.savefig(out)
                    plt.close()
                    

        return groups
        

    def __groupSensitivity(s, gp, group):
        """Calculate impact of meta-parameter-groups"""

        #Copy original parameter array
        cPars = copy.deepcopy(gp)

        a, b, r = [],[],[]
        #Iterate meta-pars
        for ppx in s._pp:

            #Manipulate parameter array
            for pos in group:
                cPars[pos] = (s._right[pos]-s._left[pos]) * ppx + s._left[pos]

            #Give parameters to function
            alpha, beta, r2 = s._function(cPars)

            #Save performance
            a.append(alpha)
            b.append(beta)
            r.append(r2)
            
        return [a, b, r]
    
        
    def optimize(s, pars):
        """Start meta-model calibration"""

        #Test if meta-parameter-groups are defined
        if len(s.g_a) > 0 or len(s.g_b) >0 or len(s.g_r) > 0:
            #Groups are defined
            if sum([len(s.g_a), len(s.g_b), len(s.g_r)]) != s._nP:
                ##Number of parameter assigned to groups does not match number of parameters
                raise IndexError('Number of given meta-group-indexes does not match number of parameters')
            elif len(s.g_a)*len(s.g_b)*len(s.g_r) == 0:
                ##One or more group is emtpy
                raise IndexError('A meta-parameter group is empty')            
        else:
            #Groups are empty
            print '... assorting parameter groups'
            s.g_a, s.g_b, s.g_r = singleSensitivity(r_except = True, printSummary = False, parNames = None, printGraphics = False)

        
        #Repeat calibration maxIter-times
        for xy in range(0, s._maxIter):

            #Instance storages
            s.__a, s.__b, s.__r = [],[],[]
            s.__ALPHA, s.__BETA, s.__R2 = [],[],[]

            #Successively optimise performance
            for j in range(0, 4):

                #j=0 -> beta, j=1 -> alpha, j=2 -> r2, j=3 -> readjustment
                if j < 3:
                    #Perform Group-sensitivity and adjust meta-parameter
                    if j == 0:
                        s.__b = s.__groupSensitivity(pars, s.g_b)[1]
                        px = max(np.interp(1., s.__b, s._pp), 0.)
                        posis = s.g_b
                    elif j == 1:
                        s.__a = s.__groupSensitivity(pars, s.g_a)[0]
                        px = max(np.interp(1., s.__a, s._pp), 0.)
                        posis = s.g_a
                    elif j == 2:
                        s.__r = s.__groupSensitivity(pars, s.g_r)[2]
                        px = s._pp[np.argmax(s.__r)]
                        posis = s.g_r

                    #Manipulate parameters
                    for pos in posis:
                        pars[pos] = (s._right[pos]-s._left[pos]) * px + s._left[pos]

                    #Give parameters to model
                    s.__alpha, s.__beta, s.__r2 = s._function(pars)

##                    print s.__alpha, s.__beta, s.__r2
                    
                    s.__ALPHA.append(s.__alpha)
                    s.__BETA.append(s.__beta)
                    s.__R2.append(s.__r2)

                #Readjustment with oversteering
                else:
                    ##Actual performance
                    dac = 1-s.__alpha
                    dbc = 1-s.__beta
                    drc = 1-s.__r2
                    ##Best performance so far
                    dab = np.min(np.abs(np.subtract(1, np.array(s.__ALPHA))))
                    dbb = np.min(np.abs(np.subtract(1, np.array(s.__BETA))))
                    drb = np.min(np.abs(np.subtract(1, np.array(s.__R2))))
                    ##Difference
                    deltas = [dbc - dbb, dac - dab, drc - drb]
                    ##KGE so far
                    s.__kge = 1-np.sqrt(np.power(1-s.__alpha,2)+np.power(1-s.__beta,2)+np.power(1-s.__r2,2))
                    ##Save parameters before readjustment
                    save_pars = copy.deepcopy(pars)
                    save_params = [copy.deepcopy(s.__alpha), copy.deepcopy(s.__beta), copy.deepcopy(s.__r2)]
                    
                    ##Calculate oversteering
                    targ = [deltas[0]+1, deltas[1]+1, 1]
                    ##Repeat adjustment including oversteering
                    for j in range(0, 2):
                        if j == 0:
                            s.__b = s.__groupSensitivity(pars, s.g_b)[1]
                            px = max(np.interp(1., s.__b, s._pp), 0.)
                            posis = s.g_b
                        elif j == 1:
                            s.__a = s.__groupSensitivity(pars, s.g_a)[0]
                            px = max(np.interp(1., s.__a, s._pp), 0.)
                            posis = s.g_a
                        elif j == 2:
                            s.__r = s.__groupSensitivity(pars, s.g_r)[2]
                            px = s._pp[np.argmax(s.__r)]
                            posis = s.g_r

                        for pos in posis:
                            pars[pos] = (s._right[pos]-s._left[pos]) * px + s._left[pos]

                        s.__alpha, s.__beta, s.__r2 = s._function(pars)

                    #Test if KGE has improved
                    if s.__kge > 1-np.sqrt(np.power(1-s.__alpha,2)+np.power(1-s.__beta,2)+np.power(1-s.__r2,2)):
                        #False: Restore previous parameters and performance values
                        pars = copy.deepcopy(pars)
                        s.__alpha = copy.deepcopy(save_params[0])
                        s.__beta = copy.deepcopy(save_params[1])
                        s.__r2 = copy.deepcopy(save_params[2])
                    else:
                        #True: Update KGE
                        s.__kge = 1-np.sqrt(np.power(1-s.__alpha,2)+np.power(1-s.__beta,2)+np.power(1-s.__r2,2))

        return pars

    def last_opt_val(s):
        """Return last KGE value of optimisation"""
        return s.__kge
    
    def last_opt_components(s):
        """Return last KGE-componentes of optimisation"""
        return [s.__alpha, s.__beta, s.__r2]
    
