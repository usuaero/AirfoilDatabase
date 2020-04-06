"""
Generic Multivariable Polynomial Fit Using Least Squares Regression With
Full Control Of Interaction Terms And Constraint Capabilities

This module contains functions to calculate the polynomial coefficients for 
an arbitrary order polynomial curve fit to a dataset with an arbitrary
number of independent variables.

Routine Listings
-----------------

multivariablePolynomialFit : function for calculating a curve fit to data
    with an arbitrary number of independent variables

multivariablePolynomialFunction : function for calculating a polynomial with
    an arbitrary number of independent variables

multivariableR2 : function calculating the coefficient of determination
    value, or R^2 value, of a polynomial fit of an arbitrary number of
    independent variables

multivariableRMS : function calculating an RMS (root, mean, squared) error
    and a custom RMSN (root, mean, squared, normalized) error where
    normalized means the error is divided by the mean of the absolute value
    of the dependent variables, for a multidimensional polynomial function

compose_j : function used by the multivariable series of functions that
    composes the n values of the independent variables into the counter j,
    this function can also be used for the nhat values and i

decompose_j : function used by the multivariable seris of functions that
    decomposes the counter j into the different n values for the independent
    variables, this function can also be used for i and the nhat values.
"""
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td


def multivariablePolynomialFit(Nvec, xx, yy ,interaction=True, sym=[], sym_same=[], sym_diff=[], zeroConstraints=[], constraints=[], percent=False, weighting=None, verbose=True):
    """
    inputs
    
        Nvec = list with a length V equal to the number of independent
            variables. The ith element values are integers of polynomial
            order of the ith independent variable
        x = numpy matrix of size k by V, where k is the total number of
            points in the dataset. The ith column represents the ith
            independent variable values with the rows representing the
            different data points
        y = list with a length of k with the dependent variable values
        interaction = boolean value with default set to False. This variable
            determines whether or not interaction terms are included in the
            fit function. If set to True, interaction terms up the max order
            for each independent variable are included, i.e. if Nvec = [3,2]
            then the highest interaction term included is x_1^3*x_2^2.
            Specific interaction terms can be omitted using the constraints
            input
        sym = optional list that defaults as an empty list. If used, the
            length should be V and each element should contain a boolean,
            True or False. The ith element determines if the ith independent
            variable is symmetric either even or odd, which is determined by
            the order given in Nvec. This will also remove the cooresponding
            interaction terms if they are enabled.
        sym_same = optional list that defaults as an empty list. If used,
            the entries in the list should be tuples with two integers. The
            integers represent the independent variables that the "same"
            symmetry condition will be applied. The "same" symmetry forces
            all interaction terms
        sym_diff = optional list that defaults as an empty list. 
        zeroConstraints = an optional list that defaults as an empty list.
            Entries in the list contain integer tuples of length V. The
            integer values represent the powers of the independent variables
            whose coefficient will be forced to 0 before the best fit
            calculations are performed, allowing the user to omit specific
            interaction terms or regular polynomial terms
        constraints = an optional list that defaults to an empty list.
            Entries in the list contain tuples of length 2. The first entry 
            is a list of integers that represent the powers of the
            independent variables whose coefficient will then be forced to
            be equal to the second entry in the tuple, which should be a
            float.
        percent = boolean value with default set to False. When set to True
            the least squares is performed on the percent error squared.
            This option should not be used if y contains any zero or near
            zero values, as this might cause a divide by zero error.
        weighting = optional callable function that defaults to None. If
            given, weighting should be a function that takes as arguments:
            x, y, and p where x and y are the independent and dependent
            variables defined above and p is the index representing a
            certain data point. weighting should return a 'weighting factor'
            that determines how important that datapoint is. Returning a '1'
            weights the datapoint normally.
    
    returns
    
        a = list of the polynomial coefficients. 'a' has a length equal to
            the products of the n_vec elements plus one.
            i.e.: (n_vec0+1)*(n_vec1+1)*...
        r2 = the coefficient of determination, also referred to as the R^2
            value. A value of 1 means the polynomial given by 'a' fits the
            data perfectly
    """
    # input variables and setup calculations
    ########################################################################
    x = np.copy(xx)
    y = np.copy(yy)
    # calculate number of dimensions used
    if type(Nvec) != list:
        Nvec = [Nvec]
    V = len(Nvec)
    # calculate the number of points in dataset
    k = len(y)
    # check for inconsistencies in dimensions from input variables
    if len(x.shape) == 1:
        x = np.transpose([x])
    if x.shape[1] != V: raise ValueError('Dimensions for V don\'t match between n_vec and x. Lenght of n_vec and number of columns of x should equal the number of independent variables used, V.')
    if x.shape[0] != k: raise ValueError('Number of rows between x and y don\'t agree! The number of rows should be the total number of points, k, of the dataset.')
    # calculate the length of the coefficient list to be returned
    # J = 1
    # for n in Nvec:
        # J *= n + 1
    J = calcJ(Nvec)
    # if sym wasn't given, initialize it to False values
    if type(sym) != list: sym = [sym]
    if sym == []:
        sym = [False] * V
    elif len(sym) != V:
        raise ValueError('Length of sym doesn\'t match the number of dimensions, V.')
    # create active list
    ########################################################################
    # set active to empty list
    active = []
    # loop through j values
    for j in range(J):
        # calculate the n values
        n = decompose_j(j, Nvec)
        # check if n is a constraint then continue on to the next j
        if tuple(n) in zeroConstraints: continue
        # check if j is an interaction term and if interactions aren't allowed then continue on to the next j
        if sum(n) != max(n) and not interaction: continue
        # initialize flag variable to false
        flag = False
        # loop through the sym list to find the symmetry constraints
        for count,symm in enumerate(sym):
            # check if flag has been tripped, then continue to the next j if it has
            if flag: break
            # check for a symmetry constraint
            if symm:
                # check if the order of the count-th independent variable is even
                if Nvec[count]%2 == 0:
                    # check if the n value of the count-th independent variable is odd
                    if n[count]%2 == 1:
                        flag = True
                # this else block means the order of the count-th independent variable is odd
                else:
                    # check if the n value of the count-th independent variable is even
                    if n[count]%2 == 0:
                        flag = True
        # if the flag has been tripped, skip to the next j value
        if flag: continue
        # loop through sym_same constraints
        for val in sym_same:
            # check if the n values from both variables given in val are even, then trip flag
            if n[val[0]]%2 == 0 and n[val[1]]%2 == 0:
                flag = True
            # check if the n values from both variables given in val are odd, then trip flap
            if n[val[0]]%2 == 1 and n[val[1]]%2 == 1:
                flag = True
        # loop through sym_diff constraints
        for val in sym_diff:
            # check if the n values from both variables given in val are even and odd, then trip flag
            if n[val[0]]%2 == 0 and n[val[1]]%2 == 1:
                flag = True
            # check if the n values from both variables given in val are odd and even, then trip flap
            if n[val[0]]%2 == 1 and n[val[1]]%2 == 0:
                flag = True
        # if flag hasn't been tripped, append j value onto the active list
        if not flag: active.append(j)
    #create constraints list
    ########################################################################
    con = {}
    for n,val in constraints:
        j = compose_j(n, Nvec)
        con['{}'.format(j)] = val
    #create A matrix and b vector
    ########################################################################
    # initialize A matrix
    A = np.zeros( ( len(active),len(active) ) )
    # initialize b vector
    b = np.zeros( len(active) )
    # setup progress bar for display
    if verbose: prog = oneLineProgress(len(active), msg='Building Matrix Eq for MultiPolyFit')
    # loop through i values
    for ii,i in enumerate(active):
        # calculate the nhat values
        nhat = decompose_j(i, Nvec)
        # loop through the j values
        for jj,j in enumerate(active):
            # calculate the n values
            n = decompose_j(j, Nvec)
            # calcualte Aij entry
            #####################
            
            if str(i) in con:
                if i == j:
                    A[ii,jj] = 1.
                else:
                    A[ii,jj] = 0.
            else:
                # initialize summ to 0
                summ = 0.
                # loop through points in dataset
                for p in range(1,k+1):
                    # initialize product series variable to 1
                    if y[p-1] != None:
                        prod = 1.
                    else:
                        prod = 0.
                    # loop through dimensions
                    for v in range(1,V+1):
                        # multiply the term onto the product series
                        prod *= x[p-1,v-1] ** (n[v-1] + nhat[v-1])
                    #==================================================
                    # add weighting factor
                    if callable(weighting):
                        prod *= weighting(x, y, p-1)
                    #==================================================
                    # add the product series variable to the summation
                    if percent:
                        if y[p-1] != None: summ += prod/abs(y[p-1])
                    else:
                        summ += prod
                # set Aij to the finalized summation
                A[ii,jj] = summ
        # calculate bi entry
        ####################
        if str(i) in con:
            b[ii] = con[str(i)]
        else:
            # initialize summation variable to 0
            summ = 0.
            # loop through points in the dataset
            for p in range(1,k+1):
                # initialize the product series variable to 1
                if y[p-1] != None:
                    prod = 1.
                else:
                    prod = 0.
                # loop through the dimensions
                for v in range(1,V+1):
                    # multiply the term onto the product series
                    prod *= x[p-1,v-1] ** nhat[v-1]
                #==================================================
                # add weighting factor
                if callable(weighting) and y[p-1] != None:
                    prod *= weighting(x, y, p-1)
                #==================================================
                # add the term onto the summation
                if percent:
                    summ += prod
                else:
                    if y[p-1] != None: summ += y[p-1] * prod
            # set bi to the finalized summation
            b[ii] = summ
        if verbose: prog.display()
    #solve Aa=b equation
    ########################################################################
    if verbose: print('solving the Aa=b equation')
    a = np.linalg.solve(A,b)
    #input the missing 0 coefficients into 'a' so that it can be used with the multidimensional_poly_func
    ########################################################################
    for i in range(J):
        if not i in active:
            a = np.insert(a,i,0.)
            active = np.insert(active,i,0)
    #calculate R^2 value
    ########################################################################
    r = multivariableR2(a, Nvec, x, y, verbose=verbose)       #r2_2D(a,M,x,y,z)
    #return values
    return a, r

def multivariablePolynomialFunction(a, Nvec, x):
    """
    Multivariable Polynomial Function
    
    inputs:
    
        a = list of the polynomial coefficients. 'a' has a length equal to
            the products of the Nvec elements plus one.
            i.e.: (Nvec0+1)*(Nvec1+1)*...
        Nvec = list with a length V equal to the number of independent
            variables. The ith element values are integers of polynomial
            order of the ith independent variable
        ??x = numpy matrix of size k by V, where k is the total number of
            points in the dataset. The ith column represents the ith
            independent variable values with the rows representing the
            different data points
    
    returns:
    
        f = value of the multivariable polynomial function for the given
            independent variables
    """
    # initialize summation to 0
    f = 0.
    # calculate total number of datapoints
    k = len(a)
    # calculate total number of dimensions
    V = len(x)
    # loop through the datapoints
    for j in range(k):
        # calculate the n values
        n = decompose_j(j, Nvec)
        # initialize the product series variable to 1
        prod = 1.
        # loop through the dimensions
        for v in range(V):
            # multiply onto the product series the term
            prod *= x[v] ** n[v]
        # add onto the summation the proper term
        f += a[j] * prod
    # return the finalized summation value
    return f

def multivariableR2(a, Nvec, xx, yy, verbose=True):
    """
    Routine to calculate the R^2 value of a multivariable polynomial fit to
    a dataset
    
    inputs:
    
        a = array of polynomial coefficients
        Nvec = list of integers representing the polynomial order of the
            independent variables
        xx = numpy matrix of size k by V, where k is the total number of
            points in the dataset. The ith column represents the ith
            independent variable values with the rows representing the
            different data points
        yy = list with a length of k with the dependent variable values
    
    returns:
    
        R2 = R^2 value, or the coefficient of determination
    """
    # ensure x and y are in the proper format
    x = np.copy(xx)
    y = np.copy(yy)
    ynew = np.copy([temp for temp in y if temp != None])
    # calculate k
    k = len(ynew)
    # calculate mean y value
    y_ = sum(ynew) / float(k)
    # calculate the SSt value
    SSt = sum( (ynew - y_) ** 2. )
    # initialize the f array
    f = []
    # loop through the datapoints
    if verbose: prog = oneLineProgress(len(y), msg='Determining R^2 of the fit')
    for i in range(len(y)):
        if y[i] == None:
            if verbose: prog.display()
            continue
        # calculate the f value from the polynomial function
        f.append( multivariablePolynomialFunction(a,Nvec,x[i,:]) )
        if verbose: prog.display()
    f = np.copy(f)
    # calculate the SSr term
    SSr = sum( (ynew - f) ** 2. )
    # calculate and return the R^2 value
    return 1. - SSr / SSt

def multivariableRMS(raw_x, raw_y, a, Nvec, verbose=True):
    """
    Routine to calculate the RMS and RMSN errors of a multivariable
    polynomial fit to a dataset
    
    inputs:
    
        raw_x = array of independent variables values from the dataset of
            shape (k,V) where k is the total number of points in the dataset
            and V is the number of independent variables
        raw_y = array of dependent variable values from the dataset
        a = array of polynomial coefficients
        Nvec = list of integers representing the polynomial order of the
            independent variables
    
    returns:
    
        RMS = root, mean, squared error
        RMSN = root, mean, squared, normalized error
    """
    x = np.copy( raw_x )
    y = np.copy( raw_y )
    avg = np.mean(abs(y))
    k = len(x[:,0])
    func = np.zeros(k)
    e = np.zeros(k)
    e_per = np.zeros(k)
    if verbose: prog = oneLineProgress(k, msg='Determining RMS of the fit')
    for i in range(k):
        func[i] = multivariablePolynomialFunction(a, Nvec, x[i])
        e[i] = (y[i] - func[i]) ** 2.
        e_per[i] = ((y[i] - func[i])/avg) ** 2.
        if verbose: prog.display()
    return np.sqrt(np.mean(e)), np.sqrt(np.mean(e_per))

def compose_j(n, Nvec):
    """
    Eq. 4 in Poly Fits Derivation. Routine to compose the j counter from
    the n values. Can also be used as Eq. 10 with the i counter and the nhat
    values.
    
    inputs:
    
        n = list of integer values representing the independent variables'
            exponents for the jth term in the multidimensional polynomial
            function (Eq. 2)
        Nvec = list of integers representing the polynomial order of the
            independent variables
    
    returns:
    
        j = integer representing the column of the A matrix or the jth
            polynomial coefficient
    """
    # calculate V
    V = len(Nvec)
    # initialize j to 0
    j = 0
    # loop through independent variables
    for v in range(1,V+1):
        # initialize product series to 1
        prod = 1
        # loop through w values for product series
        for w in range(v+1,V+1):
            # multiply on the term to the product series
            prod *= Nvec[w-1] + 1
        # add on term onto j
        j += n[v-1] * prod
    return j

def decompose_j(j, Nvec):
    """
    Eq. 5 in Poly Fits Derivation. Routine to decompose the j counter into
    the n values. Can also be used as Eq. 1 with the i counter and the nhat
    values.
    
    inputs:
    
        j = integer representing the column of the A matrix or the jth
            polynomial coefficient
        Nvec = list of integers representing the polynomial order of the
            independent variables
    
    returns:
    
        n = list of integer values representing the independent variables'
            exponents for the jth term in the multidimensional polynomial
            function (Eq. 2)
    """
    # calculate V
    V = len(Nvec)
    # initialize n values to nothing
    n = [[]]*V
    # loop through the n values that need to be solved, starting at the highest and working down
    for v in range(V,0,-1):
        # initialize the denomenator product series to 1
        denom = 1
        # loop through the w values needed for the product series
        for w in range(v+1,V+1):
            # multiply on the terms for the denomenator product series
            denom *= Nvec[w-1] + 1
        # initialize the summation variable to 0
        summ = 0
        # loop through the u values necessary for the summation
        for u in range(v+1,V+1):
            # initialize the product series variable inside the summation to 1
            prod = 1
            # loop through the s values needed for the product series that is inside the summation
            for s in range(u+1,V+1):
                # multiply on the term for the product series that is inside of the summation
                prod *= Nvec[s-1] + 1
            # add on the needed term to the summation series
            summ += n[u-1] * prod
        # finally calculate the n value cooresponding to this v
        n[v-1] = int(round( ((j-summ)/denom)%(Nvec[v-1]+1) ))
    return n

############################################################################
############################################################################
############################################################################

def calcJ(Nvec):
    J = 1
    for n in Nvec:
        J *= n + 1
    return J

def kDecompose(k, V):
    t = 1
    ###################################################
    ## find the category
    c = 0
    vals = [0] * V
    while k > t:
        c += 1
        m = [c] * (V-1)
        vals[0] = calcJ(m)
        for j in range(V-1):
            m[j] -= 1
            vals[j+1] = calcJ(m)
        t += sum(vals)
    if c == 0:
        return [0]*V
    ####################################################
    ## find the subcategory
    for sc in range(V-1,-1,-1):
        t -= vals[sc]
        if k > t:
            break
    ####################################################
    ## determine n
    # initialize n
    n = [None]*V
    n[sc] = c
    # create mx based on the sc, then decompose to get m
    mx = [c]*(V-1)
    for i in range(sc):
        mx[i] -= 1
    m = decompose_j(k-t-1, mx)
    # set m values into n and return
    j = -1
    for i in range(V):
        if i != sc:
            j += 1
            n[i] = m[j]
    return n

def kCompose(n):
    ########################################
    V = len(n)
    if V == 1: return calcJ(n)
    mx = max(n)
    if mx == 0: return 1
    k = 1
    ## calculate lower number sets
    for i in range(1,mx):
        m = [i] * (V-1)
        k += calcJ(m)
        for j in range(V-1):
            m[j] -= 1
            k += calcJ(m)
    ## calculate location in current number set
    for i in range(V):
        M = [mx]*(V-1)
        for j in range(i):
            M[j] -= 1
        if n[i] != mx:
            k += calcJ(M)
        else:
            m = [n[j] for j in range(V) if j != i]
            k += compose_j(m, M) + 1
            return k
    raise ValueError('Unable to compose n into k: current k value {}'.format(k))

class oneLineProgress():
    
    def __init__(self, Total, msg='', showETR=True):
        self.total = Total
        self.msg = msg
        self.count = 0
        self.showETR = showETR
        self.start = dt.now()
        self.rollTimer = dt.now()
        self.rollCount = -1
        self.rollDelta = 0.2
        self.display()
    
    def increment(self):
        self.count += 1
    
    def decrement(self):
        self.count -= 1
    
    def __str__(self):
        pass
    
    def __len__(self):
        l = len(str(self))
        self.decrement()
        return l
    
    def Set(self, count):
        self.count = count
    
    def display(self):
        rolling = '-\\|/'
        rollDelta = (dt.now()-self.rollTimer).total_seconds()
        
        p2s = False
        if rollDelta >= self.rollDelta or self.rollCount == -1:
            p2s = True
            self.rollTimer = dt.now()
            self.rollCount += 1
            if self.rollCount >= len(rolling):
                self.rollCount = 0
        
        perc = self.count / self.total * 100.
        self.increment()
        
        if not p2s and perc < 100.: return
        
        s = '\r' + ' '*(len(self.msg)+50) + '\r'
        s += self.msg + ' '*4
        
        # j = 0
        for i in range(10):
            if perc >= i*10:
                j = i
        
        if perc < 100.:
            s += u'\u039e'*j + rolling[self.rollCount] + '-'*(9-j)
        else:
            s += u'\u039e'*10
        
        # for i in range(1,11):
            # if i*10 <= perc:
                # s += u'\u039e'
            # else:
                # s += '-'
        s += ' '*4 + '{:7.3f}%'.format(perc)
        if not self.showETR:
            if perc >= 100.: s += '\n'
            print(s, end='')
            return
        
        if perc <= 0:
            etr = '-:--:--.------'
            s += ' '*4 + 'ETR = {}'.format(etr)
        elif perc >= 100.:
            s += ' '*4 + 'Run Time {}'.format(dt.now()-self.start) + '\n'
        else:
            time = (dt.now()-self.start).total_seconds()
            etr = td(seconds=time / perc * 100. - time)
            s += ' '*4 + 'ETR = {}'.format(etr)
        print(s, end='')
        return

def zSort(v, *W, ascend=True, verbose=True, msg='Sorting the arrays'):
    k = len(v)
    for w in W:
        if len(w) != k: raise ValueError('All arrays need to be the same length in zSort')
    c = []
    if verbose: prog = oneLineProgress(sum([i for i in range(k)])+len(W), msg=msg)
    for m in range(k):
        for j in range(k-1,m,-1):
            i = j-1
            
            if (ascend and v[j] < v[i]) or (not ascend and v[j] > v[i]):
                c.append(j)
                temp = v[j]
                v[j] = v[i]
                v[i] = temp
            if verbose: prog.display()
    
    for w in W:
        for j in c:
            i = j-1
            temp = w[j]
            w[j] = w[i]
            w[i] = temp
        if verbose: prog.display()

def isClose(x, y, tol=1.e-12):
    return y-tol <= x and x <= y+tol


def autoPolyFit(X, y, max_order=6, tol=1.e-12, sigma=None, sigma_multiplier=1., verbose=True):
    '''
    autoPolyFit function performs a mutivariable polynomial curve
    fit to a dataset and automatically determines which terms in the
    polynomial to use based on a balance between a goodness of the fit and a
    predictive capabilities measure that attempts to make the model compact.
    
    inputs:
        X : numpy array of shape (N,m). X consists of all the independent
            variables in the dataset. N is the number of data points in the
            set and m is the number of independent variables
        y : list or numpy array with length N. y is the dependent variable
            values cooresponding to the independent variables in X
        max_order : optional integer. gives the max order of polynomial for
            any one of the independent varialbes to try. defaults to 12
        tol : optional float. Gives the cut-off value for any polynomial
            coefficient to not be included in the final results. If a
            coefficient has an absolute value below tol, it won't be
            included. defaults to 1e-12
        sigma : optional float. value used to determine the trade off
            between how good of a fit to perform and how many terms to keep.
            defaults to None, which causes the function to calculate sigma
            automatically using the mean squared of the difference of the
            independent variable values with respect to the mean independent
            variable value of the dataset
        sigma_multiplier : optional float. term multiplied onto sigma to
            change it's value. Allows using a multiple of the automatically
            determined sigma value. Defaults to 1.
    
    returns:
        list : a list of the polynomial coefficients in a form used by the
            poly_fits module
        list : a list of the max polynomial orders for each independent
            variable. The length of this list is therefore m. This list is
            the 'Nvec' object used in the poly_fits module
        float : the coefficient of determination, R^2 value, representing
            the goodness of the fit
    '''
    ## number of independent variables
    m = X.shape[1]
    ## max range of polynomials to try
    Nvec = tuple([max_order]*m)
    ## number of datapoints
    N = len(y)
    ## number of p functions
    K = kCompose(Nvec)
    ###################################################################
    ##           determine the orthogonal p functions
    if verbose: prog = oneLineProgress(K-1, msg='Determining the orthogonal p functions')
    ## initialize the P matrix
    P = np.zeros((N, K))
    P[:,0] = 1.
    ## loop thru k values
    for k in range(2,K+1):
        ## determine the set associated with k
        n = decompose_j(k-1, Nvec)
        ## find pkhat and mu
        mu = None
        for i in range(m):
            nhat = n[:]
            if nhat[i] > 0:
                nhat[i] -= 1
                khat = compose_j(nhat, Nvec) + 1
                if khat < k:
                    mu = i
                    break
        if mu == None: raise ValueError('Unable to find previous khat set')
        pkhat = P[:,khat-1]
        xmu = X[:,mu]
        ## calculate pk first term
        temp = xmu * pkhat
        phik = sum(n)
        pk = temp[:]
        ## loop thru summation in eq 18
        for j in range(1,k):
            ## check if value is needed
            phij = sum(decompose_j(j-1, Nvec))
            if phik - phij <= 2:
                ## calculate gamma
                pj = P[:,j-1]
                gamma = np.dot(pj, temp) / np.dot(pj, pj)
                pk -= gamma * pj
        ## add pk to P
        P[:,k-1] = pk[:]
        if verbose: prog.display()
    #################################################################
    ##              sort the p functions by effectiveness
    order = [i for i in range(K)]
    ranks = [None] * K
    for i in range(K):
        pj = P[:,i]
        pjdot = np.dot(pj, pj)
        ajhat = np.dot(pj, y) / pjdot
        ranks[i] = ajhat ** 2. * pjdot
    zSort(ranks, order, ascend=False, msg='Sorting the p functions by effectivenss', verbose=verbose)
    Pordered = np.zeros((N,K))
    for i,o in enumerate(order):
        Pordered[:,i] = P[:,o]
    P = Pordered[:,:]
    ###################################################################
    ##          determine how many of the orthogonal p functions to use
    if verbose: prog = oneLineProgress(K, msg='Determining number of p functions to use')
    PSEold = None
    foundMin = False
    if sigma == None:
        yavg = sum(y) / N
        sigma = sum([(i - yavg)**2. for i in y]) / N
    sigma *= sigma_multiplier
    
    for n in range(1,K+1):
        Phat = P[:,:n]
        ahat = np.matmul(np.matmul(np.linalg.inv(np.matmul(Phat.transpose(), Phat)), Phat.transpose()), y)
        yhat = np.matmul(Phat, ahat)
        MSE = np.dot(y - yhat, y - yhat) / N
        PSEnew = MSE + sigma * n / N
        if verbose: prog.display()
        if PSEold == None or PSEnew <= PSEold:
            PSEold = PSEnew
        else:
            foundMin = True
            nn = n-1
            P = Phat[:,:nn]
            order = order[:nn]
            if verbose: 
                prog.Set(K)
                prog.display()
            break
    if not foundMin:
        raise ValueError('Unable to find minimum PSE')
    ###################################################################
    ##              final coefficients and polynomial size
    if verbose: prog = oneLineProgress(4+nn, msg='Determining final coefficients and polynomial size')
    b = np.zeros((nn,nn))
    for k in range(1,nn+1):
        j = k - 1
        pj = P[:,j]
        w = np.ones((N,k))
        for i in range(k):
            n = decompose_j(order[i], Nvec)
            for ii,e in enumerate(n):
                w[:,i] *= X[:,ii] ** e
        vals = np.matmul(np.matmul(np.linalg.inv(np.matmul(w.transpose(), w)), w.transpose()), pj)
        b[j,:k] = vals[:]
        if verbose: prog.display()
    
    A = [np.dot(P[:,i],y)/np.dot(P[:,i],P[:,i]) for i in range(nn)]
    if verbose: prog.display()
    
    c = [np.dot(A,b[:,i]) for i in range(nn)]
    js = [decompose_j(order[i], Nvec) for i in range(nn)]
    if verbose: prog.display()
    
    js = [js[i] for i in range(nn) if not isClose(c[i], 0., tol=tol)]
    c = [i for i in c if not isClose(i, 0., tol=tol)]
    if verbose: prog.display()
    
    nvec = [None]*m
    for i in range(m):
        nvec[i] = max([j[i] for j in js])
    JJ = 1
    for n in nvec:
        JJ *= n+1
    
    a = [0.] * JJ
    for j in range(JJ):
        n = decompose_j(j, nvec)
        for i in range(len(c)):
            if n == js[i]:
                a[j] = c[i]
    if verbose: prog.display()
    
    return a, nvec, multivariableR2(a, nvec, X, y, verbose=verbose)

