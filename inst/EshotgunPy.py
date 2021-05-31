import os
import numpy as np
import GPy as gp
import pygmo as pg
import scipy
import cma
import nlopt
import warnings
from pyDOE2.doe_lhs import lhs

# The following code, besides the function "def callShotgun(Xtr, Ytr, f_lb, f_ub, q=10, epsilon = 0.1):", was developed by
# George De Ath, Richard M. Everson, Jonathan E. Fieldsend, and Alma A. M. Rahat. 2020.
# e-shotgun : e-greedy Batch Bayesian Optimisation. In Genetic and Evolutionary Computation Conference (GECCO ’20), July 8–12, 2020, Cancún, Mexico. 
# ACM, New York, NY, USA, 9 pages. 
# https://doi.org/10.1145/3377930.3390154
# https://github.com/georgedeath/eshotgun

def aquisition_DIRECT(model, aq_func, cf, aq_kwargs={}):
    dim = model.X.shape[1]
    cf_info = CF_INFO(cf)

    def f(X):
        # if there's a constraint function, evaluate it
        if cf_info.not_valid(X):
            return cf_info.bad_value

        X = X.reshape(-1, dim)

        pred_mean, pred_var = model.predict(X, full_cov=False)
        pred_std = np.sqrt(pred_var)

        # negative as optimisers are minimising
        return -np.squeeze(aq_func(pred_mean, pred_std, **aq_kwargs))

    return f


def aquisition_LBFGSB(model, aq_func, cf, aq_kwargs={}):
    dim = model.X.shape[1]
    cf_info = CF_INFO(cf)

    def f(X):
        # if there's a constraint function, evaluate it
        if cf_info.not_valid(X):
            return cf_info.bad_value

        X = X.reshape(-1, dim)

        pred_mean, pred_var = model.predict(X, full_cov=False)
        pred_std = np.sqrt(pred_var)

        # negative as optimisers are minimising
        return -np.squeeze(aq_func(pred_mean, pred_std, **aq_kwargs))

    return f


def aquisition_CMAES(model, aq_func, cf=None, aq_kwargs={}):
    dim = model.X.shape[1]
    cf_info = CF_INFO(cf)

    def f(X):
        # if X is a list, assume that each list element is either
        # a float or a numpy array -> convert to (N, ndim) numpy array
        if isinstance(X, list):
            X = np.reshape(np.array(X), (len(X), dim))

        # else must be a numpy array so we have to assume that it is (N, ndim)
        elif len(X.shape) != 2:
            X = np.atleast_2d(X)

        pred_mean, pred_var = model.predict(X, full_cov=False)
        pred_std = np.sqrt(pred_var)

        # negative as optimisers are minimising
        aq_res = -aq_func(pred_mean, pred_std, **aq_kwargs)
        aq_res = aq_res.ravel().tolist()

        # evaluate constraint function for each decision vector
        for i in range(len(aq_res)):
            if cf_info.not_valid(X[i, :]):
                aq_res[i] = cf_info.bad_value.flat[0]

        if len(aq_res) == 1:
            return aq_res[0]

        return aq_res

    return f


class CF_INFO:
    def __init__(self, cf):
        self.cf = cf
        self.got_cf = cf is not None
        self.bad_value = np.array(np.inf)

    def not_valid(self, X):
        return self.got_cf and (not self.cf(X))


def minimise_with_DIRECT(f, lb, ub, maxeval=5000, cf=None, ftol_abs=1e-15):
    dim = lb.size

    # define a direct optimisation instance
    opt = nlopt.opt(nlopt.GN_DIRECT_L_RAND, dim)
    opt.set_min_objective(f)

    # set the lower and upper bounds
    opt.set_lower_bounds(lb)
    opt.set_upper_bounds(ub)

    # set max evaluations and function tolerance - this latter option
    # has a big performance influence when optimising in a small region
    opt.set_maxeval(maxeval)
    opt.set_ftol_abs(ftol_abs)

    # perform the optimisation
    xopt = opt.optimize(np.random.uniform(lb, ub))

    return xopt


def minimise_with_CMAES(f, lb, ub, maxeval=5000, cf=None, ftol_abs=1e-15):
    # set the options
    cma_options = {'bounds': [list(lb), list(ub)],
                   'tolfun': ftol_abs,
                   'maxfevals': maxeval,
                   'verb_disp': 0,
                   'verb_log': 0,
                   'verbose': -1,
                   'CMA_stds': np.abs(ub - lb),
                   }

    if cf is None:
        x0 = lambda: np.random.uniform(lb, ub)

    else:
        def inital_point_generator(cf, lb, ub):
            def wrapper():
                while True:
                    x = np.random.uniform(lb, ub)
                    if np.all(x >= lb) and np.all(x <= ub) and cf(x):
                        return x
            return wrapper

        class feas_func:
            def __init__(self, cf):
                self.cf = cf
                self.c = 0

            def __call__(self, x, f):
                if self.c > 10000:
                    return True

                is_feas = self.cf(x)

                if not is_feas:
                    self.c += 1

                return is_feas

        is_feasible = feas_func(cf)
        cma_options['is_feasible'] = is_feasible
        x0 = inital_point_generator(cf, lb, ub)

    # ignore warnings about flat fitness (i.e starting in a flat EI location)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')

        # run CMA-ES with bipop (small and large population sizes) for up to
        # 9 restarts (or until it runs out of budget)
        xopt, _ = cma.fmin2(f, x0=x0, sigma0=0.25, options=cma_options,
                            bipop=True, restarts=9)
        warnings.resetwarnings()

    return xopt


def minimise_with_LBFGSB(f, lb, ub, maxeval=5000, cf=None, ftol_abs=1e-15):
    dim = lb.size

    # number of optimisation runs and *estimated* number of L-BFGS-B
    # function evaluations per run; note this was calculate empirically and
    # may not be true for all functions.
    N_opt_runs = 10
    fevals_assumed_per_run = 100

    N_LHS_samples = maxeval-(N_opt_runs*fevals_assumed_per_run)
    if N_LHS_samples <= N_opt_runs:
        N_LHS_samples = N_opt_runs

    # initially perform a grid search using LHS (maximin) for N_LHS_samples
    x0_points = lhs(dim, samples=N_LHS_samples, criterion='m')
    x0_points = x0_points * (ub - lb)[np.newaxis, :] + lb[np.newaxis, :]

    fx = f(x0_points).ravel()

    # select the top N_opt_runs to evaluate with L-BFGS-B
    x0_points = x0_points[np.argsort(fx)[:N_opt_runs], :]

    # Find the best optimum by starting from n_restart different random points.
    bounds = [(l, b) for (l, b) in zip(lb, ub)]

    # storage for the best found location (xb) and its function value (fx)
    xb = np.zeros((N_opt_runs, dim))
    fx = np.zeros((N_opt_runs, 1))

    # ensure we're using a good stopping criterion
    # ftol = factr * numpy.finfo(float).eps
    factr = ftol_abs / np.finfo(float).eps

    for i, x0 in enumerate(x0_points):
        xb[i, :], fx[i, :], d = scipy.optimize.fmin_l_bfgs_b(f,
                                                             x0=x0,
                                                             bounds=bounds,
                                                             approx_grad=True,
                                                             factr=factr)

    best_idx = np.argmin(fx.flat)
    return xb[best_idx, :]


def minimise_with_LBFGSB_once(f, lb, ub, x0, maxeval=5000, ftol_abs=1e-15):
    bounds = [(l, b) for (l, b) in zip(lb, ub)]

    # ensure we're using a good stopping criterion
    # ftol = factr * numpy.finfo(float).eps
    factr = ftol_abs / np.finfo(float).eps

    xb, fx, d = scipy.optimize.fmin_l_bfgs_b(f, x0=x0, bounds=bounds,
                                             approx_grad=True, factr=factr)

    return xb


def NSGA2_pygmo(model, fevals, lb, ub, cf=None):
    """Finds the estimated Pareto front of a GPy model using NSGA2 [1]_.

    Parameters
    ----------
    model : GPy.models.gp_regression.GPRegression
        GPy regression model on which to find the Pareto front of its mean
        prediction and standard deviation.
    fevals : int
        Maximum number of times to evaluate a location using the model.
    lb : (D, ) numpy.ndarray
        Lower bound box constraint on D
    ub : (D, ) numpy.ndarray
        Upper bound box constraint on D
    cf : callable, optional
        Constraint function that returns True if it is called with a
        valid decision vector, else False.

    Returns
    -------
    X_front : (F, D) numpy.ndarray
        The F D-dimensional locations on the estimated Pareto front.
    musigma_front : (F, 2) numpy.ndarray
        The corresponding mean response and standard deviation of the locations
        on the front such that a point X_front[i, :] has a mean prediction
        musigma_front[i, 0] and standard deviation musigma_front[i, 1].

    Notes
    -----
    NSGA2 [1]_ discards locations on the pareto front if the size of the front
    is greater than that of the population size. We counteract this by storing
    every location and its corresponding mean and standard deviation and
    calculate the Pareto front from this - thereby making the most of every
    GP model evaluation.

    References
    ----------
    .. [1] Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan.
       A fast and elitist multiobjective genetic algorithm: NSGA-II.
       IEEE Transactions on Evolutionary Computation 6, 2 (2001), 182–197.
    """
    # internal class for the pygmo optimiser
    class GPY_WRAPPER(object):
        def __init__(self, model, lb, ub, cf, evals):
            # model = GPy model
            # lb = np.array of lower bounds on X
            # ub = np.array of upper bounds on X
            # cf = callable constraint function
            # evals = total evaluations to be carried out
            self.model = model
            self.lb = lb
            self.ub = ub
            self.nd = lb.size
            self.got_cf = cf is not None
            self.cf = cf
            self.i = 0  # evaluation pointer

        def get_bounds(self):
            return (self.lb, self.ub)

        def get_nobj(self):
            return 2

        def fitness(self, X):
            X = np.atleast_2d(X)
            f = model_fitness(X, self.model, self.cf, self.got_cf,
                              self.i, self.i + X.shape[0])
            self.i += X.shape[0]
            return f

    # fitness function for the optimiser
    def model_fitness(X, model, cf, got_cf, start_slice, end_slice):
        valid = True

        # if we select a location that violates the constraint,
        # ensure it cannot dominate anything by having its fitness values
        # maximally bad (i.e. set to infinity)
        if got_cf:
            if not cf(X):
                f = [np.inf, np.inf]
                valid = False

        if valid:
            mu, sigmaSQR = model.predict(X, full_cov=False)
            # note the negative sigmaSQR here as NSGA2 is minimising
            # so we want to minimise the negative variance
            f = [mu.flat[0], -np.sqrt(sigmaSQR).flat[0]]

        # store every point ever evaluated
        model_fitness.X[start_slice:end_slice, :] = X
        model_fitness.Y[start_slice:end_slice, :] = f

        return f

    # get the problem dimensionality
    D = lb.size

    # NSGA-II settings
    POPSIZE = D * 100
    N_GENS = int(np.ceil(fevals / POPSIZE))
    TOTAL_EVALUATIONS = POPSIZE * N_GENS

    nsga2 = pg.algorithm(pg.nsga2(gen=1,
                                  cr=0.8,       # cross-over probability.
                                  eta_c=20.0,   # distribution index (cr)
                                  m=1 / D,        # mutation rate
                                  eta_m=20.0))  # distribution index (m)

    # preallocate the storage of every location and fitness to be evaluated
    model_fitness.X = np.zeros((TOTAL_EVALUATIONS, D))
    model_fitness.Y = np.zeros((TOTAL_EVALUATIONS, 2))

    # problem instance
    gpy_problem = GPY_WRAPPER(model, lb, ub, cf, TOTAL_EVALUATIONS)
    problem = pg.problem(gpy_problem)

    # initialise the population
    population = pg.population(problem, size=POPSIZE)

    # evolve the population
    for i in range(N_GENS):
        population = nsga2.evolve(population)

    # indices non-dominated points across the entire NSGA-II run
    front_inds = pg.non_dominated_front_2d(model_fitness.Y)

    X_front = model_fitness.X[front_inds, :]
    musigma_front = model_fitness.Y[front_inds, :]

    # convert the standard deviations back to positive values; nsga2 minimises
    # the negative standard deviation (i.e. maximises the standard deviation)
    musigma_front[:, 1] *= -1

    return X_front, musigma_front

def estimate_L(model, xj, lengthscale, lb, ub):
    """
    Estimate the Lipschitz constant of f by taking maximizing the
    norm of the expectation of the gradient of *f*.

    Adapated from GPyOpt:
        GPyOpt/core/evaluators/batch_local_penalization.py

    """
    def df(x, model):
        x = np.atleast_2d(x)
        dmdx, _ = model.predictive_gradients(x)

        # simply take the norm of the expectation of the gradient
        res = np.sqrt((dmdx * dmdx).sum(1))
        
        # bfgs (scipy 1.5.0) expects shape (d,) rather than (1,d)
        if x.shape[0] == 1:
            res = res[0]
            
        return -res

    # generate bounds, box constraint on xj
    n_dim = xj.size

    # centred on xj, lengthscale wide in each dimension, subject to the domain
    df_lb = np.maximum(xj - lengthscale, lb)
    df_ub = np.minimum(xj + lengthscale, ub)

    # scipy bounds
    bounds = list(zip(df_lb, df_ub))

    # generate some samples in the box around xj and evaluate their gradient
    samples = np.random.uniform(df_lb, df_ub, size=(500, n_dim))
    samples = np.vstack([samples, model.X])

    samples_df = df(samples, model)

    # select the starting point as that with the largest (negative) gradient
    x0 = samples[np.argmin(samples_df)]

    xopt, minusL, _ = scipy.optimize.fmin_l_bfgs_b(df,
                                                   x0,
                                                   bounds=bounds,
                                                   args=(model,),
                                                   maxiter=2000,
                                                   approx_grad=True)

    L = -np.squeeze(minusL).item()

    if L < 1e-7:
        L = 10  # to avoid problems in cases in which the model is flat.

    return L


def calculate_ball_radius(xj, model, lb, ub):
    # r_j \leq ||\mu(x_j) - M|| / L   +  \gamma * \sigma(x_j) / L
    # where L = estimated Lipschitz constant locally
    #           within a lengthscale of x_j
    # gamma = 1
    # M = min Ytr

    # calculates the radius of the ball to sample in, centred on xj
    ls = model.kern.lengthscale[0]

    # locally estimate the Lipshitz constant as the largest gradient within
    # a lengthscale of xj and within the problem domain
    L = estimate_L(model, xj, ls, lb, ub)

    # estimate the best function value to be the best seen function value
    M = np.min(model.Y)

    # gamma: Asynchronous Batch Bayesian Optimisation
    #        with Improved Local Penalisation
    gamma = 1

    mu_xj, sigmaSQR_xj = model.predict(np.atleast_2d(xj))
    sigma_xj = np.sqrt(sigmaSQR_xj)

    rj = (np.abs(mu_xj - M) + (gamma * sigma_xj)) / L

    return np.squeeze(rj)


def egreedy_shotgun(model, f_lb, f_ub, feval_budget, q, cf,
                    epsilon, pf=False, aq_func=lambda mu, sigma: -mu):

    n_dim = f_lb.size

    # epsilon of the time, randomly choose a point in space
    if np.random.uniform() < epsilon:
        if pf:
            # calculate the pareto front and randomly select a location on it
            X_front, musigma_front = NSGA2_pygmo(model, feval_budget,
                                                 f_lb, f_ub, cf)
            xj = X_front[np.random.choice(X_front.shape[0]), :]

        else:
            xj = np.random.uniform(f_lb, f_ub)

    # else find the point that maximises an acquisition function
    else:
        # we can only use CMA-ES on 2 or more dimensional functions
        if n_dim > 1:
            opt_acq_func, opt_caller = aquisition_CMAES, minimise_with_CMAES

        # else use L-BFGS-B
        else:
            opt_acq_func, opt_caller = aquisition_LBFGSB, minimise_with_LBFGSB

        # cma-es wrapper for it (this just negates it)
        f = opt_acq_func(model, aq_func, cf)

        # minimise the negative acquisition function with cma-es
        xj = opt_caller(f, f_lb, f_ub, feval_budget, cf=cf)

    # radius of the ball around xj in which to sample new locations
    rj = calculate_ball_radius(xj, model, f_lb, f_ub)

    # maximum ball radius = half the size of norm of the domain
    rj = np.minimum(rj, np.linalg.norm(f_ub - f_lb) / 2)

    Xnew = [xj]

    # sample new locations, x_i ~ N(xj, rj), where N is the normal distribution
    while len(Xnew) < q:
        Xtest = np.random.normal(loc=xj, scale=rj)

        # ensuring batch location lies within the problem domain
        while True:
            # ndarray of indices of locations out of the domain
            bad_inds = np.flatnonzero(np.logical_or(Xtest < f_lb,
                                                    Xtest > f_ub))

            # if they're all in bounds, break the checking loop
            if bad_inds.size == 0:
                break

            # isotropic scaling so we can individually sample components
            if bad_inds.size > 0:
                Xtest[bad_inds] = np.random.normal(loc=xj[bad_inds], scale=rj)

        Xnew.append(Xtest)

    Xnew = np.array(Xnew)

    return Xnew


def egreedy_shotgun_v2(model, f_lb, f_ub, feval_budget, q, cf,
                       epsilon, pf=False, aq_func=lambda mu, sigma: -mu):

    n_dim = f_lb.size

    # always select the (estimated) best point in the acquisition function
    # we can only use CMA-ES on 2 or more dimensional functions
    if n_dim > 1:
        opt_acq_func, opt_caller = aquisition_CMAES, minimise_with_CMAES

    # else use L-BFGS-B
    else:
        opt_acq_func, opt_caller = aquisition_LBFGSB, minimise_with_LBFGSB

    # CMA-ES wrapper for it (this just negates it)
    f = opt_acq_func(model, aq_func, cf)

    # minimise the negative acquisition function with CMA-ES
    xj = opt_caller(f, f_lb, f_ub, feval_budget, cf=cf)

    Xnew = [xj]

    # decide how many random samples to place
    r = np.random.uniform(0, 1, size=q - 1)
    n_random = np.count_nonzero(r < epsilon)

    # generate the remaining points randomly
    if n_random > 0:
        if pf:
            # calculate the Pareto front
            X_front, _ = NSGA2_pygmo(model, feval_budget, f_lb, f_ub, cf)

            # check the front is large enough
            n_front_samples = np.minimum(n_random, X_front.shape[0])

            # if not then uniformly sample the extra points needed
            for _ in range(n_random - n_front_samples):
                Xnew.append(np.random.uniform(f_lb, f_ub))

            # randomly select locations on the front to evaluate
            inds = np.random.choice(X_front.shape[0], size=n_front_samples,
                                    replace=False)
            for ind in inds:
                Xnew.append(X_front[ind, :])

        else:
            for _ in range(n_random):
                Xnew.append(np.random.uniform(f_lb, f_ub))

    # sample the remaining points via the shotgun blast approach
    if len(Xnew) < q:
        # radius of the ball around xj in which to sample new locations
        rj = calculate_ball_radius(xj, model, f_lb, f_ub)

        # maximum ball radius = half the size of norm of the domain
        rj = np.minimum(rj, np.linalg.norm(f_lb - f_ub) / 2)

        # sample new locations, x_i ~ N(xj, rj), where N is the normal
        # distribution
        while len(Xnew) < q:
            Xtest = np.random.normal(loc=xj, scale=rj)

            # ensuring that they lie within the problem domain
            if np.logical_and(np.all(Xtest >= f_lb), np.all(Xtest <= f_ub)):
                Xnew.append(Xtest)

    Xnew = np.array(Xnew)

    return Xnew


# force numpy, blas, and mkl to use a max number of threads
n_threads = 4
os.environ["MKL_NUM_THREADS"] = "{:d}" .format(n_threads)
os.environ["NUMEXPR_NUM_THREADS"] = "{:d}" .format(n_threads)
os.environ["OMP_NUM_THREADS"] = "{:d}" .format(n_threads)
os.environ["OPENBLAS_NUM_THREADS"] = "{:d}" .format(n_threads)
os.environ["VECLIB_MAXIMUM_THREADS"] = "{:d}" .format(n_threads)

def build_and_fit_GP(Xtr, Ytr):
    # create a gp model with the training data and fit it
    kernel = gp.kern.Matern52(input_dim=Xtr.shape[1], ARD=False)
    model = gp.models.GPRegression(Xtr, Ytr, kernel, normalizer=True)
    
    model.constrain_positive('')
    (kern_variance, kern_lengthscale,
     gaussian_noise) = model.parameter_names()
     
    model[kern_variance].constrain_bounded(1e-6, 1e6, warning=False)
    model[kern_lengthscale].constrain_bounded(1e-6, 1e6, warning=False)
    model[gaussian_noise].constrain_fixed(1e-6, warning=False)
    model.optimize_restarts(optimizer='lbfgs',
                            num_restarts=10,
                            num_processes=1,
                            verbose=False)
    return model
  

def callShotgun(Xtr, Ytr, f_lb, f_ub, q=10, epsilon = 0.1):
    model = build_and_fit_GP(Xtr, Ytr)
    
    if type(f_lb) is float:
      Xnew = egreedy_shotgun_v2(model, np.array([f_lb]), np.array([f_ub]), 1000, q, None, epsilon)
    else:
      Xnew = egreedy_shotgun_v2(model, f_lb, f_ub, 1000, q, None, epsilon)
    return Xnew
    
