# converted from MLfit_GLM.m
# [TODO] add packages to use

def MLfit_GLM(gg, Stim, optimArgs=None):
    """
    #  [gg,neglogli,H,Xstruct] = MLfit_GLM(gg,Stim,optimArgs)
    #
    #  Computes the ML estimate for GLM params, using grad and hessians.
    #  Assumes basis for temporal dimensions of stim filter

    #  Inputs:
    #     gg = param struct
    #     Stim = stimulus
    #     optimArgs = cell array of optimization params (optional)

    #  Outputs:
    #     ggnew = new param struct (with ML params);
    #  neglogli = negative log-likelihood at ML estimate
    #         H = Hessian of negative log-likelihood at ML estimate
    #   Xstruct = structure with design matrices for spike-hist and stim terms
    """

    # Set optimization parameters
    algopts=getFminOptsForVersion(version)
    if nargin > 2:
        opts=optimset(algopts[:],optimArgs[:])
    else:
        opts=optimset(algopts[:])

    # --- Create design matrix extract initial params from gg ----------------
    prs0,Xstruct=setupfitting_GLM(gg,Stim,nargout=2)

    # --- Set loss function --------------------------------------------------
    #if isequal(Xstruct.nlfun,expfun) or isequal(Xstruct.nlfun,exp):
    lfunc=lambda prs=None: Loss_GLM_logli_exp(prs,Xstruct)
    #else:
    #    lfunc=lambda prs=None: Loss_GLM_logli(prs,Xstruct)

    # --- minimize negative log likelihood --------------------
    prsML,neglogli=fminunc(lfunc,prs0,opts)

    # Compute Hessian if desired
    # if nargout > 2:
    #     neglogli,__,H=Loss_GLM_logli(prsML,Xstruct)

    # Put returned vals back into param structure ------
    gg=reinsertFitPrs_GLM(gg,prsML,Xstruct)
    # #----------------------------------------------------

# Optional debugging code
# #----------------------------------------------------

    # # ------ Check analytic gradients and Hessians -------
#  HessCheck(lfunc,prs0,opts);
#  HessCheck_Elts(@Loss_GLM_logli, [1 12],prs0,opts);
#  tic; [lival,J,H]=lfunc(prs0); toc;
