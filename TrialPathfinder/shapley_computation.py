from .utils import *
    
def shapley_computation(cohort, features, drug_treatment, drug_control, name_rules, tolerance=0.001, iter_max=1000, covariates_cont=[], covariates_cat=[], thresh_censor=None, name_DrugName='DrugName', name_StartDate='StartDate', name_OutcomeDate='OutcomeDate', name_LastVisitDate='LastVisitDate', indicator_miss='Missing', random_seed=1001, verbose=0):
    '''
    Emulate trial
    Return HR, confidence interval, data fit for the cox model
    '''  
    np.random.seed(random_seed)
    
    def get_HR(name_rules, data_survival):
        '''Return HR given criteria and survival data'''
        data_cox = generate_trial_cox(cohort, data_survival, drug_treatment, drug_control, 
                                                  name_rules=name_rules, name_DrugName=name_DrugName,
                                                  covariates_cont=covariates_cont, covariates_cat=covariates_cat)
        HR, _ = cox(data_cox)
        return HR
    
    def compute_SEM(dHRs):
        '''Compute the standard error of the Monte Carlo mean'''
        dHRs = np.array(dHRs)
        SEM = np.mean([np.std(dHRs[:, i])/np.sqrt(dHRs.shape[0]) for i in range(dHRs.shape[1])])
        return SEM
    
    # Generate survival information.
    data_survival = generate_survival_data(features.copy(), cohort.name_PatientID,
                                           covariates=covariates_cont+covariates_cat, thresh_censor=thresh_censor,
                                           name_DrugName=name_DrugName, name_StartDate=name_StartDate, 
                                           name_OutcomeDate=name_OutcomeDate, name_LastVisitDate=name_LastVisitDate,
                                           indicator_miss=indicator_miss) 

    # HR for Empty set and full set
    HR_empty = get_HR([], data_survival)
    HR_full = get_HR(name_rules, data_survival)
    dHRs = []

    # Shapley Computation
    n_rules = len(name_rules)
    name_rules = np.array(name_rules)
    for m in range(iter_max):
        dHR = np.zeros([n_rules])
        idx = np.random.permutation(n_rules)
        HRs = [HR_empty]
        for i_rule in range(1, n_rules):
            name_rules_subset = name_rules[idx][:i_rule]
            HR = get_HR(name_rules_subset, data_survival)
            HRs.append(HR)
        HRs.append(HR_full)
        dHR[idx] = np.array([HRs[i]-HRs[i-1] for i in range(1, len(HRs))])
        dHRs.append(dHR)
        # Convergence checking
        SEM = compute_SEM(dHRs)
        if verbose:
            print('Shapley Computation Iteration %d | SEM = %.4f' % (m, SEM))
        if (m>0) and (SEM < tolerance):
            print('Stopping criteria satisfied!')
            break
    if m == (iter_max-1):
        print('Maximum iteration reached!')
    shapley_value = np.mean(dHRs, axis=0)

    return shapley_value
    