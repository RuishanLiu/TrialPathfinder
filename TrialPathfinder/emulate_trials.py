from .utils import *
    
def emulate_trials(cohort, features, drug_treatment, drug_control, name_rules, covariates_cont=[], covariates_cat=[], thresh_censor=None, name_DrugName='DrugName', name_StartDate='StartDate', name_OutcomeDate='OutcomeDate', name_LastVisitDate='LastVisitDate', indicator_miss='Missing'):
    '''
    Emulate trial
    Return HR, confidence interval, data fit for the cox model
    '''  
    # Generate survival information.
    data_survival = generate_survival_data(features.copy(), cohort.name_PatientID,
                                           covariates=covariates_cont+covariates_cat, thresh_censor=thresh_censor,
                                           name_DrugName=name_DrugName, name_StartDate=name_StartDate, 
                                           name_OutcomeDate=name_OutcomeDate, name_LastVisitDate=name_LastVisitDate,
                                           indicator_miss=indicator_miss) 
    
    # Generate data fit for Cox proportional hazards model.
    data_cox = generate_trial_cox(cohort, data_survival, drug_treatment, drug_control, 
                                  name_rules=name_rules, name_DrugName=name_DrugName,
                                  covariates_cont=covariates_cont, covariates_cat=covariates_cat)

    HR, CI = cox(data_cox)
    
    return HR, CI, data_cox
    