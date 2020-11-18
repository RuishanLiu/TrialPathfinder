import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import datetime
import os, sys
import matplotlib.pyplot as plt 

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors

from lifelines import CoxPHFitter, KaplanMeierFitter

    
def emulate_trials(cohort, features, drug_treatment, drug_control, name_rules, covariates_cont=[], covariates_cat=[], ps_method='IPTW', thresh_censor=None, name_DrugName='DrugName', name_StartDate='StartDate', name_OutcomeDate='OutcomeDate', name_LastVisitDate='LastVisitDate', verbose=False):
    '''
    Emulate trial
    Return HR, confidence interval, data fit for the cox model
    '''
    MISSING = 'missing'
    name_PatientID = cohort.name_PatientID
    
    def generate_survival_data(df_raw, covariates=None, thresh_censor=None, dateofevent='dateofdeath'): 
        '''
        Generate survival data from features. Add column event - whether died or censored.
        '''
        df = df_raw[[name_PatientID, name_DrugName]].copy()

        df['duration_recomputed'] = MISSING

        # Add event. True - died.  False - censored.
        df.loc[df.index, 'event'] = False
        ids_death = df_raw[dateofevent]!=MISSING
        df.loc[ids_death, 'event'] = True
        # death
        inds = ids_death
        d_days = pd.to_datetime(df_raw[dateofevent].loc[inds]) - pd.to_datetime(df_raw[name_StartDate].loc[inds])
        df.loc[inds, 'duration_recomputed'] = [d_day.days for d_day in d_days]
        # not death and has last_vist
        inds = ~ids_death & (df_raw[name_LastVisitDate]!=MISSING)
        d_days = pd.to_datetime(df_raw[name_LastVisitDate].loc[inds]) - pd.to_datetime(df_raw[name_StartDate].loc[inds])
        df.loc[inds, 'duration_recomputed'] = [d_day.days for d_day in d_days]

        if thresh_censor is not None:
            inds_longer = df['duration_recomputed']>thresh_censor
            df.loc[inds_longer, 'duration_recomputed'] = thresh_censor
            df.loc[inds_longer, 'event'] = False

        # Add covariates
        if covariates is not None:
            df = pd.merge(df, df_raw[[name_PatientID]+covariates].copy(), how='left', on=name_PatientID)
        return df


    def generate_trial_data(df_raw, arm_exp, arm_control, covariates_cont=[], covariates_cat=[]):
        '''
        Generate data for each trial. 
        arm_exp: treatment names in experiment
        arm_control: treatment names in control
        Function:
            - One-hot for categorical covariates and standardize continuous covariates
            - Add a column of treatment. 1-treatment; 0-control
        Return data with columns [duration, event, treatemnt, all the covariates]
        '''
        df = df_raw[['duration_recomputed', 'event']].copy()

        # Add a column of treatment.  Treated - 1; Control - 0; Not included - -1
        df.loc[df.index, 'treatment'] = -1
        for linename in arm_exp:
            df.loc[df_raw[name_DrugName]==linename, 'treatment'] = 1
        for linename in arm_control:
            df.loc[df_raw[name_DrugName]==linename, 'treatment'] = 0

        # Add covariates
        df = pd.merge(df, df_raw[covariates_cont+covariates_cat].copy(), left_index=True, right_index=True)

        # Subset which are included in this arm
        df = df.loc[df['treatment']!=-1]

        # One-hot encoding Categorical Variables
        if covariates_cat != []:
            data_categorical = df[covariates_cat].copy()
            df = df.drop(covariates_cat, axis=1)
            one_hot = pd.get_dummies(data_categorical, drop_first=True)
            df = pd.merge(df, one_hot, left_index=True, right_index=True)

        # Standardize conticuous one and set missing to be 0
        for col in covariates_cont:
            idx_notmiss = (df[col]!=MISSING)
            df.loc[~idx_notmiss, col] = 0
            if np.sum(idx_notmiss) > 0:
                data_continuos = np.array(df.loc[idx_notmiss, col])
                data_continuos = StandardScaler().fit_transform(data_continuos.reshape(-1, 1)).reshape(-1)
                df.loc[idx_notmiss, col] = data_continuos

        return df

    def generate_cox_data(data_trial, ps_method='IPTW', verbose=1):
        '''
        ps_method = 'IPTW' or 'Match'
        Return 
        - 'IPTW_stabilized' 'IPTW_unstabilized': [duration, event, treatemnt, stabilized IP weights]
        - 'Match': [duration, event, treatemnt] after propensity score matching
        '''
        df = data_trial.iloc[:, :3]
        model = LogisticRegression(solver='lbfgs', n_jobs=-1, class_weight='balanced')
        X = np.array(data_trial.iloc[:, 3:])
        y = np.array(data_trial['treatment'])
        model.fit(X, y)
        p_treated = float(np.sum(y==1))/y.shape[0]
        propensity_scores = model.predict_proba(X)[:, 1]

        sample_weight = np.ones(y.shape[0])
        sample_weight[y==0] = p_treated / (1-p_treated) 
        if verbose:
            print('Score for Logistic Model: %.4f (sample balanced)' % (model.score(X, y, sample_weight=sample_weight)))

        if 'IPTW' in ps_method:
            # Stabilized IP Weights
            df.loc[df.index, 'weights'] = 0
            IP_treated = 1 / propensity_scores
            IP_untreated = 1 / (1 - propensity_scores)
            if 'stabilized' in ps_method:
                IP_treated = p_treated * IP_treated
                IP_untreated = (1 - p_treated) * IP_untreated          
            df.loc[df['treatment']==1, 'weights'] = IP_treated[df['treatment']==1]
            df.loc[df['treatment']==0, 'weights'] = IP_untreated[df['treatment']==0]
        elif ps_method == 'Match':
            # Propensity Matching
            idx_untreated = df.loc[df['treatment']==0].index
            ps_untreated = 1-propensity_scores[df['treatment']==0]
            ps_treated = propensity_scores[df['treatment']==1]
            # Matching
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(ps_untreated.reshape(-1, 1))
            _, indices = nbrs.kneighbors(ps_treated.reshape(-1, 1))
            indices = indices.reshape(-1)
            df = pd.concat([df.loc[df['treatment']==1].copy(), df.loc[idx_untreated[indices]]])

        return df

    def generate_trial_cox(cohort, data_survival, arm_exp, arm_control, name_rules=[], ps_method='IPTW', covariates_cont=[], covariates_cat=[], verbose=1):
        '''
        Generate Trial data and Cox data from cohort and survival data
        ps_method: 'IPTW' or 'Match'
        If name_rules[0] == -1: the excluded cohort
        '''
        if len(name_rules) > 0:
            if name_rules[0] != -1:
                cohort.selection(name_rules)
                data_cohort = data_survival.loc[data_survival[cohort.name_PatientID].isin(cohort.select_id)]
            else:
                cohort.selection(name_rules[1:])
                patients_id = cohort.select_id
                data_cohort = data_survival.loc[~data_survival[name_PatientID].isin(patients_id)]
        else:
            data_cohort = data_survival

        # Data for this Trial
        data_trial = generate_trial_data(data_cohort, arm_exp, arm_control, covariates_cont, covariates_cat)

        # Data For Cox
        data_cox = generate_cox_data(data_trial, ps_method, verbose=verbose)

        if verbose:
            print('Total data shape:', data_survival.shape, 'Selected cohort data shape:', data_cohort.shape)
            print('Selected Data for the trial: %d | Experiment: %d | Control: %d' % (data_trial.shape[0], np.sum(data_trial['treatment']==1), np.sum(data_trial['treatment']==0)))
            print('Final Data for Cox model | Experiment: %d | Control: %d' % (np.sum(data_cox['treatment']==1), np.sum(data_cox['treatment']==0)))
            print('Final Data for Cox model - Event (not censored) | Experiment: %d | Control: %d' % (np.sum((data_cox['treatment']==1)&(data_cox['event']==True)), (np.sum((data_cox['treatment']==0)&(data_cox['event']==True)))))

        return data_cox, data_trial

    def cox(data_cox, data_trial, ps_method='IPTW'):
        '''
        Run Cox Proportional Hazard Model
        Return: Hazard Ratio
        ps_method: 'IPTW' or 'Match'
        '''
        # Run Cox
        cph = CoxPHFitter()
        if 'IPTW' in ps_method:
            cph.fit(data_cox, 'duration_recomputed', 'event', weights_col='weights', robust=True)
        elif ps_method == 'Match':
            cph.fit(data_cox, 'duration_recomputed', 'event', robust=True)

        HR = cph.hazard_ratios_['treatment']
        return HR, cph
    
    
    covariates = covariates_cont + covariates_cat
    
    data_survival = generate_survival_data(features.copy(), covariates=covariates, thresh_censor=thresh_censor, dateofevent=name_OutcomeDate) 
    data_survival = data_survival.loc[data_survival['duration_recomputed']!=MISSING]
    
    HR, cph = cox(data_cox, data_trial, ps_method=ps_method)
    HR_conf = np.exp(np.array(cph.confidence_intervals_).reshape(-1))
    
    return HR, HR_conf, data_cox
    