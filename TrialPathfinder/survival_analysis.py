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

def select_feature(cohort, features):
    '''
    return the features for selected cohort
    '''
    name_patientid = cohort.name_patientid
    return features.loc[features[name_patientid].isin(cohort.select_id)]

string_startdate = 'startdate' # 'startdate'  'startdate_drug'
def generate_survival_data(df_raw, covariates_col=None, thresh_censor=None, dateofevent='dateofdeath'): # 'dateofdeath' 'progressiondate_toend'
    '''
    Generate survival data from features. Add column event - whether died or censored.
    '''
    df = df_raw[['patientid', 'linename']].copy()
    
    df['duration_recomputed'] = MISSING

    # Add event. True - died.  False - censored.
    df.loc[df.index, 'event'] = False
    ids_death = df_raw[dateofevent]!=MISSING
    df.loc[ids_death, 'event'] = True
    # death
    inds = ids_death
    d_days = pd.to_datetime(df_raw[dateofevent].loc[inds]) - pd.to_datetime(df_raw[string_startdate].loc[inds])
    df.loc[inds, 'duration_recomputed'] = [d_day.days for d_day in d_days]
    # not death and has last_vist
    inds = ~ids_death & (df_raw['last_visit']!=MISSING)
    d_days = pd.to_datetime(df_raw['last_visit'].loc[inds]) - pd.to_datetime(df_raw[string_startdate].loc[inds])
    df.loc[inds, 'duration_recomputed'] = [d_day.days for d_day in d_days]
    
    if thresh_censor is not None:
        inds_longer = df['duration_recomputed']>thresh_censor
        df.loc[inds_longer, 'duration_recomputed'] = thresh_censor
        df.loc[inds_longer, 'event'] = False

    # Add covariates
    if covariates_col is not None:
        df = pd.merge(df, df_raw[['patientid']+covariates_col].copy(), how='left', on='patientid')
    return df


def generate_trial_data(df_raw, arm_exp, arm_control, continuous_col=[], categorical_col=[]):
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
        df.loc[df_raw['linename']==linename, 'treatment'] = 1
    for linename in arm_control:
        df.loc[df_raw['linename']==linename, 'treatment'] = 0

    # Add covariates
    df = pd.merge(df, df_raw[continuous_col+categorical_col].copy(), left_index=True, right_index=True)

    # Subset which are included in this arm
    df = df.loc[df['treatment']!=-1]

    # One-hot encoding Categorical Variables
    if categorical_col != []:
        data_categorical = df[categorical_col].copy()
        df = df.drop(categorical_col, axis=1)
        one_hot = pd.get_dummies(data_categorical, drop_first=True)
        df = pd.merge(df, one_hot, left_index=True, right_index=True)

    # Standardize conticuous one and set missing to be 0
    for col in continuous_col:
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

def generate_trial_cox(cohort, data_survival, arm_exp, arm_control, name_rules=[], ps_method='IPTW', continuous_col=[], categorical_col=[], verbose=1):
    '''
    Generate Trial data and Cox data from cohort and survival data
    ps_method: 'IPTW' or 'Match'
    If name_rules[0] == -1: the excluded cohort
    '''
    if len(name_rules) > 0:
        if name_rules[0] != -1:
            cohort.selection(name_rules)
            data_cohort = select_feature(cohort, data_survival)
        else:
            cohort.selection(name_rules[1:])
            patients_id = cohort.select_id
            data_cohort = data_survival.loc[~data_survival['patientid'].isin(patients_id)]
    else:
        data_cohort = data_survival

    # Data for this Trial
    data_trial = generate_trial_data(data_cohort, arm_exp, arm_control, continuous_col, categorical_col)

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
    