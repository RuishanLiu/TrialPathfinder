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

from lifelines import CoxPHFitter

def generate_survival_data(df_raw, name_PatientID, covariates=None, thresh_censor=None, name_DrugName='DrugName', name_StartDate='StartDate', name_OutcomeDate='OutcomeDate', name_LastVisitDate='LastVisitDate', indicator_miss='Missing'): 
    '''
    Generate survival data from features. 
    Add column event - whether outcome recorded (event=1) or censored (event=0).
    Add column duration - time from StartDate to OutcomeDate (event=1) or LastVisitDate (event=0)
    '''
    
    df = df_raw[[name_PatientID, name_DrugName]].copy()
    df['duration'] = indicator_miss

    # Add event. 
    df.loc[df.index, 'event'] = 0
    ids_death = df_raw[name_OutcomeDate]!=indicator_miss
    df.loc[ids_death, 'event'] = 1
    # death
    inds = ids_death
    d_days = pd.to_datetime(df_raw[name_OutcomeDate].loc[inds]) - pd.to_datetime(df_raw[name_StartDate].loc[inds])
    df.loc[inds, 'duration'] = [d_day.days for d_day in d_days]
    # not death and has last_vist
    inds = ~ids_death & (df_raw[name_LastVisitDate]!=indicator_miss)
    d_days = pd.to_datetime(df_raw[name_LastVisitDate].loc[inds]) - pd.to_datetime(df_raw[name_StartDate].loc[inds])
    df.loc[inds, 'duration'] = [d_day.days for d_day in d_days]
    df = df.loc[df['duration']!=indicator_miss]
    
    if thresh_censor is not None:
        inds_longer = df['duration']>thresh_censor
        df.loc[inds_longer, 'duration'] = thresh_censor
        df.loc[inds_longer, 'event'] = 0

    # Add covariates
    if covariates is not None:
        df = pd.merge(df, df_raw[[name_PatientID]+covariates].copy(), how='left', on=name_PatientID)
    return df


def generate_trial_data(df_raw, drug_treatment, drug_control, covariates_cont=[], covariates_cat=[], name_DrugName='DrugName'):
    '''
    Generate data for each trial. 
        - One-hot for categorical covariates and standardize continuous covariates
        - Add a column of treatment. 1-treatment; 0-control
    Return data with columns [duration, event, treatemnt, all the covariates]
    '''
    df = df_raw[['duration', 'event']].copy()

    # Add a column of treatment.  Treated - 1; Control - 0; Not included - -1
    df.loc[df.index, 'treatment'] = -1
    for linename in drug_treatment:
        df.loc[df_raw[name_DrugName]==linename, 'treatment'] = 1
    for linename in drug_control:
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

    # Standardize conticuous one and set indicator_miss to be 0
    for col in covariates_cont:
        idx_notmiss = (df[col]!=indicator_miss)
        df.loc[~idx_notmiss, col] = 0
        if np.sum(idx_notmiss) > 0:
            data_continuos = np.array(df.loc[idx_notmiss, col])
            data_continuos = StandardScaler().fit_transform(data_continuos.reshape(-1, 1)).reshape(-1)
            df.loc[idx_notmiss, col] = data_continuos

    return df

def generate_cox_data(data_trial):
    '''
    Generate data fit for Cox model (IPTW).
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

    # IP Weights
    df.loc[df.index, 'weights'] = 0
    IP_treated = 1 / propensity_scores
    IP_untreated = 1 / (1 - propensity_scores)
    # Stabilization
    IP_treated = p_treated * IP_treated
    IP_untreated = (1 - p_treated) * IP_untreated          
    df.loc[df['treatment']==1, 'weights'] = IP_treated[df['treatment']==1]
    df.loc[df['treatment']==0, 'weights'] = IP_untreated[df['treatment']==0]

    return df

def generate_trial_cox(cohort, data_survival, drug_treatment, drug_control, name_rules=[], covariates_cont=[], covariates_cat=[], name_DrugName='DrugName'):
    '''
    Generate Trial data and Cox data from cohort and survival data
    '''
    if len(name_rules) > 0:
        cohort.selection(name_rules)
        data_cohort = data_survival.loc[data_survival[cohort.name_PatientID].isin(cohort.select_id)]
    else:
        data_cohort = data_survival

    # Detailed Data for this Trial
    data_trial = generate_trial_data(data_cohort, drug_treatment, drug_control, covariates_cont, covariates_cat, name_DrugName)

    # Data For Cox Model
    data_cox = generate_cox_data(data_trial)

    return data_cox

def cox(data_cox):
    '''
    Run Cox Proportional Hazard Model
    Return: Hazard Ratio, Confidence Interval
    '''
    # Run Cox
    cph = CoxPHFitter()
    cph.fit(data_cox, 'duration', 'event', weights_col='weights', robust=True)

    HR = cph.hazard_ratios_['treatment']
    CI = np.exp(cph.confidence_intervals_.values.reshape(-1))
    return HR, CI
    