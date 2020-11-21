# Trial PathFinder

Python library for Trial Pathfinder, an AI framework to systematically evaluate clinical trial eligibility criteria. Functions provided by this package: encoding eligibility criteria, emulating existing trials under combinations of eligibility criteria, evaluating individual eligibiliey rule and suggesting data-driven criteria.

*Working paper. More background and information will be provided when the paper is ready.*

This repository contains tutorials and examples to use the library; see [tutorials/tutorial.ipynb](https://github.com/RuishanLiu/TrialPathfinder/blob/master/tutorial/tutorial.ipynb).


# Installation

The package TrialPathfinder is available on PyPI -
```shell
pip install TrialPathfinder
```
    
We also provide the option for manual installation: download this Github repository and run
```shell
cd TrialPathfinder/
python setup.py install --user
```

# Quick Guidance

Here we give a quick guidance of using TrialPathfinder. More details see [tutorial.ipynb](https://github.com/RuishanLiu/TrialPathfinder/blob/master/tutorial/tutorial.ipynb).

```python
import TrialPathfinder as tp

###### Encode Eligibility Criteria #####

# Create cohort selection object
cohort = tp.cohort_selection(patientids, name_PatientID='PatientID')
# Add the data tables needed in the eligibility criterion
cohort.add_table(name_table1, table1)
# Add individual eligibility criterion
cohort.add_rule(rule1)

###### Emulate Existing Trials and Survival Analysis ######

# Given a combination of eligibility rules names_rules (an empty list name_rules=[] indicates fully-relaxed criteria)).
HR, CI, data_cox = tp.emulate_trials(cohort, features, drug_treatment, drug_control, name_rules)

###### Evalute Individual Criterion ######

# Return the Shapley values for each rule in names_rules
shapley_values = tp.shapley_computation(cohort, features, drug_treatment, drug_control, names_rules)

###### Criteria Relaxation - Data-driven Criteria ######

# Select all the rules with Shapley value less than 0.
names_rules_relax = names_rules[shapley_values<0.]
# Survival analysis on the data-driven criteria
HR, CI, data_cox = tp.emulate_trials(cohort, features, drug_treatment, drug_control, name_rules_relax)
```

# Documentation


We highly recommend reading the [tutorial.ipynb](https://github.com/RuishanLiu/TrialPathfinder/blob/master/tutorial/tutorial.ipynb).


## 1. Data Requirements

TrialPathfinder reads tables in Pandas dataframe structure (pd.dataframe) as default. The date information should be read as datetime (use function pd.to_datetime to convert if not).

**1. Features**:
- <font color=darkblue>*Patient ID*</font>
- Treatment Information
    - <font color=darkblue>*Drug name*</font>.
    - <font color=darkblue>*Start date*</font>.
    - <font color=darkblue>*Date of outcome*</font>. For example, if overall survival (OS) is used as metric, the date of outcome is the date of death. If progression-free survival (PFS) is used as metric, the date of outcome is the date of progression.
    - <font color=darkblue>*Date of last visit*</font>. The patient's last record date of visit, used for censoring.
- <font color=darkblue>*Covariates (optional)*</font>: adjusted to emulate the blind assignment, used by Inverse probability of treatment weighting (IPTW) or propensity score matching (PSM). Some examples: age, gender, composite race/ethnicity, histology, smoking status, staging, ECOG, and biomarkers status.

**2. Tables used by eligibility criteria.**
- Use the same Patient ID as the features table.


## 2. Stadards of encoding eligibility criteria

We built a computational workflow to encode the description of eligibility criteria in the protocols into standardized instructions which can be parsed by Trial Pathfinder for cohort selection use. 

**1. Basic logic.**

- Name of the criteria is written in the first row.
- A new statement starts with “#inclusion” or “#exclusion” to indicate the criterion’s type. Whether to include patients who have missing entries in the criteria: “(missing include)” or “(missing exclude)”. The default choice is including patients with missing entries. 
- Data name format: “Table[‘featurename’]”. For example, “demographics[‘birthdate’]” denotes column date of birth in table demographics.
- Equation: ==, !=, <, <=, >, >=. 
- Logic: AND, OR.
- Other operations: MIN, MAX, ABS.
- Time is encoded as “DAYS(80)”: 80 days; “MONTHS(4)”: 4 months; “YEARS(3)”: 3 years.

*Example: criteria "Age" - include patients more than 18 years old when they received the treatment.*

> Age \
\#Inclusion \
features['StartDate'] >= demographics['BirthDate'] + @YEARS(18> 


**2. Complex rule with hierachy.**
- Each row is operated in sequential order
    - The tables are prepared before the last row. 
    - The patients are selected at the last row. 

*Example: criteria "Platelets" - include patients whose platelet count ≥ 100 x 10^3/μL*. \
To encode this criterion, we follow the procedure: 
1. Prepare the lab table: 
    1. Pick the lab tests for platelet count
    2. The lab test date should be within a -28 to +0 window around the treatment start date
    3. Use the record closest to the treatment start date to do selection.
2. Select patients: lab value larger than 100 x 10^3/μL.
> Platelets \
\#Inclusion \
(lab['LabName'] == 'Platelet count') \
(lab['TestDate'] >= features['StartDate'] - @DAYS(28) ) AND (lab['TestDate'] <= features['StartDate']) \
MIN(ABS(lab['TestDate'] - features['StartDate'])) \
lab['LabValue'] >= 100 
