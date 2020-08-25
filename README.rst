TrialPathfinder
===================

Python library for Trial Pathfinder, an AI framework to systematically evaluate clinical trial eligibility criteria. 

Functions provided by this package: encoding eligibility criteria, emulating existing trials under combinations of eligibility criteria, evaluating individual eligibiliey rule and suggesting data-driven criteria.

Working paper. More background and information will be provided when the paper is ready.

Installation
--------------------

Download this repository and run

.. code-block:: shell-session
    python setup.py


Import
--------------------

.. code-block:: python

    from TrialPathfinder import cohort_selection, emulate_trials, survival_analysis, shapley_computation
    
   
Usage
-------------------------------
This repository provides examples to use the library. See the python notebooks in the :code:`examples` folder (will be uploaded soon when the de-identification of Flatiron data structrue is ready).

Basic usage of this package:

.. code-block:: python

    from TrialPathfinder import cohort_selection, emulate_trials, survival_analysis, shapley_computation
    
    ###### Encode Eligibility Criteria #####
    # Create cohort selection object
    cohort = cohort_selection(patientids, name_patientid='patientid')
    # Add the data tables needed in the eligibility criterion
    cohort.add_table(name_table1, table1)
    # Add individual eligibility criterion
    cohort.add_rule(rule1)
    
    ###### Emulate Existing Trials ######
    # Given a combination of eligibility rules names_rules (an empty list [] indicates fully-relaxed criteria), process patients features for survival analysis (features is pandas Dataframe by default).
    data_survival = emulate_trials(cohort, features, names_rules)
    
    ###### Survival Analysis ######
    HR, confidence_interval, p_value = survival_analysis(data_survival)
    
    ###### Evalute Individual Criterion ######
    # Return the Shapley values for each rule in names_rules
    shapley_values = shapley_computation(cohort, features, drug_treatment, drug_control, names_rules)
    
    ###### Criteria Relaxation - Data-driven Criteria ######
    # Select all the rules with Shapley value less than -0.01
    names_rules_relax = names_rules[shapley_values<-0.01]
    # Survival analysis on the data-driven criteria
    data_survival_relax = emulate_trials(cohort, features, drug_treatment, drug_control, names_rules_relax)
    HR, confidence_interval, p_value = survival_analysis(data_survival_relax)
    
    
Optional Parameters
-------------------------------
**Encode Eligibility Criteria**: when adding rules with names already exist in the object :code:`cohort_selection`, people can choose to overwrite (parameter :code:`force=True`) or not add (parameter :code:`force=False`).

**Emulate Existing Trials**: :code:`emulate_trials` processes patients features given the criteria combination for survival analysis. We offer choice to do specific processing for continuous features (parameter :code:`cols_continuous`) and categorical features (parameter :code:`cols_categorical`).

**Shapley Analysis**: :code:`shapley_computation` estimates Shapley value by Monte Carlo sampling. The default iteriations is 100. To change the number of iterations, use parameter :code:`n_iter=100`.