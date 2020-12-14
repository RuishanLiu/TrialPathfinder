# Tutorial for TrialPathFinder


[tutorial/tutorial.ipynb](https://github.com/RuishanLiu/TrialPathfinder/blob/master/tutorial/tutorial.ipynb) provides a detailed tutorial and example to use the library TrialPathFinder.

### Pipeline

- Explaination about data requirement
- Stadards of encoding eligibility criteria and examples
- How to use TrialPathfinder (*tp*) with application to synthetic data
    - Patient selection by *tp.cohort_selection()*
    - Survival analysis by *tp.emulate_trials()*
    - Shapley computation by *tp.shapley_computation()*


### Synthetic Data

The synthetic example data used by the tutorial are provided in directory [tutorial/data](https://github.com/RuishanLiu/TrialPathfinder/tree/master/tutorial/data). 
- Eligibility criteria (*'criteria.csv'*) have five rules: Age, Histology_Squamous, ECOG, Platelets, Bilirubin.
- The features (*'features.csv'*) contain the required treatment information and three covariates (gender, race, ecog). 
- Two tables (*'demographics.csv'* and *'lab.csv'*) are used by its eligibility criteria.