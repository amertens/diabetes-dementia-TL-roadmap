# diabetes-dementia-TL-roadmap

File structure: 

- Cohort_creation/Code/
  - 01_diabetes_cohort.R- creates general diabetes cohort based only on exposure criteria
  - 02_dementia_exclusions.R - creates a file of additional exclusion variables based on dementia dx, insulin prior to index
  - [03_data-structure.R--not created yet, will be the data structure code from Zeyi]
- Data/
  - reference- this is where I put reference files used to create the cohort-- dx codes, rx codes, anything else we'll need to import in and reference
  - tmp- temporary data files (things deleted upon each re-running of code)output
  - output/finalized datasets
- Data summaries/ - this is where we have Thomas's data distributions
- Example code/ - has an ltmle example and some cleaning code form previous cleaning
- data_raw/ - I'm putting the raw heaven simulated data here
- presentation/ - everything from the most recent talk

