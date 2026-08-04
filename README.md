Data and code for Verma et al: Parasite-induced hormonal manipulation associated with severe male reproductive costs in wild primate hosts

There are three R code files and four data files for this paper. The code specifies which file is required.

### Breakdown:
1. Plasma hormone figures and analyses: plasma_hormones.csv & plasma_hormones_figs_analyses.Rmd

- DRT_ID: Individual ID
- Mean: Mean plasma hormone concentration of two replicates
- Std.Dev: Standard Deviation of plasma hormone concentration
- CV: Coefficient of Variation across replicates
- assay_type: Testosterone or estradiol assay
- Sex: Sex of individual
- Age: Age category of individual (Adult or Subadult)
- male_status: Social status of male (bachelor, follower, or natal)
- positive: Infection status (0/1)
- Body_mass_kg: Body mass in kg
- cyst: Cyst presence (0/1)
- infection: Recoded infection variable
- assay_type2: Recoded assay type with unit
- repro.state.2: Female reproductive state (immature, cycling, pregnant, or unknown)

2. Fecal androgen metabolite figures and analyses: fecal_androgen_metabolites.csv & fecal_androgen_metabolites_figs_analyses.Rmd

- code: Unique individual ID
- sample_year: Year of faecal sample collection
- sex: Sex of individual
- CV: Coefficient of Variation across faecal androgen metabolite replicates
- T_ng_g: Faecal androgen metabolite concentration in ng/g
- age_at_sample: Age of individual at time of faecal sample collection
- cyst: Cyst presence (0/1)
- positive: Infection status (0/1)
- sample_type: Sample type used to infer infection status (urine/plasma)
- male_status: Social status of male (bachelor, follower, or leader)
- season: Season at time of faecal sample collection (Hot Dry/Cold Dry/Cold Wet)
- rain_90: 90-day rainfall (mm)

3. Tenure survival analyses and figures: tenures.csv, taenia_infection.csv & survival_analysis_figs.Rmd
