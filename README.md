# ICL_Lucy_Cornes_MRes_Macaque_Code
Analysis of associations between maternal traits and the expression of same-sex sexual behaviour (SSB) in rhesus macaques.
# README: Rhesus Macaque Maternal Effects on Same-Sex Sexual Behaviour Study

## Project Overview

This repository contains R scripts for analysing the associations between maternal traits, female fecundity (measured by offspring count) and the expression of same-sex sexual behaviour (SSB) in male rhesus macaque offspring. All IDs have been randomised from their original ID.

## Research Questions

1. **H1: Maternal Care → Maternal Fecundity**: Do mothers with higher maternal care rates have more offspring?
2. **H2A: Social Connectedness → Maternal Fecundity**: Do more socially connected mothers have more offspring?
3. **H2B: Maternal Care → Sons' SSB**: Do sons of mothers with higher maternal care exhibit more SSB?
4. **H2C: Social Connectedness → Sons' SSB**: Do sons of more socially connected mothers exhibit more SSB?
5. **H3: Facial Pigmentation → Sons' SSB**: Do sons of mothers with redder faces exhibit more SSB?
6. **Supplementary**: Direct relationship between mothers' SSB and fecundity

## Workflow Overview

### Phase 1: Data Preparation (Scripts 01-02)
### Phase 2: Behavioural Categorisation (Scripts 03-05)
### Phase 3: PCA Analysis (Scripts 6a-6b)
### Phase 4: Main Hypotheses Testing (Scripts 6C-9)
### Phase 5: Supplementary Analyses (Script 10a-b)

---

## Script Details

### **01_Female_Data_Clean.R**
**Purpose**: Initial data cleaning and duration scaling
- Imports raw focal follow data (`All_Focals_Combined.csv`)
- Removes unnecessary columns and cleans text fields
- Scales focal durations from 7200s to 9000s for standardisation for David Score
- **Outputs**: 
  - `New_All_Focals_Cleaned.csv`
  - `focal_data_scaled.csv`
  - `focal_duration_check.csv`

### **02_Behaviour_Time_Summary.R**
**Purpose**: Creates behavioural summaries and integrates demographic data
- **Inputs**: `focal_data_scaled.csv`, `Infant_count.csv`, `Offspring Count.csv`
- Summarises time spent on each behaviour per individual
- Adds age and offspring count information
- Performs data quality checks
- **Outputs**:
  - `behaviour_summary.csv`
  - `female_behaviour.csv`

### **03_Facial_Pigmentation_Extraction.R**
**Purpose**: Extracts facial redness measures from photos
- Processes facial photos to extract redness metrics
- Applies lighting corrections
- Creates multiple redness indices
- **Outputs**: `facial_pigmentation_data.csv`
- **Note**: Requires manual photo directory path setup and manual addition of dates

### **04_Behaviour_categorisation.R**
**Purpose**: Categorises behaviours into maternal care and social connectedness
- **Inputs**: `focal_data_scaled.csv`, demographic files
- Defines 14 maternal care categories and 12 social connectedness categories
- Calculates rates per hour for each category
- Filters to mothers with infants <1 year old
- **Outputs**:
  - `Updated_female_beh_rate.csv`
  - `maternal_care_analysis_data.csv`

### **05_VIF_Analysis.R**
**Purpose**: Tests for multicollinearity in behavioural categories
- **Inputs**: `Updated_female_beh_rate.csv`
- Calculates Variance Inflation Factor (VIF) scores
- Assesses correlations between behavioural categories
- **Outputs**: VIF tables and correlation assessments

### **06_Female_David_Score.R**
**Purpose**: Calculates dominance rankings using David's Score
- **Inputs**: `focal_data_scaled.csv`
- Analyses aggressive interactions
- Computes David's Scores by social group
- **Outputs**: `Female_DS_byGroup.csv`

### **7A_Maternal_PCA.R**
**Purpose**: Principal Component Analysis for maternal care behaviours
- **Inputs**: `maternal_behaviour_rate.csv`, `Female_DS_byGroup.csv`
- Creates maternal care PCA scores
- Integrates with David's Scores (dominance rankings)
- Tests PCA scores against offspring count
- **Outputs**: `new_maternal_data_with_PCA_and_DavidScore.csv`

### **7B_Social_PCA.R**
**Purpose**: PCA for social connectedness behaviours
- **Inputs**: `social_connectedness.csv`, `Female_DS_byGroup.csv`
- Creates social connectedness PCA scores
- Tests against reproductive success
- **Outputs**: `social_data_with_PCA_and_DavidScore.csv`

### **7C_Maternal_Fecundity.R**
**Purpose**: H1 - Tests maternal care behaviours against maternal fecundity
- **Inputs**: `New_All_Focals_Cleaned.csv`, demographic files
- Stepwise model selection for 14 maternal care categories
- Uses Poisson/Negative Binomial models with age offset
- **Outputs**: `maternal_care_analysis_data_updated.csv`

### **7D_Social_Fecundity.R**
**Purpose**: H2A - Tests social connectedness against maternal fecundity
- Similar approach to 6C but for social behaviours
- **Outputs**: `social_connectedness.csv`

### **8A_PCA_maternal_SSB.R**
**Purpose**: Tests maternal PCA scores against sons' SSB
- **Inputs**: PCA results, male behaviour data, mother-son linkages
- Mixed-effects models with random intercepts for mothers
- Weighted by proportion of sons sampled
- **Note**: Uses Tweedie family for SSB distribution

### **8B_Maternal_independent_SSB.R**
**Purpose**: H2B - Individual maternal care behaviours → sons' SSB
- **Inputs**: `maternal_behaviour_rate.csv`, male behavioural data
- Tests each of 14 maternal care behaviours independently
- Standardises all predictors for effect size comparison
- Uses Tweedie mixed-effects models
- **Key Results**: Muzzle contact and retrieving offspring significant
- **Outputs**: 
  - `maternal_care_SSB_analysis_data_tweedie.csv`
  - `maternal_care_SSB_results_summary_tweedie.csv`

### **8C_maternal_multivariate_SSB.R**
**Purpose**: Multivariate analysis of maternal care behaviours
- Tests all maternal behaviours simultaneously in one model
- Backwards stepwise elimination for model selection
- Controls for intercorrelations between behaviours
- **Key Results**: Muzzle contact remains significant in multivariate model

### **8D_PCA_social_SSB.R**
**Purpose**: Tests social connectedness PCA scores against sons' SSB
- Parallel analysis to 8a but for social behaviours

### **8E_Social_connectedness_independent_ssb.R**
**Purpose**: H2C - Individual social behaviours → sons' SSB
- **Inputs**: `social_connectedness.csv`, male behavioural data
- Tests each social behaviour independently
- **Key Results**: Initiating contact behaviour significant
- **Outputs**: Standardised analysis results

### **8F_social_multivariate_SSB.R**
**Purpose**: Multivariate analysis of social connectedness behaviours
- Parallel to maternal multivariate analysis
- **Key Results**: Initiating contact and teeth chatter significant

### **9_Facial_Colouration_SSB.R**
**Purpose**: H3 - Tests maternal facial redness against sons' SSB
- **Inputs**: `redness_data.csv`, male behavioural data
- Multiple redness indices tested
- **Note**: Requires facial pigmentation data from script 03

### **10_SI_SSB_Fecundity_Mothers.R**
**Purpose**: Supplementary analysis - mothers' own SSB → fecundity
- **Inputs**: `mother_master_sheet2025.csv`
- Tests direct relationship between female SSB and offspring count

### **10_SIHistogram_Time.R**
**Purpose**: Creates temporal distribution plots for focal follows
- **Inputs**: `focal_count.csv`
- Visualises observation timing patterns
- **Outputs**: `focal_follow_histograms.png`

---

## Key Data Files Required

### Input Files (must be provided):
- `All_Focals_Combined.csv` - Raw focal follow data
- `Infant_count.csv` - Infant age information  
- `Offspring Count.csv` - Maternal reproductive history
- `male_behaviour_master_sheet2025.csv` - Male SSB scores
- `mother_master_sheet2025.csv` - Mother-son linkages
- `focal_count.csv` - Observation timing data

### Optional Files:
- Photo directory for facial pigmentation analysis
- `redness_data.csv` - Pre-extracted facial redness data

### Additional Data Files

#### Facial Pigmentation Data
The repository includes processed facial pigmentation data (`facial_pigmentation_data.csv`) extracted from photographs of study females. This dataset contains:

- **Individual IDs**: Anonymised identifiers matching those used across all other datasets
- **Redness metrics**: Multiple measures of facial redness including raw RGB values, red ratios, and lighting-corrected measures
- **Photo metadata**: Image quality indicators, brightness levels, and exposure metrics
- **Photo dates**: Manually added dates indicating when each photograph was taken

**Important note**: The photo dates in this dataset were manually entered after the automated facial pigmentation extraction process. Each photograph was individually dated to ensure accurate temporal alignment with behavioural observations and reproductive data.

#### Facial Pigmentation Analysis
The facial pigmentation data allows for investigation of potential relationships between facial colouration and:
- Reproductive success (offspring count)
- Behavioural patterns
- Maternal characteristics

The dataset includes both raw redness measures and lighting-corrected values to account for variations in photographic conditions across different images.

---

## Statistical Approach

### Model Families Used:
- **Tweedie**: For SSB outcomes (handles zeros and continuous values)
- **Poisson/Negative Binomial**: For offspring counts
- **Linear Mixed Models**: For normally distributed outcomes

### Key Features:
- **Random Effects**: Mother ID for clustering
- **Weights**: Proportion of sons sampled per mother
- **Standardisation**: All behavioural predictors standardised (mean=0, SD=1)
- **Age Offsets**: For reproductive analyses

### Effect Size Interpretation:
- Results show percent change per 1 standard deviation increase
- Standardisation allows comparison across different behaviours
- 95% confidence intervals provided for all effects

---

## Key Findings Summary

### Significant Effects:
1. **Muzzle contact** behaviour → increased sons' SSB
2. **Retrieving offspring** behaviour → decreased sons' SSB  
3. **Initiating contact** (social) → increased sons' SSB
4. **Teeth chatter** (social) → increased sons' SSB

### Model Validation:
- All models include comprehensive diagnostics
- DHARMa residual testing
- Multicollinearity assessments (VIF)
- Model comparison using AIC

---

## Usage Instructions

### Prerequisites:
```r
# Required R packages
install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "glmmTMB", 
                   "DHARMa", "corrplot", "factoextra", "car", "MASS", "statmod"))
```

### Execution Order:
1. **Data Preparation**: Run scripts 01-02 sequentially
2. **Behavioural Processing**: Run scripts 03-05 (03 optional if no photos)
3. **Hypothesis Testing**: Run scripts 6a-8E based on research questions
4. **Supplementary**: Run scripts 9-10 as needed

### File Dependencies:
- Each script lists its required input files at the top
- Output files from earlier scripts become inputs for later analyses
- Some scripts can be run independently if intermediate files exist

---

## Troubleshooting

### Common Issues:
1. **Missing Files**: Check that previous scripts completed successfully
2. **ID Matching**: Ensure consistent ID formats across datasets
3. **Memory**: Large datasets may require increased R memory allocation
4. **Packages**: Some advanced packages may need manual installation

### Model Convergence:
- If Tweedie models fail, try restarting R and ensure packages are up to date

---

## Output Interpretation

### Effect Sizes:
- **Small Effect**: <10% change
- **Medium Effect**: 10-30% change  
- **Large Effect**: >30% change

### Statistical Significance:
- **p < 0.05**: Significant
- **p < 0.01**: Highly significant
- **p < 0.001**: Very highly significant

### Model Quality:
- **AIC**: Lower values indicate better models
- **R²**: Proportion of variance explained
- **Diagnostic p-values >0.05**: Good model fit

---



*Last updated: [Date]*
*Version: 1.0*
