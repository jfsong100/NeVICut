# NeVI-Cut: Neural Variational Inference for Cut-Bayes

**NeVI-Cut** implements neural variational inference for cut-Bayes. It supports multiple flow families (RQ-NSF(AR),RQ-NSF(C),UMNN) and includes reproducible simulations and real-data case studies.

Authors: Jiafang Song, Sandipan Pramanik, Abhirup Datta

Acknowledgement: This work was developed with support from the National Institute of Environmental Health Sciences (NIEHS) (grant R01ES033739); Gates Foundation Grants  INV-034842 and INV-070577; Johns Hopkins University Institute for Data Intensive Engineering and Science Seed Funding; and K99 NIH Pathway to Independence Award (1K99HD114884-01A1).


## Setup
### Option A: Conda environment (recommended)
```bash
conda env create -f environment.yml
conda activate nevi-cut
```

### Option B: Pip install via setup.py:
```bash
pip install .
# or for development
pip install -e .
```

## Repository Layout
```
cutBayesFlow/          
experiments/            
  simulations/
    theory_justification.ipynb
    biased_normal.ipynb
    propensity_score/
  real_data/
    HPV/
    COMSA/
```

### Repository Description

#### Core Implementation
- **`cutBayesFlow/`**
  - `model.py`: `CutBayesFlow` with flow options
    - **NSF-AR**: Rational-Quadratic Neural Spline Flow, autoregressive (Durkan et al.,2019)
    - **NSF-C**: Rational-Quadratic Neural Spline Flow, coupling (Durkan et al.,2019)
     - **UMNN**: Unconstrained Monotonic Neural Network (Wehenkel and Louppe, 2019)
  - `train.py`: training loop

#### Simulation Studies
- **`experiments/simulations/`**
  - **`Theory_justification.ipynb`**: demonstration of NeVI-Cut with UMNN and RQ-NSF(AR). (Paper §5.1)
  - **`biased_normal.ipynb`**: misspecified prior experiment. (Paper §5.2)
  - **`propensity_score/`**: synthetic upstream propensity-score posteriors and downstream NeVI-Cut fits under outcome-model misspecification. (Paper §5.3)

#### Real-Data Applications
- **`experiments/real_data/HPV/`**
  - Association between HPV prevalence and cervical cancer incidence (motivated by Maucort-Boulch et al., 2008).  
  - Upstream posteriors via conjugate Beta models (`stage 1/`).  
  - Downstream analysis via NeVI-Cut (`hpv.ipynb`). 
  - See Paper § 6.1 for full analysis.

- **`experiments/real_data/COMSA/`**
  - `comsa_data`: **posterior samples of confusion matrices** and **downstream cause-of-death counts** (COMSA–Mozambique).
        **downstream cause-of-death counts** 
  - Upstream data: confusion-matrix posteriors are publicly available at the [CCVA-Misclassification-Matrices repository](https://github.com/sandy-pramanik/CCVA-Misclassification-Matrices).  
  - Downstream data: counts come from *CCVA Outputs for Publicly Available Verbal Autopsy (VA) Data from COMSA–Mozambique*, which include results from three computer-coded verbal autopsy (CCVA) algorithms, **EAVA**, **InSilicoVA**, **InterVA**. Outputs for EAVA are generated using the EAVA R package, while InSilicoVA and InterVA outputs are produced with the openVA R package.
  - Downstream analysis via NeVI-Cut (`NeVI_Cut.ipynb`). 
  - Downstream analysis via Parametric Variational Inference (`Parametric_VI.ipynb`). 
  - See Paper § 6.2 for the full COMSA application.

## Citation
Please cite the following paper when you use NeVI-Cut:


## References
- Durkan, C., et al. (2019). *Neural Spline Flows.* *Advances in Neural Information Processing Systems,* 32.
- Wehenkel, A., & Louppe, G. (2019). *Unconstrained Monotonic Neural Networks.* *Advances in Neural Information Processing Systems,* 32.  
- Plummer, M. (2015). *Cuts in Bayesian Graphical Models.* *Statistics and Computing,* 25: 37–43.
- Maucort-Boulch, D., et al. (2008). *International Correlation Between Human Papillomavirus Prevalence and Cervical Cancer Incidence.* *Cancer Epidemiology Biomarkers & Prevention,* 17 (3): 717–720.
- Wilson, E., et al. (2025). **EAVA**: Deterministic Verbal Autopsy Coding with Expert Algorithm Verbal Autopsy. [R package](https://doi.org/10.32614/CRAN.package.EAVA).  
- Li, Z. R., et al. (2024). **openVA**: Automated Methods for Verbal Autopsy. [R package](https://cran.r-project.org/web/packages/openVA/index.html).  
- Pramanik, S., et al. (2025). *Country-Specific Estimates of Misclassification Rates of Computer-Coded Verbal Autopsy Algorithms.* *medRxiv.*  
- Pramanik, S., et al. (2025). *Modeling Structure and Country-Specific Heterogeneity in Misclassification Matrices of Verbal Autopsy-Based Cause-of-Death Classifiers.* *Annals of Applied Statistics,* 19 (2): 1214–1239.  
- McCormick, T. H., et al. (2016). *Probabilistic Cause-of-Death Assignment Using Verbal Autopsies.* *JASA,* 111 (515): 1036–1049.  
- Byass, P., et al. (2012). *Strengthening Standardised Interpretation of Verbal Autopsy Data: the InterVA-4 Tool.* *Global Health Action,* 5, 19281.  
- Macicame, I., et al. (2023). *Countrywide Mortality Surveillance for Action in Mozambique.* *American Journal of Tropical Medicine and Hygiene,* 108 (Suppl 5): 5–16.  
- Kalter, H. D., et al. (2015). *Direct Estimates of National Neonatal and Child Cause-Specific Mortality Proportions in Niger.* *Journal of Global Health,* 5.
