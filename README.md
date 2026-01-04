# Interactive Visualisation of Populist Leadership Effects

This Shiny app visualizes causal estimates of the effects of populist leaders on economic growth and distributional outcomes (e.g., inequality) using the traditional Synthetic Control Method (SCM) and the Augmented Synthetic Control Method (ASCM). The outcome model employs matrix completion.  

Users can select specific populist episodes, outcomes, and donor pools from three databases (PLE, GDP, or combined). For each case, the app displays either:

1. Observed versus counterfactual outcome trajectories  
2. Estimated treatment effect series  

Available options include:

- **95% jackknife+ confidence bands** ([Barber et al., 2021, Ann. Statist.](https://arxiv.org/abs/1905.02928))  
- **Backtrace (3 years) and Backtrace (5 years):** in-time placebo exercises that counterfactually advance the treatment date by 3 or 5 years  
- **Leave-One-Donor-Out permutations:** robustness check omitting one donor at a time  
- **Donor-Placebo permutations:** placebo treatments assigned counterfactually to donor units (not available for counterfactual series)  

For the case-level perspective, the app supports a **single-panel view** or **Dual Plot View Mode** (side-by-side panels). In Dual Plot View, users can choose:

- **Shared Mode:** identical options across panels  
- **Separate Mode:** independent options per panel  

The **Pooled Perspective** button accesses an aggregate, cross-case perspective, allowing users to subset results by region and date and to summarize treatment effects using mean or median diagnostic measures.

Cases and options are restricted to those meeting pre-estimation requirements; availability varies by database and data coverage. Case names follow the scheme **ISO3 country code + treatment year**. Detailed information on individual populist cases is in the Appendix of the accompanying paper.  

---

## Data

The app uses a pre-processed dataset (`ShinyPopDat.qs`, ~85 MB) that is too large to store directly in the GitHub repository.  

The data file is hosted as a **release asset**:

- Download link: [ShinyPopDat.qs v1.0](https://github.com/fsimghub/populism-shiny/releases/download/v1.0/ShinyPopDat.qs)  

The app automatically downloads the file on first run if it is not already present in the `data/` directory.

---

## Installation & Running Locally

```r
# Clone the repository
git clone https://github.com/fsimghub/populism-shiny.git
cd populism-shiny

# Run the app. It will automatically check for and install any missing R packages
R -e "shiny::runApp(launch.browser = TRUE)"

```
## Citation

If you use this app or its outputs in academic work, please cite:

Frank Simmen (2026). Does Populism Cause Economic Outcomes?, Working Paper.

## License

This project is licensed under the MIT License. Free to use, modify, and distribute with attribution. No warranty provided.

© 2025 Frank Simmen
