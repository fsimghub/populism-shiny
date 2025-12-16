# Shiny Population Data App

A Shiny application for exploring and visualizing population data.

## Data

The app uses a pre-processed dataset (`ShinyPopDat.qs`, ~93 MB) that is too large to store directly in the GitHub repository.

The data file is hosted as a release asset:

- Download link: [https://github.com/fsimghub/populism-shiny/releases/download/v1.0/ShinyPopDat.qs](https://github.com/fsimghub/populism-shiny/releases/download/v1.0-data/ShinyPopDat.qs)

The app automatically downloads the file on first run if it is not already present in the `data/` directory.

## Installation & Running Locally

```r
# Clone the repository
git clone https://github.com/fsimghub/populism-shiny.git
cd populism-shiny

# Install dependencies (if needed)
install.packages(c("shiny", "qs", ...))  # add any other required packages

# Run the app
shiny::runApp()
```

## License

This project is licensed under the MIT License.
