# State-Space Chamber Modeling for Particle Transport

This repository contains the MATLAB codebase used for developing physics-based grey-box state-space (chamber) models to simulate particle transport in fusion reactors, specifically for the MAST-U and TCV tokamaks. 

The codebase facilitates the extraction of experimental data, the iterative fitting of model parameters (such as confinement times and ionization times), and the rigorous validation and sensitivity analysis of the resulting models.

## General Overview

The workflow is divided into three main stages, each handled by a dedicated script:
1. **Data Retrieval (`GetData.m`)**: Extracts and formats experimental data from the reactor database.
2. **Model Identification (`Fitting.m`)**: Fits the grey-box state-space models to the experimental data using batched optimization.
3. **Validation & Analysis (`ModelTesting.m`)**: Tests the fitted models against validation data and performs sensitivity analyses (Monte Carlo and Bode).

### Prerequisites
* MATLAB (with the System Identification Toolbox).
* **`mastu_rstbx`**: An external toolbox required for extracting MAST-U data. Ensure this is added to your MATLAB path.

---

## How to Run the Main Scripts

### 1. `GetData.m`
**Purpose**: Extracts raw experimental data (Valve inputs, Interferometry, D-alpha, Fast Ion Gauges) from the MAST-U or TCV databases and saves it as formatted `.mat` files for system identification.

**Usage**:
1. Open `GetData.m`.
2. Set the `shotnumber` (or list of shots) you wish to extract.
3. Run the script. The output will be saved in the `Data/` directory as `[shotnumber]measuredTest.mat`.

Note this script is designed to be run on the DIFFER Rekenserver.

### 2. `Fitting.m`
**Purpose**: The core script for identifying the model parameters. It minimizes the error between the experimental data and the model's simulated response using `greyest`.

**Usage**:
1. Open `Fitting.m`.
2. Under "Settings", specify the `data_filenames` you generated in Step 1.
3. Choose the chamber model class (e.g., `settings.chamber_model = FourChamberModelUniversal();`).
4. Configure the `batch_strategy`. You can fit parameters in distinct stages (e.g., fitting gain parameters first, then time constants) to avoid local minima.
5. Run the script. It will align and preprocess the data, perform the batched fitting, and save the resulting system (usually to `Temp/Fitted_Models.mat`).

### 3. `ModelTesting.m`
**Purpose**: Loads a previously fitted model and evaluates its performance on specific datasets. It generates time-domain comparisons, frequency-domain (Bode) plots, and performs sensitivity analyses.

**Usage**:
1. Open `ModelTesting.m`.
2. Specify the `data_filenames` to use for validation.
3. Load your fitted model (e.g., `load("Temp\Fitted_Models.mat")`).
4. Toggle the analysis flags as needed (e.g., `settings.do_monte_carlo = true`, `settings.do_sensitivity_bode = true`).
5. Run the script to generate performance reports and visualizations.

---

## Chamber Models: Usage and Modification

The physics of the particle transport is defined in object-oriented classes (e.g., `FourChamberModelUniversal.m`, `ThreeChamberModelUniversal.m`). These classes act as grey-box templates for the MATLAB System Identification Toolbox.

### Understanding the Class Structure
If you open `FourChamberModelUniversal.m`, you will see it inherits from the `handle` class and contains two main sections: `properties` and `methods`.

#### Properties
This is where you configure the model architecture and initial guesses:
* **Architecture Settings**: Flags like `midplane_injection`, `direct_input`, and `ionisation_splitting` alter how the A and B matrices are constructed.
* **`default_parameters`**: A cell array containing the name and initial guess for every physical parameter (e.g., `'Ionization time', 0.003`).
* **`default_margin_factors`**: Defines the upper and lower bounds for the optimizer. Setting a margin to `1` fixes that parameter so it will not be fitted.

#### Methods
The methods define the mathematics and data handling:
* **`getMatrices(...)`**: This is the core function called by the `idgrey` estimator. It takes the physical parameters as inputs and returns the state-space matrices: $A$, $B$, $C$, and $D$. 
* **`preProccesData(raw_data, settings)`**: Defines how raw signals are mapped to the model's inputs and outputs. It handles alignment, detrending, and normalizations specific to the chosen reactor and diagnostics.

### How to Create a New Chamber Model
If you need to model a new geometry (e.g., a 5-Chamber model), follow these steps:
1. **Duplicate an existing model**: Copy `FourChamberModelUniversal.m` and rename it to `FiveChamberModel.m`.
2. **Update Properties**: Change the `name` property. Add or remove variables in `default_parameters` and `default_margin_factors` to match your new physics. Ensure the order of parameters matches exactly.
3. **Update State-Space Matrices**: Modify the `getAMatrix`, `getBMatrix`, `getCMatrix`, and `getDMatrix` functions. Ensure the matrix dimensions match your new state count (e.g., 5x5 for the A matrix).
4. **Update `preProccesData`**: If your new model requires a different diagnostic as an output, map that signal from the `raw_data` to the `allignedData` array. 
5. **Implement**: In `Fitting.m`, simply change the initialization to `settings.chamber_model = FiveChamberModel();`.