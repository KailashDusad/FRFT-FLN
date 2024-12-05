# Nonlinear System Identification using Adaptive Filters

This repository contains MATLAB implementations and results for nonlinear system identification using various adaptive filter algorithms, including:

- **Fractional Fourier Transform Functional Link Network (FRFT-FLN)**
- **Trigonometric Functional Link Network (TFLN)**
- **Adaptive Exponential Functional Link Network (AEFLN)**

The repository includes the following:

- MATLAB code for FRFT-FLN-based nonlinear system identification.
- MATLAB code comparing FRFT-FLN, AEFLN, and TFLN for various nonlinear systems.
- Result images showcasing the performance of these adaptive filters.

---

## Repository Contents

### Code Files
1. **`FRFT_FLN.m`**
   - Implements the FRFT-FLN-based adaptive filter for nonlinear system identification.
   - Key features:
     - Functional expansion using fractional order.
     - Adaptive fractional order (`alpha`) update mechanism.
     - Handles multiple independent trials for ergodicity.
   - Generates Mean Squared Error (MSE) plots for learning curves.

2. **`Comparison.m`**
   - Compares the performance of FRFT-FLN, AEFLN, and TFLN adaptive filters.
   - Key features:
     - Implements each algorithm with the respective functional expansion.
     - Evaluates performance across different nonlinear systems.
     - Produces MSE plots for all algorithms.

---

### Results
The repository includes result images demonstrating the performance of the algorithms:
  - **Learning Curve (MSE) of FRFT-FLN**
  - **Comparison of MSE for FRFT-FLN, AEFLN, and TFLN**

---

## Usage Instructions

### Prerequisites
   - Signal Processing Toolbox (for generating additive white Gaussian noise).

### Running the Code
1. Clone the repository:
   ```bash
   git clone https://github.com/KailashDusad/FRFT-FLN.git
   cd FRFT-FLN
2. Open MATLAB & Run FRFT_FLN.m for individual FRFT-FLN evaluation:
    - Modify parameters such as the number of inputs, noise level, and learning rates as needed.
    - Outputs a learning curve plot showing MSE over iterations.
3. Run Comparison.m to evaluate and compare the three filters:
    - Outputs a combined plot comparing MSE values of FRFT-FLN, AEFLN, and TFLN.



# Key Algorithms and Concepts
1. Fractional Fourier Transform Functional Link Network (FRFT-FLN)

    - Incorporates fractional order parameters in functional expansion.
    - Adaptive updates for fractional order and weights.
    - Enhances modeling of complex nonlinear systems.

2. Trigonometric Functional Link Network (TFLN)

    - Uses trigonometric functional expansions (sin and cos).
    - Simpler approach for nonlinear system modeling.

3. Adaptive Exponential Functional Link Network (AEFLN)

    - Employs exponential functional expansions with an adaptive parameter.
    - Suitable for systems with varying degrees of nonlinearity.

# Results Overview
 ## FRFT-FLN Learning Curve:
   - Demonstrates faster convergence and lower steady-state error compared to TFLN and AEFLN.
## Comparison Results:
  - FRFT-FLN outperforms AEFLN and TFLN in terms of MSE for most nonlinear systems.
