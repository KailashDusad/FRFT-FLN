# Fraction Order Functional Link Network (FRFT-FLN) Based Adaptive Filter

This repository contains MATLAB code for nonlinear system identification using the **Fraction Order Functional Link Network (FRFT-FLN)**. The FRFT-FLN algorithm is a novel approach that combines fractional-order functional expansions with adaptive filtering to model and identify nonlinear systems effectively.

---

## Repository Contents

### Code
- **`FRFT_FLN.m`**  
  - Implements the FRFT-FLN-based adaptive filter for nonlinear system identification.
  - Key features include:
    - Fractional-order functional expansions.
    - Adaptive updates for fractional order (`alpha`) and weights.
    - Noise-tolerant learning with Mean Squared Error (MSE) tracking.
    - Smoothing of the learning curve using a moving average filter.
   
 - **`Comparison.m`**
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

## Key Features of FRFT-FLN

1. **Fractional Order Functional Expansion Block (FEB):**
   - Expands the input buffer using sinusoidal and exponential terms parameterized by a fractional order (`alpha`).
   - Enables flexible and detailed representation of input-output relationships in nonlinear systems.

2. **Adaptive Fractional Order (`alpha`):**
   - Dynamically updated using a gradient-based rule to optimize system performance during learning.

3. **Robust Weight Update Rule:**
   - Employs gradient descent to update filter weights, ensuring fast convergence and low steady-state error.

4. **Noise Resilience:**
   - Handles noisy environments with an additive Gaussian noise model for input signals.

5. **Learning Curve Analysis:**
   - Tracks the MSE across iterations, demonstrating the algorithm's convergence behavior.

---

## Usage Instructions

### Prerequisites
- Signal Processing Toolbox (required for generating additive white Gaussian noise).

### Running the Code
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/FRFT-FLN.git
   cd FRFT-FLN
   
2. Open MATLAB & Run FRFT_FLN.m for individual FRFT-FLN evaluation:
    - Modify parameters such as the number of inputs, noise level, and learning rates as needed.
    - Outputs a learning curve plot showing MSE over iterations.
3. Run Comparison.m to evaluate and compare the three filters:
    - Outputs a combined plot comparing MSE values of FRFT-FLN, AEFLN, and TFLN.

---

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

---
## Author

- **[Kailash Dusad]**  
- **[kailash.dusad@iitgn.ac.in]**  

