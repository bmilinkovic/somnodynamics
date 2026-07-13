# py_ssdi: Python Implementation of Dynamical Independence

This package implements Dynamical Independence computation for linear state-space systems, as described in:

L. Barnett and A. K. Seth, *Dynamical independence: discovering emergent macroscopic processes in complex dynamical systems*, [arXiv:2106.06511 [nlin.AO]](https://arxiv.org/abs/2106.06511), 2021.

Please cite the above if you use this code in your work.

This is a Python implementation of the original MATLAB [SSDI-1 toolbox](https://github.com/lcbarnett/ssdi).

## Installation

### Using Conda (recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/py_ssdi.git
cd py_ssdi

# Create and activate the conda environment
conda env create -f environment.yml
conda activate py_ssdi

# Install the package in development mode
pip install -e .
```

### Using pip

```bash
pip install py_ssdi
```

## Package Structure

The package is organized into the following modules:

- `py_ssdi.models`: Contains the implementation of state-space and VAR models
  - `StateSpaceModel`: Class for state-space model representation and manipulation
  - `VARModel`: Class for vector autoregressive model representation and manipulation
- `py_ssdi.metrics`: Contains implementations of various dynamical independence metrics
  - `dynamical_dependence`: Function to compute dynamical dependence
  - `causal_emergence`: Function to compute causal emergence
- `py_ssdi.utils`: Utility functions for optimization, transformation, and numerical stability
- `py_ssdi.visualization`: Tools for visualizing models, projections, and metrics
- `py_ssdi.core`: Core mathematical functions and algorithms

## Usage Examples

### Creating and Analyzing Models

```python
import numpy as np
from py_ssdi.models import StateSpaceModel, VARModel
from py_ssdi.metrics import dynamical_dependence, causal_emergence

# Create a random VAR model
n = 5  # dimension
r = 3  # model order
var_model = VARModel.create_random(n, r, 0.9)  # spectral radius 0.9

# Convert VAR model to state-space form
ss_model = var_model.to_state_space()

# Create a random state-space model directly
ss_model2 = StateSpaceModel.create_random(n, 3*n, 0.9)

# Compute dynamical dependence for a projection
L = np.random.randn(n, 2)  # random projection to 2D
L, _ = np.linalg.qr(L)  # orthonormalize

dd = dynamical_dependence(ss_model, L)
ce = causal_emergence(ss_model, L)

print(f"Dynamical Dependence: {dd}")
print(f"Causal Emergence: {ce}")
```

### Optimizing Projections

```python
from py_ssdi.utils.optimization import optimize_projection

# Start with a random projection
initial_L = np.random.randn(n, 2)
initial_L, _ = np.linalg.qr(initial_L)

# Optimize the projection to maximize causal emergence
optimal_L, optimal_ce = optimize_projection(
    ss_model, 
    initial_L, 
    metric='causal_emergence',
    method='gradient_ascent',
    n_iterations=100
)

print(f"Initial Causal Emergence: {causal_emergence(ss_model, initial_L)}")
print(f"Optimized Causal Emergence: {optimal_ce}")
```

### Visualization

```python
from py_ssdi.visualization import plot_causal_network, plot_spectral_density

# Plot the causal network of the model
plot_causal_network(ss_model)

# Plot the spectral density
plot_spectral_density(ss_model)

# Compare original and projected models
from py_ssdi.models import project_model
projected_model = project_model(ss_model, optimal_L)

# Plot both models' spectral densities for comparison
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plot_spectral_density(ss_model, ax=plt.gca(), title="Original Model")

plt.subplot(1, 2, 2)
plot_spectral_density(projected_model, ax=plt.gca(), title="Projected Model")

plt.tight_layout()
plt.show()
```

## Features

- State-space model creation and analysis
- VAR model creation and analysis
- Dynamical Independence metrics calculation
- Causal Emergence calculation and optimization
- Model transformation and projection
- Visualization tools for causal networks and spectral properties
- Integrated testing framework

## Dependencies

- NumPy: For numerical computations
- SciPy: For linear algebra and optimization
- Matplotlib: For visualization
- Control: For state-space system analysis
- NetworkX: For network representation and algorithms
- Pandas: For data manipulation
- Scikit-learn: For additional optimization algorithms

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License
