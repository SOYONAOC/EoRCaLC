# EoRCaLC - Epoch of Reionization Calculator

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CUDA](https://img.shields.io/badge/CUDA-12.x-green.svg)](https://developer.nvidia.com/cuda-toolkit)

**EoRCaLC** is a Python package for calculating ionization diffusion and power spectrum analysis during the Epoch of Reionization (EoR). It features GPU acceleration support for high-performance cosmological simulations.

## Features

- 🚀 **GPU Acceleration**: High-performance calculations using CuPy (CUDA 12.x)
- 📊 **Power Spectrum Analysis**: Matter power spectrum calculations via CAMB
- 🌌 **Mass Functions**: Halo mass function computations with custom transfer functions
- 🔬 **Ionization Physics**: Ionization diffusion and recombination modeling
- ⚡ **Optimized Algorithms**: Efficient numerical integration and interpolation

## Installation

### Prerequisites

- Python 3.8 or higher
- CUDA 12.x compatible GPU (for GPU acceleration)
- CUDA Toolkit 12.x installed

### Basic Installation

```bash
pip install -e .
```

### GPU Support (CUDA 12.x)

For GPU acceleration with CUDA 12.x:

```bash
pip install -e ".[gpu]"
```

Or install CuPy separately:

```bash
pip install cupy-cuda12x
```

### Development Installation

```bash
pip install -e ".[dev]"
```

## Dependencies

### Core Dependencies
- **numpy** (≥1.20.0): Numerical computing
- **scipy** (≥1.6.0): Scientific computing and integration
- **astropy** (≥4.0): Astronomical calculations and units
- **camb** (≥1.3.0): Cosmological matter power spectrum
- **massfunc** (≥0.0.10): Halo mass function calculations
- **cupy-cuda12x** (≥12.0.0): GPU acceleration (optional)
- **filelock** (≥3.0.0): File locking for concurrent access
- **cosfunc**: Cosmological functions (local dependency)

## Usage

### Basic Example - CPU Version

```python
from eorcalc.ioninti import Ion

# Initialize ionization calculator
z = 7.0  # Redshift
ion_calc = Ion(
    z=z,
    fesc=0.2,          # Escape fraction
    A2byA1=0.1,        # Transfer function parameter
    ktrans=200,        # Transition scale [Mpc^-1]
    alpha=2.0,         # Power law index
    beta=0.0           # Secondary index
)

# Calculate ionization properties
# ... your calculations here ...
```

### GPU-Accelerated Version

```python
from eorcalc.ioninti_gpu import Ion

# Same interface but with GPU acceleration
ion_calc = Ion(
    z=7.0,
    fesc=0.2,
    A2byA1=0.1,
    ktrans=200,
    alpha=2.0,
    beta=0.0
)

# GPU-accelerated computations
# ... your calculations here ...
```

### Ionization Diffusion

```python
from eorcalc.iondiff import Ion

# Full ionization diffusion calculations
ion_diff = Ion(
    z=7.0,
    fesc=0.2,
    kakaka=0.7e-28,     # Recombination coefficient
    xi_ion=10**25.6,    # Ionizing efficiency
    A2byA1=0.1,
    ktrans=200,
    alpha=2.0,
    beta=0.0
)

# Perform diffusion calculations
# ... your calculations here ...
```

### Power Spectrum Analysis

```python
from eorcalc.powerspec import MassFunctions

# Initialize cosmological calculator
cosmo = MassFunctions(
    A2byA1=0.1,
    kMpc_trans=200,
    alpha=2.0,
    beta=0.0
)

# Calculate mass functions and power spectra
z = 7.0
M = 1e10  # Solar masses
sigma = cosmo.sigma_M(M, z)
mass_func = cosmo.massfunc(M, z)
```

## Module Overview

### `eorcalc.ioninti`
CPU-based ionization calculations with standard numerical methods.

### `eorcalc.ioninti_gpu`
GPU-accelerated ionization calculations using CuPy for enhanced performance.

### `eorcalc.iondiff`
Ionization diffusion modeling including recombination effects.

### `eorcalc.powerspec`
CAMB-based matter power spectrum and halo mass function calculations.

### `eorcalc.special`
Special functions and utilities:
- `qion_sb99`: Ionizing photon production from Starburst99 models
- `interp1d_gpu`: GPU-accelerated interpolation
- `xim`: Ionization fraction calculations
- `fstar`: Star formation efficiency

## Project Structure

```
EoRCaLC/
├── eorcalc/
│   ├── __init__.py
│   ├── iondiff.py       # Ionization diffusion (CPU)
│   ├── ioninti.py       # Ionization initial (CPU)
│   ├── ioninti_gpu.py   # Ionization (GPU accelerated)
│   ├── powerspec.py     # Power spectrum & mass functions
│   └── special.py       # Special functions & utilities
├── pyproject.toml       # Project configuration
├── LICENSE             # MIT License
└── README.md           # This file
```

## Performance Tips

1. **GPU Memory**: Ensure sufficient GPU memory for large-scale simulations
2. **Batch Processing**: Process multiple redshifts in batches for efficiency
3. **Caching**: Enable CAMB result caching for repeated calculations
4. **Precision**: Use `float32` on GPU for memory-intensive tasks if precision allows

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{eorcalc2025,
  author = {Hajime Hinata},
  title = {EoRCaLC: Epoch of Reionization Calculator},
  year = {2025},
  url = {https://github.com/SOYONAOC/IonDiff}
}
```

## Contact

- **Author**: Hajime Hinata
- **Email**: onmyojiflow@gmail.com
- **Repository**: https://github.com/SOYONAOC/IonDiff
- **Issues**: https://github.com/SOYONAOC/IonDiff/issues

## Acknowledgments

This package builds upon:
- **CAMB**: Cosmological Boltzmann code
- **CuPy**: NumPy-compatible array library for GPU
- **Astropy**: Community Python library for Astronomy
- **massfunc**: Halo mass function library

## Version History

### v0.1.0 (2025-10-02)
- Initial release
- CPU and GPU implementations
- Power spectrum analysis
- Ionization diffusion modeling
