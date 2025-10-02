# EoRCaLC - Epoch of Reionization Calculator

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CUDA](https://img.shields.io/badge/CUDA-12.x-green.svg)](https://developer.nvidia.com/cuda-toolkit)

**EoRCaLC** is a GPU-accelerated Python package for simulating reionization fields and calculating ionization evolution during the Epoch of Reionization (EoR). It combines high-performance GPU computing with sophisticated cosmological physics to model the ionization history of the universe.

## Features

- 🚀 **GPU Acceleration**: High-performance calculations using CuPy (CUDA 12.x)
- 🌌 **Reionization Field Simulation**: Full 3D ionization field evolution with recombination
- 📊 **Power Spectrum Analysis**: Matter power spectrum calculations via CAMB
- 🔬 **Ionization Physics**: Advanced ionization/recombination modeling with mini-halos
- ⚡ **Multi-scale Smoothing**: Efficient top-hat filtering at multiple scales
- 📈 **Optical Depth Calculation**: Automatic Thomson scattering optical depth computation
- 💾 **Data Management**: Binary field I/O and CSV output for analysis

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
pip install cupy-cuda12x
```

For CUDA 11.x users:

```bash
pip install -e ".[gpu-cuda11]"
```

### Development Installation

```bash
pip install -e ".[dev]"
```

## Dependencies

### Core Dependencies
- **numpy** (≥1.20.0): Numerical computing
- **scipy** (≥1.6.0): Scientific computing and integration
- **pandas** (≥1.3.0): Data analysis and CSV output
- **astropy** (≥4.0): Astronomical calculations and units
- **camb** (≥1.3.0): Cosmological matter power spectrum
- **massfunc** (≥0.0.10): Halo mass function calculations
- **cosfunc** (≥0.1.0): Cosmological functions (n_H, dtdz, etc.)
- **xxiop** (≥0.1.0): Optical depth calculations
- **cupy-cuda12x** (≥12.0.0): GPU acceleration
- **filelock** (≥3.0.0): File locking for concurrent access

## Module Overview

### `eorcalc.reion_field`
Main module for reionization field simulation with GPU acceleration.

**Main Function**: `reionization_calculator()`
- Simulates 3D ionization field evolution across redshift
- Includes recombination effects and mini-halo feedback
- Multi-scale smoothing with top-hat filters
- Calculates optical depth and saves results

### `eorcalc.ioninti_gpu`
GPU-accelerated ionization physics calculations.

**Class**: `Ion`
- Source term calculations (nion_interp, nion_st)
- Mini-halo absorption (nxi_interp, nxi_st)
- IGM recombination (dnrec_dz_path)
- GPU-optimized numerical integration

### `eorcalc.ioninti`
CPU-based ionization calculations (alternative to GPU version).

### `eorcalc.iondiff`
Ionization diffusion modeling with star formation rate density.

### `eorcalc.powerspec`
Power spectrum and mass function calculations.

**Class**: `MassFunctions`
- Matter power spectrum via CAMB
- Halo mass functions (Sheth-Tormen, Press-Schechler, EPS)
- Custom transfer functions with scale-dependent modifications

### `eorcalc.special`
Utility functions for data I/O and GPU interpolation.
- `load_binary_data()`: Load density fields from binary files
- `TopHat_filter()`: GPU-accelerated top-hat filtering
- `xHII_field_update()`: Ionization field updates
- `interp1d_gpu()`: GPU interpolation
- `fstar()`, `xim()`: Star formation efficiency and mini-halo feedback

## Usage

### Reionization Field Simulation

```python
from eorcalc import reionization_calculator

# Run full reionization simulation
optical_depth = reionization_calculator(
    fesc=0.2,              # Escape fraction
    A2byA1=1.0,            # Transfer function amplitude ratio
    kMpc_trans=1e6,        # Transition scale [Mpc^-1]
    alpha=0.0,             # Power law index
    beta=0.0,              # Secondary index
    label='MH',            # Output label
    DIM=256,               # Grid dimension
    box_length=800,        # Box size [Mpc]
    save_on=True          # Save ionization fields
)

print(f"Optical depth: {optical_depth:.4f}")
```

**Output**:
- Ionization fields saved to `reionf/{label}/rf_{z:.2f}.npy`
- Ionization fraction history saved to `csvfile/{label}.csv`
- Prints redshift evolution and optical depth

### GPU-Accelerated Ionization Physics

```python
from eorcalc.ioninti_gpu import Ion
import cupy as cp

# Initialize at redshift z=7
ion = Ion(
    z=7.0,
    fesc=0.2,
    A2byA1=1.0,
    ktrans=1e6,
    alpha=0.0,
    beta=0.0
)

# Load density field (GPU array)
delta_field = cp.random.randn(256, 256, 256)

# Calculate source term
m_grid = ion.cosmo.rhom * (800/256)**3
source = ion.nion_interp(m_grid, delta_field)

# Calculate mini-halo absorption
minihalo = ion.nxi_interp(m_grid, delta_field)

# IGM neutral hydrogen
igm = ion.n_HI(delta_field)

print(f"Mean source: {cp.mean(source):.2e}")
print(f"Mean mini-halo: {cp.mean(minihalo):.2e}")
```

### Power Spectrum and Mass Functions

```python
from eorcalc.powerspec import MassFunctions

# Initialize with custom transfer function
cosmo = MassFunctions(
    A2byA1=1.0,      # Amplitude ratio
    kMpc_trans=1e6,  # Transition scale
    alpha=0.0,       # Power index
    beta=0.0         # Secondary index
)

# Calculate halo mass function at z=7
z = 7.0
M = 1e10  # Solar masses

# Sheth-Tormen mass function
dndm_st = cosmo.dndmst(M, z)

# Press-Schechter mass function
dndm_ps = cosmo.dndmps(M, z)

# Excursion set peaks with environmental dependence
deltaV = 0.5
Mv = 1e12
dndm_eps = cosmo.dndmeps(M, Mv, deltaV, z)

print(f"ST dndm: {dndm_st:.2e} Mpc^-3")
print(f"PS dndm: {dndm_ps:.2e} Mpc^-3")
```

### Data I/O

```python
from eorcalc.special import load_binary_data, TopHat_filter
import cupy as cp

# Load density field from binary file
delta_field = load_binary_data(
    'ktrans1e2/updated_smoothed_deltax_z006.00_256_800Mpc',
    DIM=256
)

# Convert to GPU array
delta_gpu = cp.asarray(delta_field)

# FFT for filtering
delta_ffted = cp.fft.rfftn(delta_gpu, norm="forward")

# Apply top-hat filter at R=10 Mpc
delta_smoothed = TopHat_filter(
    delta_ffted, 
    R=10, 
    DIM=256, 
    box_length=800
)

print(f"Original variance: {cp.var(delta_gpu):.4f}")
print(f"Smoothed variance: {cp.var(delta_smoothed):.4f}")
```

## Physics Implementation

### Reionization Model

The package implements a semi-numerical reionization model:

1. **Source Term**: 
   $$n_{\rm ion} = f_{\rm esc} \cdot Q_{\rm ion} \cdot f_* \cdot \frac{\rho_{\rm b}}{\rho_{\rm m}} \cdot \int_{M_{\rm min}}^{M_{\rm max}} \frac{dn}{dM} \, M \, dM$$

2. **Mini-halo Absorption**:
   $$n_{\xi} = \frac{\rho_{\rm b}}{\rho_{\rm m}} \cdot \int_{M_{\rm J}}^{M_{\rm min}} \xi_{\rm ion}(M) \cdot \frac{dn}{dM} \cdot M \, dM$$

3. **Recombination**:
   $$\frac{dn_{\rm rec}}{dz} = C_{\rm HII} \cdot x_{\rm HE} \cdot \alpha_A \cdot n_{\rm HI} \cdot Q_{\rm HII} \cdot (1+z)^3 \cdot \frac{dt}{dz}$$

4. **Ionization Criterion**:
   $$\bar{f}_{\rm esc} \cdot n_{\rm ion}(R) > n_{\rm HI}(R) + \bar{n}_{\xi} \cdot n_{\xi}(R)$$

### Multi-scale Smoothing

The code performs filtering at multiple scales from cell size to 50 Mpc:
- Logarithmically spaced smoothing radii
- Top-hat filters in Fourier space
- Ionization propagates from large to small scales

## Directory Structure

```
EoRCaLC/
├── eorcalc/                 # Main package
│   ├── __init__.py         # Package initialization
│   ├── reion_field.py      # Main reionization simulator
│   ├── ioninti_gpu.py      # GPU ionization physics
│   ├── ioninti.py          # CPU ionization physics
│   ├── iondiff.py          # Ionization diffusion
│   ├── powerspec.py        # Power spectrum & mass functions
│   └── special.py          # Utility functions
├── pyproject.toml          # Package configuration
├── README.md               # This file
└── LICENSE                 # MIT License
```

## Input Data Format

The package expects density field data in binary format:
- Filename pattern: `updated_smoothed_deltax_z{z:06.2f}_{DIM}_{box_length}Mpc`
- Binary format: Little-endian float32
- Shape: (DIM, DIM, DIM)
- Values: Overdensity δ = ρ/ρ̄ - 1

## Output Data

### Ionization Fields
- Location: `reionf/{label}/rf_{z:.2f}.npy`
- Format: NumPy array (CuPy saved)
- Shape: (DIM, DIM, DIM)
- Values: Ionized fraction (0 to 1)

### Ionization History
- Location: `csvfile/{label}.csv`
- Columns: `z` (redshift), `ionf` (ionization fraction)
- Sorted by decreasing redshift

## Performance Notes

- **GPU Memory**: Requires ~2-4 GB VRAM for 256³ grids
- **Speed**: ~10-100x faster than CPU for large grids
- **Multi-GPU**: Not currently supported
- **Precision**: Uses float32 for GPU, float64 for critical calculations

## Examples

See the Jupyter notebook `cece.ipynb` for complete examples and visualization.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{eorcalc,
  author = {Hajime Hinata},
  title = {EoRCaLC: GPU-Accelerated Reionization Calculator},
  year = {2025},
  url = {https://github.com/SOYONAOC/EoRCaLC}
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- CAMB for cosmological calculations
- CuPy team for GPU acceleration tools
- The reionization modeling community

## Contact

- Author: Hajime Hinata
- Email: onmyojiflow@gmail.com
- GitHub: [SOYONAOC/EoRCaLC](https://github.com/SOYONAOC/EoRCaLC)

## Troubleshooting

### CUDA Issues
```bash
# Check CUDA version
nvcc --version

# Install matching CuPy version
pip install cupy-cuda11x  # for CUDA 11.x
pip install cupy-cuda12x  # for CUDA 12.x
```

### Import Errors
```bash
# Install missing dependencies
pip install cosfunc xxiop massfunc

# Or install all dependencies
pip install -e .
```

### Memory Errors
- Reduce grid dimension (DIM)
- Process fewer redshift snapshots at once
- Use CPU version for very large simulations
