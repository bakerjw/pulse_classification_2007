# Pulse Identification in Ground Motion Records

This repository contains the original Matlab implementation of the pulse identification algorithm described in:

**Baker, J. W. (2007). "Quantitative classification of near-fault ground motions using wavelet analysis." *Bulletin of the Seismological Society of America*, 97(5), 1486-1501.**

## Archival Notice

This code represents the original implementation from 2006-2007 and is archived here for reproducibility and historical reference. The code is preserved largely as it was at the time of publication to maintain consistency with the published results.

**Original development:** January 2006 - October 2006
**Matlab version:** Compatible with Matlab R14 (2004) and later versions with Wavelet Toolbox

## Citation

If you use this code in your research, please cite:

```bibtex
@article{baker2007quantitative,
  title={Quantitative classification of near-fault ground motions using wavelet analysis},
  author={Baker, Jack W},
  journal={Bulletin of the Seismological Society of America},
  volume={97},
  number={5},
  pages={1486--1501},
  year={2007},
  publisher={Seismological Society of America}
}
```

## Overview

This implementation uses continuous wavelet transform (CWT) analysis to:

1. **Identify velocity pulses** in ground motion records
2. **Extract pulse characteristics** including pulse period (Tp) and peak ground velocity (PGV)
3. **Classify records** as pulse-like or non-pulse-like using a logistic regression model
4. **Filter late-arriving pulses** that are less relevant for structural response

The algorithm uses Daubechies wavelets (4th order) to decompose velocity time histories and identify dominant pulse-like features characteristic of near-fault ground motions.

## Quick Start

### Prerequisites

- Matlab with Wavelet Toolbox installed
- Required Matlab functions: `cwt`, `dwt`, `wavefun`, `wavemngr`, `scal2frq`

### Running the Example

1. Ensure all files are in the same directory or add the directory to your Matlab path
2. Run the main script:

```matlab
main_analyze_ground_motion
```

This will:
- Load the example ground motion (`example_ground_motion.acc`)
- Perform wavelet analysis to identify and extract the velocity pulse
- Display classification results in the command window
- Generate a 3-panel plot showing:
  - Original velocity time history
  - Extracted pulse
  - Residual (original minus pulse)

### Expected Output

The script will output several diagnostic values:

```matlab
pulse_indicator = 1.0000    % Probability from logistic regression (0-1)
late = 0                    % Boolean: is pulse late-arriving?
PGV = 129.9                 % Peak ground velocity (cm/s)
is_pulse = 1                % Final classification (1=pulse-like, 0=non-pulse-like)
Tp = 3.3000                 % Pulse period (seconds)
```

Classification criteria:
- `is_pulse = 1` (pulse-like) when:
  - `pulse_indicator > 0.85` AND
  - `PGV > 30 cm/s` AND
  - `late = 0` (pulse is not late-arriving)

## Files

| File | Description |
|------|-------------|
| `main_analyze_ground_motion.m` | Main script demonstrating the analysis workflow |
| `fn_analyze_record.m` | Core function performing wavelet analysis and pulse classification |
| `fn_build_pulse.m` | Reconstructs the extracted pulse from wavelet coefficients |
| `fn_plot_pulse.m` | Plots original record, extracted pulse, and residual |
| `example_ground_motion.acc` | Sample ground motion acceleration record for testing |
| `Pulse_like_records.csv` | Dataset of 91 pulse-like records from the NGA West 1 database identified using this method |

## Algorithm Details

### Wavelet Analysis Process

1. **Continuous Wavelet Transform (CWT):**
   - Input: Velocity time history (integrated from acceleration)
   - Wavelet: Daubechies 4th order (`db4`)
   - Scale range: Corresponding to periods from 0.25s to 15s

2. **Pulse Extraction:**
   - Identifies maximum wavelet coefficient
   - Determines pulse scale (related to pulse period Tp)
   - Extracts up to 10 coefficients in the neighborhood of the pulse
   - Reconstructs pulse and computes residual

3. **Classification Metrics:**
   - **PGV ratio:** Ratio of residual PGV to original PGV
   - **Energy ratio:** Ratio of DWT energy in residual to original
   - **Pulse indicator:** Logistic regression combining both metrics
   - **Late arrival check:** Compares cumulative energy thresholds

### Key Parameters

Defined in `fn_analyze_record.m`:

```matlab
num_coefs = 10              % Number of wavelet coefficients to extract
wname = 'db4'               % Daubechies 4th order wavelet
Tp_min = 0.25               % Minimum pulse period (s)
Tp_max = 15                 % Maximum pulse period (s)
pulse_indicator_threshold = 0.85   % Classification threshold
PGV_threshold = 30          % Minimum PGV (cm/s) for pulse classification
```

Logistic regression coefficients (from training on manually classified records):

```matlab
pulse_indicator = 1 / (1 + exp(-(23.297 - 14.579*(PGV_resid/PGV)
                                 - 20.477*(E_resid/E_orig))))
```

## Input File Format

The example acceleration file (`example_ground_motion.acc`) follows a specific format:
- First 6 lines: Header information
- Next line contains time step (dt) after two 5-character fields
- Remaining lines: Acceleration values in g

The script integrates acceleration to velocity (cm/s) before analysis.

## Validation Dataset

`Pulse_like_records.csv` contains 91 ground motion records classified as pulse-like, including:
- Event name, year, and station
- Wavelet-based pulse period (Tp)
- Spectral velocity-based pulse period
- PGV, magnitude (Mw), and distance metrics
- Pulse indicator value

This dataset can be used to validate the implementation or compare with manual classifications.

## Modifications from Original

The code is largely unchanged from the 2006-2007 version, with only minimal modifications:
- Line ending normalization (CRLF → LF)
- Minor typo fix in comments
- Addition of this README and supporting files

## Known Limitations

- Requires Matlab Wavelet Toolbox
- File reading code is specific to the acceleration file format used
- Hardcoded classification parameters optimized for the training dataset
- Limited error handling for edge cases

## License

See LICENSE file for details.

## Contact

For questions about the original research or methodology, please refer to the published paper. For issues with this archived code, please open an issue on the GitHub repository.

---

**Archived by:** Jack Baker
**Last updated:** January 2026
