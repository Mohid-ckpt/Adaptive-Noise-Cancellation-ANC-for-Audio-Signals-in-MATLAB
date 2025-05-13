# Adaptive Noise Cancellation (ANC) for Audio Signals in MATLAB

This project implements an Adaptive Noise Cancellation (ANC) system using MATLAB and the DSP System Toolbox. The goal is to remove background noise from an audio signal, simulating real-time noise cancellation scenarios often found in telecommunication headsets or noisy environments.

## Table of Contents
1.  [Project Description](#project-description)
2.  [Core Concepts](#core-concepts)
3.  [System Diagram](#system-diagram)
4.  [Features Implemented](#features-implemented)
5.  [MATLAB Features Used](#matlab-features-used)
6.  [Setup and Installation](#setup-and-installation)
7.  [How to Run](#how-to-run)
8.  [Expected Output](#expected-output)
9.  [Parameters and Customization](#parameters-and-customization)
10. [Results and Evaluation](#results-and-evaluation)
11. [Future Work / Potential Extensions](#future-work--potential-extensions)
12. [License](#license)

## Project Description
The project demonstrates the application of adaptive filtering techniques, specifically the Least Mean Squares (LMS) algorithm, to reduce additive noise from a desired audio signal. A simulated environment is created where:
*   A clean audio signal (e.g., speech) is the desired signal.
*   A noise source is generated.
*   This noise source corrupts the clean signal via one acoustic path (Primary Path).
*   The same noise source is also picked up by a reference sensor via a different acoustic path (Reference Path).
The adaptive filter uses the reference noise to estimate and subtract the noise component from the corrupted primary signal.

## Core Concepts
Adaptive Noise Cancellation (ANC) works on the principle of using a reference input that is correlated with the noise in the primary input but uncorrelated with the desired signal.

1.  **Primary Input (`d`):** Desired signal + Noise (e.g., `clean_speech + noise_in_primary`)
2.  **Reference Input (`x`):** Noise correlated with the noise in `d` (e.g., `reference_noise`)
3.  **Adaptive Filter:** Takes `x` as input and tries to produce an estimate `y` of the noise component in `d`.
4.  **Error Signal (`e`):** `e = d - y`. This is the filter's output and, ideally, the cleaned signal. The filter adapts its coefficients to minimize the power of `e`.

## System Diagram:

![Basic-ANC-system-structure](https://github.com/user-attachments/assets/f0c8174c-8af5-4e66-a9c4-5834a72ffd32)

## Features Implemented
*   Loading a clean audio signal.
*   Generation of a synthetic noise source (`randn`).
*   Simulation of different acoustic paths for the noise affecting the primary and reference inputs using FIR filters.
*   Creation of a noisy signal by adding filtered noise to the clean signal.
*   Implementation of an LMS adaptive filter using `dsp.LMSFilter`.
*   Visualization of signals in time and frequency domains (spectrograms).
*   Performance evaluation using Signal-to-Noise Ratio (SNR) improvement.

## MATLAB Features Used
*   **DSP System Toolbox:**
    *   `dsp.LMSFilter` (for implementing the LMS adaptive filter).
    *   (Alternatively, older `adaptfilt.lms` or `adaptfilt.rls` could be used, but `dsp.LMSFilter` is current).
*   **Signal Processing Toolbox:**
    *   `audioread` (for loading audio files).
    *   `filter` (for simulating acoustic paths).
    *   `fir1` (for designing FIR filters for acoustic paths).
    *   `spectrogram` (for frequency-domain visualization).
    *   `var` (for variance calculation in SNR).
*   **Core MATLAB:**
    *   `randn` (for generating Gaussian white noise).
    *   `plot` (for time-domain visualization).
    *   Basic arithmetic and array operations.

## Setup and Installation
1.  **MATLAB:** Ensure you have MATLAB installed (preferably R2016b or newer for `dsp.LMSFilter`).
2.  **Toolboxes:**
    *   DSP System Toolbox
    *   Signal Processing Toolbox
3.  **Audio File:**
    *   The script is set up to look for `clean_speech.wav`. Place your desired clean audio file with this name in the project directory.
    *   If `clean_speech.wav` is not found, the script defaults to using MATLAB's built-in `handel.mat` audio.
4.  **Clone the Repository (Optional):**
    ```bash
    git clone <your-repository-url>
    cd <your-repository-name>
    ```

## How to Run
1.  Open MATLAB.
2.  Navigate to the project directory where `anc_project.m` (or your main script file) is located.
3.  Run the script from the MATLAB Command Window:
    ```matlab
    >> anc_project
    ```
4.  The script will:
    *   Load/generate the clean audio.
    *   Generate noise and create noisy/reference signals.
    *   Apply the LMS adaptive filter.
    *   Display plots of the clean, noisy, and cleaned signals (time-domain and spectrograms).
    *   Print SNR values (original noisy SNR, cleaned SNR, and SNR improvement) to the Command Window.

## Expected Output
The script will produce:
1.  **Command Window Output:**
    *   Status messages during execution.
    *   Calculated SNR for the noisy signal.
    *   Calculated SNR for the cleaned signal.
    *   SNR improvement in dB.
2.  **Figure 1: Time Domain Plots:**
    *   Subplot 1: Clean Signal
    *   Subplot 2: Noisy Signal
    *   Subplot 3: Cleaned Signal (output of the ANC)
3.  **Figure 2: Spectrograms:**
    *   Subplot 1: Spectrogram of Clean Signal
    *   Subplot 2: Spectrogram of Noisy Signal
    *   Subplot 3: Spectrogram of Cleaned Signal

**Example Placeholder for Plots:**

**(Imagine a screenshot here showing the time-domain plots: clean, noisy, and clearly less noisy output signal)**

**(Imagine a screenshot here showing the spectrograms: clean, noisy with a visible noise floor, and cleaned with a reduced noise floor)**

## Parameters and Customization
The `anc_project.m` script has several key parameters than can be adjusted to observe different behaviors:

*   `audioFile`: Path to the clean audio input.
*   `noise_power_relative_to_signal`: Adjusts the power of the generated noise relative to the clean signal. (e.g., `0.1` for low noise, `0.5` for moderate, `1.0` for high noise).
*   `path1_coeffs`, `path2_coeffs`: FIR filter coefficients simulating the acoustic paths. Changing these alters how the noise source is perceived by the primary and reference inputs.
*   `filter_length` (for LMS filter): The number of taps in the adaptive filter. Should be sufficient to model the impulse response of `path1_coeffs` (or the unknown primary noise path).
*   `step_size` (LMS `mu`): Crucial for convergence.
    *   Too small: Slow convergence.
    *   Too large: Instability and divergence.
    *   Typical range depends on input signal power; for `dsp.LMSFilter` with `Method = 'LMS'`, it's often small (e.g., 0.001 - 0.05). For `'Normalized LMS'`, it's typically between 0 and 2 (e.g., 0.1 - 1.0).
*   `lms_filter.Method`: Can be set to `'LMS'`, `'Normalized LMS'`, etc., offering different adaptation characteristics.

Experimenting with these parameters is highly encouraged to understand their impact on ANC performance.

## Results and Evaluation
The primary metric for evaluation is the **Signal-to-Noise Ratio (SNR) improvement**.
The script calculates:
*   **SNR of Noisy Signal:** `10 * log10(var(clean_signal) / var(actual_noise_component_in_noisy))`
*   **SNR of Cleaned Signal:** `10 * log10(var(clean_signal) / var(residual_noise_component))`
*   **SNR Improvement:** `SNR_cleaned_dB - SNR_noisy_dB`

A positive SNR improvement indicates successful noise reduction. Visual inspection of the time-domain plots and spectrograms also provides qualitative assessment of the noise cancellation.

**Typical Observations:**
*   The cleaned signal's waveform should more closely resemble the original clean signal compared to the noisy signal.
*   The spectrogram of the cleaned signal should show a reduction in the noise floor, particularly in frequency bands where the noise was prominent.

## Future Work / Potential Extensions
*   **RLS Algorithm:** Implement and compare with the Recursive Least Squares (RLS) algorithm (`dsp.RLSFilter`), which generally offers faster convergence at the cost of higher computational complexity.
*   **Normalized LMS (NLMS):** Explore the NLMS variant (available within `dsp.LMSFilter` by setting `Method = 'Normalized LMS'`) for improved stability and convergence speed, especially with varying input signal power.
*   **Different Noise Types:** Experiment with various noise types (e.g., colored noise, recorded real-world noise) instead of just `randn`.
*   **Real-Time Implementation:** Adapt the script for real-time audio processing using `audioDeviceReader` and `audioDeviceWriter` from the Audio Toolbox.
*   **Multi-Channel ANC:** Extend to scenarios with multiple noise sources and/or multiple reference microphones.
*   **Varying Acoustic Paths:** Simulate time-varying acoustic paths to test the tracking capability of the adaptive filter.
*   **Voice Activity Detection (VAD):** Integrate a VAD to pause filter adaptation during periods of silence in the desired signal, which can prevent divergence or undesirable adaptation to the clean signal itself.

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.

---
This README provides a solid foundation. Remember to replace placeholders (like `<your-repository-url>`) and add your own screenshots of the results to make it compelling for your portfolio!
