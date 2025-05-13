% --- Adaptive Noise Cancellation Project ---

% Clear workspace and command window
clear; clc; close all;

% --- 1. Setup and Load Clean Audio ---
disp('1. Loading Clean Audio...');
audioFile = 'clean_speech.wav';

% Check if the user has a specific file, otherwise use a MATLAB default
if exist(audioFile, 'file')
    [clean_signal, Fs] = audioread(audioFile);
else
    disp(['File "', audioFile, '" not found. Using MATLAB built-in "handel.mat" instead.']);
    load handel.mat; % Contains 'y' (signal) and 'Fs' (sampling rate)
    clean_signal = y;
    % 'handel' is quite long, let's take a segment
    if length(clean_signal) > 5*Fs
        clean_signal = clean_signal(1:5*Fs);
    end
end

% Ensure the signal is a column vector
if size(clean_signal, 2) > 1
    clean_signal = mean(clean_signal, 2); % Convert to mono if stereo
end

% --- MODIFICATION FOR HEADROOM ---
target_clean_peak = 0.5; % e.g., clean signal peaks at +/- 0.5, leaving 0.5 for noise
clean_signal = clean_signal / max(abs(clean_signal)) * target_clean_peak;
% --- END MODIFICATION ---

%clean_signal = clean_signal / max(abs(clean_signal)); % Normalize to [-1, 1]
N = length(clean_signal); % Signal length
t = (0:N-1)' / Fs;      % Time vector

disp(['   Loaded audio with Fs = ', num2str(Fs), ' Hz, Duration = ', num2str(N/Fs), 's']);

% --- 2. Generate Noise Source ---
disp('2. Generating Noise Source...');
% This is the original source of the noise (e.g., an engine, a fan)
noise_power_relative_to_signal = 0.2; % Adjust to make noise more or less prominent
noise_source = randn(N, 1) * sqrt(noise_power_relative_to_signal * var(clean_signal) / var(randn(N,1)));

% --- 3. Create Noisy (Primary) and Reference Signals ---
disp('3. Creating Noisy and Reference Signals...');

% Simulate acoustic paths from the noise source to the microphones

% Path 1: Noise source to the primary microphone (where speech + noise is picked up)
% Let's model this as a simple FIR filter (e.g., a delay and some coloration)
path1_coeffs = fir1(31, 0.5); % Example FIR filter
noise_in_primary = filter(path1_coeffs, 1, noise_source);
noise_in_primary = noise_in_primary / max(abs(noise_in_primary)) * max(abs(noise_source)); % Rescale



% Path 2: Noise source to the reference microphone (picks up only noise, ideally)
% This path will be different from Path 1. The adaptive filter will try to model Path 1 given Path 2.
%path2_coeffs = fir1(41, 0.45); % e.g., slightly different cutoff, or different order
path2_coeffs = fir1(41, 0.5);
%path2_coeffs = fir1(35, 0.6); % Low-pass with a higher cutoff than path1

reference_noise = filter(path2_coeffs, 1, noise_source);
reference_noise = reference_noise / max(abs(reference_noise)) * max(abs(noise_source)); % Rescale

% Align signals using cross-correlation
[cc, lags] = xcorr(reference_noise, noise_in_primary, 50);
[~, max_idx] = max(abs(cc));
delay = lags(max_idx);
if(delay > 0)
    reference_noise = [reference_noise(delay+1:end); zeros(delay,1)];
else
    reference_noise = [zeros(-delay,1); reference_noise(1:end+delay)];
end
disp(['   Applied delay correction of ', num2str(delay), ' samples']);

% Create the noisy signal (primary input for the adaptive filter)
noisy_signal = clean_signal + noise_in_primary;
% Ensure noisy_signal doesn't clip excessively, though some minor clipping is okay
if max(abs(noisy_signal)) > 1
    noisy_signal = noisy_signal / max(abs(noisy_signal));
    disp('   Warning: Noisy signal was clipped after adding noise. Consider reducing noise_power_relative_to_signal.');
end

% PSD Plot 

% figure;
% subplot(2,1,1);
% [P_noise_in_primary, F_primary] = pwelch(noise_in_primary, hamming(512), 256, 512, Fs);
% plot(F_primary, 10*log10(P_noise_in_primary));
% title('PSD of Noise in Primary Mic Signal');
% xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
% grid on;
% 
% subplot(2,1,2);
% [P_reference_noise, F_ref] = pwelch(reference_noise, hamming(512), 256, 512, Fs);
% plot(F_ref, 10*log10(P_reference_noise));
% title('PSD of Reference Mic Signal');
% xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
% grid on;
% sgtitle('Power Spectral Densities of Noise Components');

% --- 4. Implement Adaptive Filter (LMS) ---
disp('4. Implementing Adaptive Filter (LMS)...');

% Prompt user for filter length
default_length = 48; % Default value
prompt = sprintf('Enter the LMS filter length (number of taps, default = %d): ', default_length);
filter_length = input(prompt);

% Validate input
if isempty(filter_length) || ~isnumeric(filter_length) || filter_length <= 0 || mod(filter_length, 1) ~= 0
    filter_length = default_length;
    disp(['   Invalid input. Using default filter length: ', num2str(filter_length)]);
else
    disp(['   Using filter length: ', num2str(filter_length)]);
end

% --- Using Normalized LMS (NLMS) ---
% The 'StepSize' for NLMS is typically between 0 and 2.
% Common values are in the range 0.01 to 1.0.
% NLMS is less sensitive to the exact step size value than standard LMS.
nlms_step_size = 0.05; % A good starting point for NLMS. Try 0.05, 0.2 etc.
                      % This step_size is less dependent on input signal power.

lms_filter = dsp.LMSFilter('Length', filter_length, ...
    'Method', 'Normalized LMS', ... % Key change here!
    'StepSize', nlms_step_size, ...
    'LeakageFactor',0.999);

% Apply the adaptive filter
% Inputs:
%   desired_signal (noisy_signal): The signal we want to clean.
%   input_signal (reference_noise): The noise reference.
% Outputs:
%   output_signal (estimated_noise): The filter's estimate of the noise in desired_signal.
%   error_signal (cleaned_signal): desired_signal - output_signal. This is our result.
[estimated_noise, cleaned_signal] = lms_filter(reference_noise, noisy_signal);



% For older MATLAB versions (if dsp.LMSFilter is not available):
% adaptfilt_lms = adaptfilt.lms(filter_length, step_size);
% [estimated_noise, cleaned_signal] = filter(adaptfilt_lms, reference_noise, noisy_signal);

disp('   Adaptive filtering complete.');

% --- 5. Analyze Results ---
disp('5. Analyzing Results...');

% a) Time-domain plots
figure;
subplot(3,1,1);
plot(t, clean_signal);
title('Clean Signal');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, noisy_signal);
title('Noisy Signal (Primary Input)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, cleaned_signal);
title('Cleaned Signal (Filter Output Error)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;
sgtitle('Adaptive Noise Cancellation Results (Time Domain)');

% b) Spectrograms
figure;
win_length = 256; % Window length for STFT
overlap = round(win_length * 0.75); % 75% overlap
nfft = 512; % FFT points

subplot(3,1,1);
spectrogram(clean_signal, hamming(win_length), overlap, nfft, Fs, 'yaxis');
title('Spectrogram of Clean Signal');

subplot(3,1,2);
spectrogram(noisy_signal, hamming(win_length), overlap, nfft, Fs, 'yaxis');
title('Spectrogram of Noisy Signal');

subplot(3,1,3);
spectrogram(cleaned_signal, hamming(win_length), overlap, nfft, Fs, 'yaxis');
title('Spectrogram of Cleaned Signal');
sgtitle('Adaptive Noise Cancellation Results (Frequency Domain)');



% c) Signal-to-Noise Ratio (SNR)
% SNR = P_signal / P_noise = var(signal) / var(noise)
% Note: MATLAB's snr() function calculates power differently (sum of squares).
% For a more traditional variance-based SNR:
% SNR_dB = 10 * log10(var(signal) / var(noise_component))

% Calculate the actual noise component in the noisy signal
actual_noise_component_in_noisy = noisy_signal - clean_signal;

% Calculate the residual noise component in the cleaned signal
residual_noise_component = cleaned_signal - clean_signal;
% Alternative for residual noise if clean_signal is not perfectly known:
% residual_noise_component = estimated_noise - noise_in_primary; (how well it matched)

% SNR of the noisy signal
power_clean_signal = var(clean_signal);
power_noise_in_noisy = var(actual_noise_component_in_noisy);
if power_noise_in_noisy == 0 % Avoid division by zero if no noise was added
    snr_noisy_db = inf;
else
    snr_noisy_db = 10 * log10(power_clean_signal / power_noise_in_noisy);
end

% SNR of the cleaned signal
power_residual_noise = var(residual_noise_component);
if power_residual_noise == 0 % Avoid division by zero if noise perfectly cancelled
    snr_cleaned_db = inf;
else
    snr_cleaned_db = 10 * log10(power_clean_signal / power_residual_noise);
end

disp(' ');
disp('--- Performance Evaluation ---');
fprintf('SNR of Noisy Signal: %.2f dB\n', snr_noisy_db);
fprintf('SNR of Cleaned Signal: %.2f dB\n', snr_cleaned_db);
fprintf('SNR Improvement: %.2f dB\n', snr_cleaned_db - snr_noisy_db);

% Optional: Listen to the audio
% sound(clean_signal, Fs); pause(length(clean_signal)/Fs + 0.5);
% sound(noisy_signal, Fs); pause(length(noisy_signal)/Fs + 0.5);
% sound(cleaned_signal, Fs);

disp(' ');
disp('--- End of Project ---');
