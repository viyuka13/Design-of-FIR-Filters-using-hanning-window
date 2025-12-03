# 4 C Design-of-FIR-Filters-using-hanning-window

#DESIGN OF FIR DIGITAL FILTER 

# AIM: 
          
  To generate design of low-pass and high-pass FIR digital filter using SCILAB 

# APPARATUS REQUIRED: 

  PC Installed with SCILAB 

# PROGRAM:

__LOW-PASS FIR FILTER__

```c
clc;
clear;
close;

// Low Pass FIR Filter using Hanning Window
N = 36;              // Filter length
fc = 0.3;            // Normalized cutoff frequency (fc = Fcutoff / (Fs/2))
n = 0:N-1;
alpha = (N-1)/2;

// Hanning window
w = 0.5 - 0.5*cos(2*%pi*n/(N-1));

// Ideal Low Pass Filter Impulse Response
hd = zeros(1, N);
for i = 1:N
    if (i-1) == alpha then
        hd(i) = 2*fc;
    else
        hd(i) = sin(2*%pi*fc*(i-1-alpha)) / (%pi*(i-1-alpha));
    end
end

// Multiply by window
h = hd .* w;

// Frequency Response
[H, f] = frmag(h, 512);

// Plot
figure;
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('LOW PASS FIR FILTER (Hanning Window)');

subplot(2,1,2);
plot(f, atan(imag(H), real(H)));
xlabel('Normalized Frequency');
ylabel('Phase (radians)');
title('Phase Response');

```

__HIGH-PASS FIR FILTER__

```c
clc;
clear;
close;

// High Pass FIR Filter using Hanning Window
N = 41;              // Filter length
fc = 0.1;            // Normalized cutoff frequency
n = 0:N-1;
alpha = (N-1)/2;

// Hanning window
w = 0.5 - 0.5*cos(2*%pi*n/(N-1));

// Ideal High Pass Filter Impulse Response
hd = zeros(1, N);
for i = 1:N
    if (i-1) == alpha then
        hd(i) = 1 - 2*fc;
    else
        hd(i) = -sin(2*%pi*fc*(i-1-alpha)) / (%pi*(i-1-alpha));
    end
end

// Multiply by window
h = hd .* w;

// Frequency Response
[H, f] = frmag(h, 512);

// Plot
figure;
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('HIGH PASS FIR FILTER (Hanning Window)');

subplot(2,1,2);
plot(f, atan(imag(H), real(H)));
xlabel('Normalized Frequency');
ylabel('Phase (radians)');
title('Phase Response');

```


# OUTPUT:

__LOW-PASS FIR FILTER__


<img width="809" height="411" alt="4L" src="https://github.com/user-attachments/assets/23aa70eb-2c20-4fdf-ab66-dd5b3cb75b7f" />


__HIGH-PASS FIR FILTER__


<img width="818" height="419" alt="4H" src="https://github.com/user-attachments/assets/69116754-182e-4de2-bca3-063eee63cde8" />


# RESULT

Design of low-pass and High-pass FIR digital filter using SCILAB is generated.
