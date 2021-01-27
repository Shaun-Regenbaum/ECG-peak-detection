import pandas as pd
import numpy as np
import os
from biosppy.signals import ecg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# A function to check if calculated is within range of actual given a tolerance
def inRange(calculated, actual, tolerance):
    if actual-tolerance <= calculated <= actual+tolerance:
        return 1
    else:
        return 0


# A function to bring together all the standard visualization of ECG data in regards to heartrate.
def plot_ecg(ts=None,
             raw=None,
             filtered=None,
             rpeaks=None,
             templates_ts=None,
             templates=None,
             heart_rate_ts=None,
             heart_rate=None,
             path=None,
             show=False,
             case=None):

    MAJOR_LW = 2.5
    MINOR_LW = 1.5

    fig = plt.figure()
    fig.suptitle('Shaun Regenbaum:' + case)
    gs = gridspec.GridSpec(6, 2)

    # raw signal
    ax1 = fig.add_subplot(gs[:2, 0])
    ax1.plot(ts, raw, linewidth=MAJOR_LW, label='Raw Data')
    ax1.set_title("Raw Data")
    ax1.set_ylabel('Amplitude')
    ax1.legend()
    ax1.grid()

    # filtered signal with rpeaks
    ax2 = fig.add_subplot(gs[2:4, 0], sharex=ax1)
    ymin = np.min(filtered)
    ymax = np.max(filtered)
    alpha = 0.1 * (ymax - ymin)
    ymax += alpha
    ymin -= alpha
    ax2.plot(ts, filtered, linewidth=MAJOR_LW, label='Filtered Data')
    ax2.vlines(ts[rpeaks], ymin, ymax,
               color='m',
               linewidth=MINOR_LW,
               label='R-peaks')
    ax2.set_ylabel('Amplitude')
    ax2.set_title("Filtered Data")
    ax2.legend()
    ax2.grid()

    # heart rate
    ax3 = fig.add_subplot(gs[4:, 0], sharex=ax1)
    ax3.plot(heart_rate_ts, heart_rate, linewidth=MAJOR_LW, label='Heart Rate')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Heart Rate (bpm)')
    ax3.legend()
    ax3.set_title("Heart Rate")
    ax3.grid()

    # templates
    ax4 = fig.add_subplot(gs[1:5, 1])
    ax4.plot(templates_ts, templates.T, 'm', linewidth=MINOR_LW, alpha=0.7)
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Amplitude')
    ax4.set_title('Templates')
    ax4.grid()

    # make formatting a bit better
    plt.tight_layout()

    # save to file
    if path is not None:
        fig.savefig(path, dpi=200, bbox_inches='tight')

    # show
    if show:
        plt.show()
    else:
        # close
        plt.close(fig)


# My computer's path
path = 'A:/Dropbox (GaTech)\\School Misc\\BMED 3110 Lab\\ECG-peak-detection\\ECG Data'

# Go through a loop for every file
for filename in os.listdir(path):

    # Removing the file ext:
    case = os.path.splitext(filename)[0][-2:]
    print(case)

    # Read data from CSV datafile
    dataset = pd.read_csv('ECG Data/' + filename)

    # Get sampling rate from 1 over the differnce between the first two timestamps
    hz = 1 / (dataset["Time (s)"][1] - dataset["Time (s)"][0])
    # Get the amount of time for each sample based on the last timestamp
    seconds = dataset["Time (s)"][len(dataset["Time (s)"]) - 1]

    # Process the ECG data using biosppy
    out = ecg.ecg(signal=dataset["Signal (V)"],
                  sampling_rate=hz, show=False)

    # For plotting later
    ts, ts_tmpl, templates, ts_hr, hr = out[0], out[3], out[4], out[5], out[6]

    filtered_ecg = out[1]

    # Using Hamilton's Method
    hamilton = ecg.hamilton_segmenter(filtered_ecg, sampling_rate=hz)

    # Using Engelse and Zeelenberg's Method
    engzee = ecg.engzee_segmenter(filtered_ecg, sampling_rate=hz)

    # Get the r_peaks from our processed data
    calculated_r_peaks = hamilton[0]

    # Get the actual r-peaks from the original data to compare
    actual_r_peaks = []
    for i in range(0, 23):  # The actual data is in the first 20 lines
        if dataset["Type"][i] == "N":  # Add it if it is a R peak
            actual_r_peaks.append(dataset["Sample #"][i])

    # Since we are dealing with messy data we need to define a tolerance for comparing our derived r peaks to the actual, 5 is a minimal tolerance
    comparison_tolerance = 5

    # calculating TP, FP, FN, Sensitivty:
    if len(calculated_r_peaks) == len(actual_r_peaks):
        print("Got same number of r-peaks")
        TP = 0
        FP = 0
        FN = 0
        for i in range(len(calculated_r_peaks)):
            if inRange(calculated_r_peaks[i], actual_r_peaks[i], comparison_tolerance):
                TP = TP + 1
            else:
                FN = FN + 1
                FP = FP + 1

    else:  # I took some shortcuts assuming my algorithim was just going to either miss some or add some, but never both at the same time.
        print("Not the same number of r-peaks")
        if len(calculated_r_peaks) - len(actual_r_peaks) > 0:
            FP = len(calculated_r_peaks) - len(actual_r_peaks)
            FN = 0
            TP = len(actual_r_peaks)
        else:
            FN = len(actual_r_peaks) - len(calculated_r_peaks)
            FP = 0
            TP = len(calculated_r_peaks)

    # Calculating Sensitivty
    sensitivity = TP / (TP+FN)

    # Calculating bpm from our algorithim's r peaks
    calculated_heartrate_bpm = (len(calculated_r_peaks)/seconds)*60

    # Calculating bpm from the actual r peaks
    actual_heartrate_bpm = (len(actual_r_peaks)/seconds)*60

    # Writing this all to a file for easy viewing later
    file = open("results hamilton " + case + ".csv", 'w')
    file.writelines(["Calculated R Peaks:" + str(calculated_r_peaks) + '\n',
                     "Actual R Peaks:" + str(actual_r_peaks) + '\n',
                     "Sensitivity: " + str(sensitivity) + '\n',
                     "FP: " + str(FP) + '\n',
                     "TP: " + str(TP) + '\n',
                     "FN: " + str(FN) + '\n',
                     "Calculated Heartrate(bpm): " +
                     str(calculated_heartrate_bpm) + '\n',
                     "Actual Heartrate(bpm): " + str(actual_heartrate_bpm) + '\n'])

    # We can then plot it, and save it to a file
    plot_ecg(ts=ts,
             raw=dataset["Signal (V)"],
             filtered=filtered_ecg,
             rpeaks=calculated_r_peaks,
             templates_ts=ts_tmpl,
             templates=templates,
             heart_rate_ts=ts_hr,
             heart_rate=hr,
             path="graphics hamilton " + case + ".jpg",
             show=True,
             case=case)

    file.close()
