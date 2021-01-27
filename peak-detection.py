import pandas as pd
import csv
import os
from biosppy.signals import ecg


def inRange(calculated, actual, tolerance):
    # using comaparision operator
    if actual-5 <= num <= actual+5:
        return 1
    else:
        return 0


# Go through a loop for every file
for filename in os.listdir('A:/Dropbox (GaTech)\\School Misc\\BMED 3110 Lab\\ECG-peak-detection\\ECG Data'):
    open('A:/Dropbox (GaTech)\\School Misc\\BMED 3110 Lab\\ECG-peak-detection\\ECG Data' + str(filename), 'r'):
        # Read data from CSV datafile
        dataset = pd.read_csv("ECG-peak-detection\ECG Data\ECG Data 1a.csv")

        # Get sampling rate from the differnce between the first two timestamps
        hz = dataset["Time (s)"][1]-dataset["Time (s)"][0]
        # Get the amount of time for each sample based on the last timestamp
        seconds = dataset["Time (s)"][-1]

        # Process the ECG data using biosppy
        out = ecg.ecg(signal=dataset["Signal (V)"],
                      sampling_rate=hz, show=False)

        # Get the r_peaks from our processed data
        calculated_r_peaks = out[2]

        actual_r_peaks = []
        # Get the actual r-peaks from the original data
        for i in range(0, 21):  # The actual data is in the first 20 lines
            if dataset["Type"][i] == "N":
                actual_r_peaks.append(dataset["Sample #"][i])

        # Since we are dealing with messy data we need to define a tolerance for comparing our derived r peaks to the actual, 2 is a minimal tolerance
        comparison_tolerance = 2

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

        else:
            if len(calculated_r_peaks) - len(actual_r_peaks) > 0:
                FP = len(calculated_r_peaks) - len(actual_r_peaks)
                FN = 0
                TP = len(actual_r_peaks)
            else:
                FN = len(actual_r_peaks) - len(calculated_r_peaks)
                FP = 0
                TP = len(calculated_r_peaks)

        sensitivity = TP / (TP+FN)
        heartrate_bpm = (len(calculated_r_peaks)/seconds)*60

        file = open("ECG-peak-detection\\bpm_data\\" +
                    "results " + filename, 'w')
        file.writelines(["Calculated R Peaks:" + str(calculated_r_peaks), "Actual R Peaks:" + str(actual_r_peaks), "Sensitivity: " +
                         str(sensitivity), "FP: " + str(FP), "TP: " + str(TP), "FN: " + str(FN), "Heartrate(bpm): " + str(heartrate_bpm)])
        file.close()
