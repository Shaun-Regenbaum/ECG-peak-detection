import pandas as pd
import csv
import os
from biosppy.signals import ecg


def inRange(calculated, actual, tolerance):
    # using comaparision operator
    if actual-tolerance <= calculated <= actual+tolerance:
        return 1
    else:
        return 0


path = 'A:/Dropbox (GaTech)\\School Misc\\BMED 3110 Lab\\ECG-peak-detection\\ECG Data'
# Go through a loop for every file
for filename in os.listdir(path):
    # Read data from CSV datafile

    print(os.getcwd())
    dataset = pd.read_csv('ECG Data/' + filename)

    # Get sampling rate from the differnce between the first two timestamps
    hz = 1 / (dataset["Time (s)"][1] - dataset["Time (s)"][0])
    # Get the amount of time for each sample based on the last timestamp
    seconds = dataset["Time (s)"][len(dataset["Time (s)"])-1]

    # Process the ECG data using biosppy
    out = ecg.ecg(signal=dataset["Signal (V)"],
                  sampling_rate=hz, show=True)

    filtered_ecg = out[1]

    # Using Hamilton's Method
    hamilton = ecg.hamilton_segmenter(filtered_ecg, sampling_rate=hz)

    # Using Engelse and Zeelenberg's Method
    engzee = ecg.engzee_segmenter(filtered_ecg, sampling_rate=hz)
    calculated_r_peaks = engzee[0]

    print(calculated_r_peaks)
    # Get the r_peaks from our processed data

    actual_r_peaks = []
    # Get the actual r-peaks from the original data
    for i in range(0, 23):  # The actual data is in the first 20 lines
        if dataset["Type"][i] == "N":
            actual_r_peaks.append(dataset["Sample #"][i])

    # Since we are dealing with messy data we need to define a tolerance for comparing our derived r peaks to the actual, 2 is a minimal tolerance
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

    else:
        print("Not the same number of r-peaks")
        if len(calculated_r_peaks) - len(actual_r_peaks) > 0:
            FP = len(calculated_r_peaks) - len(actual_r_peaks)
            FN = 0
            TP = len(actual_r_peaks)
        else:
            FN = len(actual_r_peaks) - len(calculated_r_peaks)
            FP = 0
            TP = len(calculated_r_peaks)

    sensitivity = TP / (TP+FN)
    calculated_heartrate_bpm = (len(calculated_r_peaks)/seconds)*60
    actual_heartrate_bpm = (len(actual_r_peaks)/seconds)*60

    file = open("results " + filename, 'w')
    file.writelines(["Calculated R Peaks:" + str(calculated_r_peaks) + '\n',
                     "Actual R Peaks:" + str(actual_r_peaks) + '\n',
                     "Sensitivity: " + str(sensitivity) + '\n',
                     "FP: " + str(FP) + '\n',
                     "TP: " + str(TP) + '\n',
                     "FN: " + str(FN) + '\n',
                     "Calculated Heartrate(bpm): " +
                     str(calculated_heartrate_bpm) + '\n',
                     "Actual Heartrate(bpm): " + str(actual_heartrate_bpm) + '\n'])
    file.close()
