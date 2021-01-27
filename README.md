# 📈ECG-peak-detection📉

## This code takes sample ECG data and attempts to detect R Peaks.

### It does this with two different algorithims:

#### [Hamilton](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2532677/)

#### [Engelse and Zeelenberg](http://www.lx.it.pt/~afred/papers/Review%20and%20Comparison%20of%20Real%20Time%20Electrocardiogram%20Segmentation%20Algorithms%20for%20Biometric%20Applications.pdf)

### It also creates a result file with the relevant data, calculated bpm, and sensitivity compared to actual.

### In addition it creates accompanying visualizations:

![Sample Visualization](https://github.com/Shaun-Regenbaum/ECG-peak-detection/blob/main/results/Engelse%20and%20Zeelenberg/graphics%201a.jpg?raw=true)

## Further info:

In my repo you will find the following structure:

- ECG-peak-detection
  - Results
    - Engelse and Zeelenberg
    - Hamilton
  - ECG Data
  - peak-detection.py
  - engelse-and-zeelenberg.py (simply for reference, not used)
  - README.md (because of Github)

I used the python library biosppy for filtering the ECG data.
I use NumPy, and Pandas for data manipulation, I use OS for saving to files, and use Matplotlib for visualizations.
You will find that everything I use laid out and my program is well documented in the peak-detection.py file.
In that file, there are two initial functions which are just for checking if the r-peaks are within range of the given r-peaks (for calculating FP and FN), and a function that creates the graphics for every case.
You can probably skip those two and go straight to the main bunch of code afterwards.

To lay it out clearly here:

I used two different methods for identifying r-peaks.

The first was Hamilton’s algorithm for segmentation which can be found here (or at least a rough variation of it):

(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2532677/)

This did not perform as well as the benchmark, thus I then continued to look for others.
I decided upon using Engelse and Zeelenberg whom proved that their algorithm outperformed Hamilton’s (among many others). Beyond this it is much easier to both implement and understand:

http://www.lx.it.pt/~afred/papers/Review%20and%20Comparison%20of%20Real%20Time%20Electrocardiogram%20Segmentation%20Algorithms%20for%20Biometric%20Applications.pdf

This time, I was 100% accurate.

Meaning: **There were no False Positive nor False Negatives, everything was a True Positive. This means that my sensitivity was 1 for every case.**

You can look at all my results in the Results folder.
