def engzee_segmenter(signal=None, sampling_rate=1000., threshold=0.48):
    """ECG R-peak segmentation algorithm.

    Follows the approach by Engelse and Zeelenberg [EnZe79]_ with the
    modifications by Lourenco *et al.* [LSLL12]_.

    Parameters
    ----------
    signal : array
        Input filtered ECG signal.
    sampling_rate : int, float, optional
        Sampling frequency (Hz).
    threshold : float, optional
        Detection threshold.

    Returns
    -------
    rpeaks : array
        R-peak location indices.

    References
    ----------
    .. [EnZe79] W. Engelse and C. Zeelenberg, "A single scan algorithm for
       QRS detection and feature extraction", IEEE Comp. in Cardiology,
       vol. 6, pp. 37-42, 1979
    .. [LSLL12] A. Lourenco, H. Silva, P. Leite, R. Lourenco and A. Fred,
       "Real Time Electrocardiogram Segmentation for Finger Based ECG
       Biometrics", BIOSIGNALS 2012, pp. 49-54, 2012

    """

    # check inputs
    if signal is None:
        raise TypeError("Please specify an input signal.")

    # algorithm parameters
    changeM = int(0.75 * sampling_rate)
    Miterate = int(1.75 * sampling_rate)
    v250ms = int(0.25 * sampling_rate)
    v1200ms = int(1.2 * sampling_rate)
    v1500ms = int(1.5 * sampling_rate)
    v180ms = int(0.18 * sampling_rate)
    p10ms = int(np.ceil(0.01 * sampling_rate))
    p20ms = int(np.ceil(0.02 * sampling_rate))
    err_kill = int(0.01 * sampling_rate)
    inc = 1
    mmth = threshold
    mmp = 0.2

    # Differentiator (1)
    y1 = [signal[i] - signal[i - 4] for i in range(4, len(signal))]

    # Low pass filter (2)
    c = [1, 4, 6, 4, 1, -1, -4, -6, -4, -1]
    y2 = np.array([np.dot(c, y1[n - 9:n + 1]) for n in range(9, len(y1))])
    y2_len = len(y2)

    # vars
    MM = mmth * max(y2[:Miterate]) * np.ones(3)
    MMidx = 0
    Th = np.mean(MM)
    NN = mmp * min(y2[:Miterate]) * np.ones(2)
    NNidx = 0
    ThNew = np.mean(NN)
    update = False
    nthfpluss = []
    rpeaks = []

    # Find nthf+ point
    while True:
        # If a previous intersection was found, continue the analysis from there
        if update:
            if inc * changeM + Miterate < y2_len:
                a = (inc - 1) * changeM
                b = inc * changeM + Miterate
                Mnew = mmth * max(y2[a:b])
                Nnew = mmp * min(y2[a:b])
            elif y2_len - (inc - 1) * changeM > v1500ms:
                a = (inc - 1) * changeM
                Mnew = mmth * max(y2[a:])
                Nnew = mmp * min(y2[a:])
            if len(y2) - inc * changeM > Miterate:
                MM[MMidx] = Mnew if Mnew <= 1.5 * \
                    MM[MMidx - 1] else 1.1 * MM[MMidx - 1]
                NN[NNidx] = Nnew if abs(
                    Nnew) <= 1.5 * abs(NN[NNidx - 1]) else 1.1 * NN[NNidx - 1]
            MMidx = np.mod(MMidx + 1, len(MM))
            NNidx = np.mod(NNidx + 1, len(NN))
            Th = np.mean(MM)
            ThNew = np.mean(NN)
            inc += 1
            update = False
        if nthfpluss:
            lastp = nthfpluss[-1] + 1
            if lastp < (inc - 1) * changeM:
                lastp = (inc - 1) * changeM
            y22 = y2[lastp:inc * changeM + err_kill]
            # find intersection with Th
            try:
                nthfplus = np.intersect1d(np.nonzero(
                    y22 > Th)[0], np.nonzero(y22 < Th)[0] - 1)[0]
            except IndexError:
                if inc * changeM > len(y2):
                    break
                else:
                    update = True
                    continue
            # adjust index
            nthfplus += int(lastp)
            # if a previous R peak was found:
            if rpeaks:
                # check if intersection is within the 200-1200 ms interval. Modification: 300 ms -> 200 bpm
                if nthfplus - rpeaks[-1] > v250ms and nthfplus - rpeaks[-1] < v1200ms:
                    pass
                # if new intersection is within the <200ms interval, skip it. Modification: 300 ms -> 200 bpm
                elif nthfplus - rpeaks[-1] < v250ms:
                    nthfpluss += [nthfplus]
                    continue
        # no previous intersection, find the first one
        else:
            try:
                aux = np.nonzero(
                    y2[(inc - 1) * changeM:inc * changeM + err_kill] > Th)[0]
                bux = np.nonzero(
                    y2[(inc - 1) * changeM:inc * changeM + err_kill] < Th)[0] - 1
                nthfplus = int((inc - 1) * changeM) + \
                    np.intersect1d(aux, bux)[0]
            except IndexError:
                if inc * changeM > len(y2):
                    break
                else:
                    update = True
                    continue
        nthfpluss += [nthfplus]
        # Define 160ms search region
        windowW = np.arange(nthfplus, nthfplus + v180ms)
        # Check if the condition y2[n] < Th holds for a specified
        # number of consecutive points (experimentally we found this number to be at least 10 points)"
        i, f = windowW[0], windowW[-1] if windowW[-1] < len(y2) else -1
        hold_points = np.diff(np.nonzero(y2[i:f] < ThNew)[0])
        cont = 0
        for hp in hold_points:
            if hp == 1:
                cont += 1
                if cont == p10ms - 1:  # -1 is because diff eats a sample
                    max_shift = p20ms  # looks for X's max a bit to the right
                    if nthfpluss[-1] > max_shift:
                        rpeaks += [np.argmax(signal[i -
                                                    max_shift:f]) + i - max_shift]
                    else:
                        rpeaks += [np.argmax(signal[i:f]) + i]
                    break
            else:
                cont = 0

    rpeaks = sorted(list(set(rpeaks)))
    rpeaks = np.array(rpeaks, dtype='int')

    return utils.ReturnTuple((rpeaks,), ('rpeaks',))
