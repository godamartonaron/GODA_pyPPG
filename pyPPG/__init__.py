from pack_ppg._ErrorHandler import _check_shape_, WrongParameter

class PPG:
    def __init__(self,s):
        """
        :param s: a struct of PPG signal:

            * s.v: a vector of PPG values
            * s.fs: the sampling frequency of the PPG in Hz
            * s.name: name of the record
            * s.v: 1-d array, a vector of PPG values
            * s.fs: the sampling frequency of the PPG in Hz
            * s.filt_sig: 1-d array, a vector of the filtered PPG values
            * s.filt_d1: 1-d array, a vector of the filtered PPG' values
            * s.filt_d2: 1-d array, a vector of the filtered PPG" values
            * s.filt_d3: 1-d array, a vector of the filtered PPG'" values
        :type s: DotMap

        """

        if s.fs <= 0:
            raise WrongParameter("Sampling frequency should be strictly positive")
        _check_shape_(s.v, s.fs)

        keys=s.keys()
        keys_list = list(keys)
        for i in keys_list:
            exec('self.'+i+' = s[i]')