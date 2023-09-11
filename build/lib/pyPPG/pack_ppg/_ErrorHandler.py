import numpy as np


def _check_shape_(signal, fs):
    assert len(signal) > fs*15, "Signal must be at list fifteen seconds"
    signal = np.array(signal)
    assert len(signal.shape) > 0, "Signal must not be empty"
    assert len(signal.shape) < 3, "Signal can be 1-dimensional array or a matrix with different signal in every column "

def _check_fragment_PRSA_(d):
    assert d > 0, "The parameter d should be strictly positive"


class WrongParameter(Exception):
    pass