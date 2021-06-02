import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mundey

import pytest

ddir = '../notebooks/'

tpffile = ddir+'tess2018349182459-s0006-0000000322899250-0126-s_tp.fits'

def test_sirius():
	tpf = mundey.mundey.mundey_tpf(tpffile)

	tpf.calibrate() # do the entire thing

	assert np.nansum(tpf.flux.value) == 4128732700.0