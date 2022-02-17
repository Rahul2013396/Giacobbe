import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe_codes'))

import PCA as PCA


spectrum = PCA.datacollect(29,2)
spectrum  = PCA.meadiannorm(spectrum)
spectrum = PCA.normalisation(spectrum)
plt.imshow(spectrum)
plt.show()


