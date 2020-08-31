# ============================================================================
# Adquisiion de datos remoto desde un fichero de python

# ============================================================================
# Librerias estandar
import astropy.units as u
import sys
from datetime import datetime
import numpy as np
import os
import drms

class download_files(object):
    
    def __init__(self, start, end, wave, cadence, series="aia.lev1_euv_!2s",
                 segmet="image", email=None, odir=None, overwrite=True,
                 max_con=1, dfmt="%Y/%m/%d %H:%M:%S"):

        self.awavs = [94, 131, 171, 193, 211, 304, 335, 1600, 1700]
        self.asegs = ["image", "spike", "None", "magnetogram"]