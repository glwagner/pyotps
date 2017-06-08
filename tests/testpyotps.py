import sys; sys.path.append('../pyotps')
import pyotps

otpspath = '/data5/glwagner/Numerics/patches/pyotps/OTPS2'

driver = pyotps.tidaldriver(otpspath=otpspath, modeltype='v1')
driver.testdrive()
