#import sys; sys.path.append('../pyotps')
import pyotps.pyotps

otpspath = '/data5/glwagner/Numerics/patches/pyotps/OTPS2'

driver = pyotps.TidalDriver(otpspath=otpspath, otpstype='v1')
driver.testdrive()
