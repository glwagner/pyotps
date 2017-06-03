import sys; sys.path.append('../pyotps')
import pyotps

otpspath_nc = '/data5/glwagner/Numerics/nestgcm/pyotps/OTPSnc'
otpspath_v1 = '/data5/glwagner/Numerics/nestgcm/pyotps/OTPS2'

driver = pyotps.tidaldriver(otpspath=otpspath_nc, modeltype='netcdf')
driver.testdrive()

driver = pyotps.tidaldriver(otpspath=otpspath_v1, modeltype='v1')
driver.testdrive()
