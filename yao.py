from math import pi
from astropy.units import hourangle, degree, rad
from astropy.time import Time
from astropy.coordinates import SkyCoord

t = 2000
t = Time(2451545+(t-2000)*365.25, scale='tt', format='jd')
A = SkyCoord('09 27 35.24270 -08 39 30.9583', obstime=t, frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
B = SkyCoord('16 29 24.45970 -26 25 55.2094', obstime=t, frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
C = SkyCoord('21 31 33.5317148 -05 34 16.232006', obstime=t,  frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
D = SkyCoord('03 46 24.2 +24 06 50', obstime=t,  frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
print(A.ra.to(rad),B.ra.to(rad),
      C.ra.to(rad),D.ra.to(rad))
print(A.dec.to(rad),B.dec.to(rad),
      C.dec.to(rad),D.dec.to(rad))
print(A.cirs.ra.to(rad),B.cirs.ra.to(rad),
      C.cirs.ra.to(rad),D.cirs.ra.to(rad))
print(A.cirs.dec.to(rad),B.cirs.dec.to(rad),
      C.cirs.dec.to(rad),D.cirs.dec.to(rad))
