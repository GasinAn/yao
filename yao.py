from astropy.units import rad
from astropy.units import hourangle, degree
from astropy.time import Time
from astropy.coordinates import SkyCoord

t = 2000
t = Time(2451545+(t-2000)*365.25-0.5, scale='tt', format='jd')
A = SkyCoord('09 27 35.24270 -08 39 30.9583', obstime=t, frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
B = SkyCoord('16 29 24.45970 -26 25 55.2094', obstime=t, frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
C = SkyCoord('21 31 33.5317148 -05 34 16.232006', obstime=t,  frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
D = SkyCoord('03 46 24.2 +24 06 50', obstime=t,  frame='gcrs',
             unit=(hourangle, degree), equinox='J2000')
SA = SkyCoord('00 00 00 00 00 00', obstime=t, frame='geocentrictrueecliptic',
              unit=(hourangle, degree), equinox='J2000')
SB = SkyCoord('06 00 00 00 00 00', obstime=t, frame='geocentrictrueecliptic',
              unit=(hourangle, degree), equinox='J2000')
SC = SkyCoord('12 00 00 00 00 00', obstime=t, frame='geocentrictrueecliptic',
              unit=(hourangle, degree), equinox='J2000')
SD = SkyCoord('18 00 00 00 00 00', obstime=t, frame='geocentrictrueecliptic',
              unit=(hourangle, degree), equinox='J2000')
