from collections import namedtuple

#E: [MPa], nu: [-], rho: [Mg/mm^3]
YarnType = namedtuple('YarnType',['E','nu','rho_density'])
yarn_prop_paper = YarnType(1400.0,0.2,1.306e-09)
yarn_prop_pink = YarnType(51.5283,0.2,8.5532e-10)
#todo: recalc yarn_prop_pink E and rho w/ lower r

r_fake = 0.08
#r_pink = 0.8059 #this is what was measured vs below is what worked
#r_pink = 0.255
#r_pink = 0.245
r_pink = 0.24

StitchSize = namedtuple('StitchSize',['lam', 'w', 'gamma', 'CO', 'h', 'delta'])
stitch_size_paper = StitchSize(3.73,1.15,0.24,0.71,1.14,0.41)
stitch_size_5 = StitchSize(11.7123,3.0417,0.6045,2.5345,4.2130,2.0000)
