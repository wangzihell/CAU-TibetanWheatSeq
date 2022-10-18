#!/usr/bin/env python

import dadi
from sys import argv

data = dadi.Spectrum.from_file(argv[1])
print(dadi.Spectrum.Watterson_theta(data))
