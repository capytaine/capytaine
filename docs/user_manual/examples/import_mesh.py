#!/usr/bin/env python

from capytaine import FloatingBody

body = FloatingBody.from_file('Cylinder.dat')
body.name = "my cylinder"
body.show()
