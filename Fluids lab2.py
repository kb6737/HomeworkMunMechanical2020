#!/usr/bin/env python
#  Copyright 2017 Keegan Joseph Brophy
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

"""Usage:
	Fluidslab2.py execute (FUNCTION)
	Fluidslab2.py -h --help
Options:
	execute (docustomexport)
"""

import math


def fluidslab2_plate(diam, dens):
	d1 = diam
	# Fill up the list with the values as ordered in the lab.
	f = [1, 2, 3, 4, 5]
	s = [58.63, 40.49, 39.16, 33.49, 28.85]
	l = [15, 14.8, 18.2, 17, 16.2]
	n = -1
	w1 = list()
	fth = list()
	flow = list()
	while n < len(f) - 1:
		n += 1
		flow.append(l[n] / s[n])
		ad = (math.pi*math.pow(d1, 2)/4)
		print math.pow(d1, 2)
		w1.append((flow[n]*0.001 / ad))
		ro = dens
		fth.append(flow[n] * ro * w1[n])
	print 'flow rate '+str(flow)
	print 'w1 '+str(w1)
	print 'fth '+str(fth)


def fluidslab2_hemisphere(diam, dens):
	d1 = diam
	# Fill up the list with the values as ordered in the lab.
	f = [1,2,3,4,5]
	s = [123.87,49.59,46.80,40.31,39.58]
	l = [16,14,15,16,18]
	n = -1
	w1 = list()
	fth = list()
	flow = list()
	while n < len(f) - 1:
		n += 1
		flow.append(l[n] / s[n])
		ad = (math.pi * math.pow(d1, 2) / 4)
		w1.append((flow[n]*0.001 / ad))
		ro = dens
		fth.append(2 * flow[n] * ro * w1[n])
	print 'flow rate '+str(flow)
	print 'w1 '+str(w1)
	print 'fth '+str(fth)


def fluidslab2_slope(diam, dens):
	d1 = diam
	# Fill up the list with the values as ordered in the lab.
	f = []
	s = []
	l = []
	n = -1
	w1 = list()
	fth = list()
	flow = list()
	while n < len(f) - 1:
		n += 1
		flow.append(l[n] / s[n])
		ad = (math.pi * math.pow(d1, 2) / 4)
		w1.append((flow[n]*0.001 / ad))
		ro = dens
		fth.append(flow[n] * ro * w1[n]*math.pow(math.cos(45), 2))
	print 'flow rate '+str(flow)
	print 'w1 '+str(w1)
	print 'fth '+str(fth)


def fluidslab2_cone(diam, dens):
	d1 = diam
	# Fill up the list with the values as ordered in the lab.
	f = []
	s = []
	l = []
	n = -1
	w1 = list()
	fth = list()
	flow = list()
	while n < len(f) - 1:
		n += 1
		flow.append(l[n] / s[n])
		ad = (math.pi * math.pow(d1, 2) / 4)
		w1.append((flow[n]*0.001 / ad))
		ro = dens
		fth.append(flow[n] * ro * w1[n]*math.pow(math.cos(45), 2))
	print 'flow rate '+str(flow)
	print 'w1 '+str(w1)
	print 'fth '+str(fth)

fluidslab2_plate(0.01, 1)
fluidslab2_hemisphere(0.01,1)