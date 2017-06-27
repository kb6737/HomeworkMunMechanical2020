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

import math


def run1(vin, lvl1, lvl2, d, u, l, k):
	lvlloss = list()
	# lAmbda = 0.3164/math.pow(re, 1./4)
	# laMbda = math.pow(2*math.log((2.51/(re*math.sqrt(lAmbda)))+ 0.27/(d/k)),-2)
	v = list()
	vdot = list()
	re = list()
	sigma = list()
	thl = list()
	lAmbda = list()
	lhvms = 2.777778 * math.pow(10, -7)
	n=-1
	g = 9.81
	for x in lvl1:
		n+=1
		lvlloss.append((lvl1[n]-lvl2[n])*0.001)
		vdot.append(vin[n]*lhvms)
		v.append((4*vdot[n])/(math.pi * d * d))
		re.append(v[n] * d / u)
		if (re[n] > 2320) & (re[n] < math.pow(10, 5)):
			lAmbda.append(0.3164/math.pow(re[n], 1./4))
			sigma.append((2*lvlloss[n]*g)/(v[n]*v[n])-(lAmbda[n]*l*d))
			thl.append((lAmbda[n]*l/d)*(v[n]*v[n]/(2*g)))
		elif re[n] > math.pow(10, 5):
			lAmbda.append(math.pow(2 * math.log((2.51 / (re[n] * math.sqrt(lAmbda[n]))) + 0.27 / (d / k)), -2))
			sigma.append((2*lvlloss[n]*g)/(v[n]*v[n])-(lAmbda[n]*l*d))
			thl.append((lAmbda[n]*l/d)*(v[n]*v[n]/(2*g)))
	print 'Loss level m '+str(lvlloss)
	print 'Flow rate m^3/s '+str(vdot)
	print 'Flow velocity m/s '+str(v)
	print "Renold's number "+str(re)
	print 'Coefficient of pipe friction '+str(lAmbda)
	print 'Length of pipe m '+str(l)
	print 'Pipe roughness m '+str(k)
	print 'Kinematic viscosity m^2/s '+str(u)
	print 'Inside diameter of pipe m '+str(d)
	print 'Resistance coefficient '+str(sigma)
	print 'theoretical loss '+str(thl)

cPID = 0.016
cPOD = 0.018
pVCPID = 0.017
pCVPOD = 0.020
# copper
t1vin = [200, 400, 600, 800, 1000, 1200, 1400]
t1lvl1 = [591, 590, 599, 610, 619, 627, 635]
t1lvl2 = [576, 552, 529, 497, 458, 415, 363]
# PVC
t2lvl1 = [573, 576, 589, 598, 618, 633, 656]
t2lvl2 = [558, 534, 503, 466, 421, 363, 399]
# knee
t3lvl1 = [511, 509, 512, 507, 499, 493, 483]
t3lvl2 = [505, 490, 474, 446, 412, 377, 324]
# elbow
t4lvl1 = [540, 540, 540, 560, 580, 590, 620]
t4lvl2 = [510, 490, 460, 430, 380, 330, 260]
u = 1.297
l = 1
k = (0.001)*1/1000
run1(t1vin, t1lvl1, t1lvl2, cPID, u, l, k)