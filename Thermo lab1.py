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
import matplotlib


def mean(numbers):
	return math.fsum(numbers)/len(numbers)


def interpolation(known1, known2, unknown1, known3, known4):
	unknown2 = known3+((unknown1-known1)*(known4-known3)/(known2-known1))
	return unknown2


def thermolab1():
	# all values in [] are for the user to use.
	dt = [5, 5, 10]
	cdt = math.fsum(dt)
	dt = mean(dt)
	p1 = [4, 4, 4]
	cp1 = math.fsum(p1)
	p1 = mean(p1)
	pamb = [1, 1, 1]
	cpamb = math.fsum(pamb)
	pamb = mean(pamb)
	t1 = [17.9, 17.8, 17.7]
	ct1 = math.fsum(t1)
	t1 = mean(t1)
	t2 = [25.7, 20.7, 21.4]
	ct2 = math.fsum(t2)
	t2 = mean(t2)
	t3 = [145.2, 145.1, 145.2]
	ct3 = math.fsum(t3)
	t3 = mean(t3)
	t4 = [116.4, 116.1, 116]
	ct4 = math.fsum(t4)
	t4 = mean(t4)
	t5 = [101.3, 101.2, 101.1]
	ct5 = math.fsum(t5)
	t5 = mean(t5)
	t6 = [32, 32.8, 33.5]
	ct6 = math.fsum(t6)
	t6 = mean(t6)
	t7 = [11.1, 11, 10.9]
	ct7 = math.fsum(t7)
	t7 = mean(t7)
	t8 = [17.8, 17.5, 17.7]
	ct8 = math.fsum(t8)
	t8 = mean(t8)
	vdotcw = [55, 55, 55]
	cvdotcw = math.fsum(vdotcw)
	vdotcw = mean(vdotcw)
	vdotgas = [55, 55, 55]
	cvdotgas = math.fsum(vdotgas)
	vdotgas = mean(vdotgas)
	vc = [72.25, 75, 135]
	cvc = math.fsum(vc)
	vc = mean(vc)
	u = [6.2, 6.1, 6]
	cu = math.fsum(u)
	u = mean(u)
	i = [158, 145, 150]
	ci = math.fsum(i)
	i = mean(i)
	ron = 201
	roc = 1000
	cpqout = 4.1855
	# maybe 2.01
	hu = 46354
	to = 273.15
	print 'dt '+str(dt)
	print 'p1 '+str(p1)
	print 'pamb '+str(pamb)
	print 't1 '+str(t1)
	print 't2 '+str(t2)
	print 't3 '+str(t3)
	print 't4 '+str(t4)
	print 't5 '+str(t5)
	print 't6 '+str(t6)
	print 't7 '+str(t7)
	print 't8 '+str(t8)
	print 'vdotcw '+str(vdotcw)
	print 'vdotgas '+str(vdotgas)
	print 'vc '+str(vc)
	print 'u '+str(u)
	print 'i '+str(i)
	print 'cdt '+str(cdt)
	print 'cp1 '+str(cp1)
	print 'cpamb '+str(cpamb)
	print 'ct1 '+str(ct1)
	print 'ct2 '+str(ct2)
	print 'ct3 '+str(ct3)
	print 'ct4 '+str(ct4)
	print 'ct5 '+str(ct5)
	print 'ct6 '+str(ct6)
	print 'ct7 '+str(ct7)
	print 'ct8 '+str(ct8)
	print 'cvdotcw '+str(cvdotcw)
	print 'cvdotgas '+str(cvdotgas)
	print 'cvc '+str(cvc)
	print 'cu '+str(cu)
	print 'ci '+str(ci)
	wt = u*i
	dhqout = cpqout*(t8-t7)
	# p0 listed but we do not have a value for p0 so assume p1
	p0 = p1
	rogass = ron*(to/(to+t1))*(pamb/p0)
	mdotgass = (vdotgas*rogass)/(3600*1000)
	pgas = mdotgass*hu
	mdotst = roc*vc/dt
	qout = mdotst*dhqout

	# pst = mdotst * (h3pp-h3p)
	# nb = pst/pgas
	# pcw = vdotcw*rocw*cp*(t8-t7)
	# pcp = mdotst*cp*(100-t6)
	# pt = pst-(pcw-pcp)
	# nth = pt/pst
	# pel = u*i
	# nel = pel/pt
	# n = pel/pgas

thermolab1()
