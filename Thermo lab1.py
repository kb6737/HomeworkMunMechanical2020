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
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


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
	vc = mean(vc)/100000000
	u = [6.2, 6.1, 6]
	cu = math.fsum(u)
	u = mean(u)
	i = [158, 145, 150]
	ci = math.fsum(i)
	i = mean(i)*0.001
	print
	print 'dt min '+str(dt)
	print 'p1 bar '+str(p1)
	print 'pamb bar '+str(pamb)
	print 't1 c '+str(t1)
	print 't2 c '+str(t2)
	print 't3 c '+str(t3)
	print 't4 c '+str(t4)
	print 't5 c '+str(t5)
	print 't6 c '+str(t6)
	print 't7 c '+str(t7)
	print 't8 c '+str(t8)
	print 'vdotcw l/h '+str(vdotcw)
	print 'vdotgas l/h '+str(vdotgas)
	print 'vc cm^3 '+str(vc)
	print 'u V '+str(u)
	print 'i mA '+str(i)
	print
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
	print
	ron = 201
	roc = 1000
	rocw = 1000
	cppcw = 4.1855
	cppcp = interpolation(30, 100, 32.766666, 4.178, 4.219)
	# maybe 201
	hu = 46354
	to = 273.15
	print 'ron kg/m^3 '+str(ron)
	print 'cppcw kj/kg*k '+str(cppcw)
	print 'cppcp kj/kg*k '+str(cppcp)
	print 'to k '+str(to)
	print 'hu kj/kg '+str(hu)
	wt = u*i
	print 'wt mW '+str(wt)
	dhqout = cppcw*(t8-t7)
	print 'dhqout kj/kg '+str(dhqout)
	# p0 listed but we do not have a value for p0 so assume p1
	p0 = p1
	rogass = ron*(to/(to+t1))*(pamb/p0)
	print 'rogass kg/m^3 '+str(rogass)
	mdotgass = (vdotgas*rogass)/(3600*1000)
	print 'mdotgass kg/s '+str(mdotgass)
	pgas = mdotgass*hu
	print 'Qin kW '+str(pgas)
	mdotst = (roc*vc/dt)*0.1667
	print 'mdotst kg/s '+str(mdotst)
	qout = mdotst*dhqout
	# slightly super heated
	h2 = 604.74
	s2 = 1.7766
	h3 = 2732.4
	s3 = 6.8959
	h4 = 2698.8
	s4 = s5 = interpolation(110, 120, 116, 7.2387, 7.1296)
	sg5 = interpolation(100, 110, 101, 7.3549, 7.2387)
	sf5 = interpolation(100, 110, 101, 1.3069, 1.4185)
	hg5 = interpolation(100, 110, 101, 2676.1, 2691.5)
	hf5 = interpolation(100, 110, 101, 419.04, 461.30)
	x = (s5-sf5)/(sg5-sf5)
	h5 = hf5 + x*(hg5-hf5)
	h6 = interpolation(32, 33, 32.766666, 129.97, 134.15)
	s6 = interpolation(32, 33, 32.766666, 0.4507, 0.4644)
	pst = mdotst*(h4-h5)
	print 'pst kW '+str(pst)
	nb = pst/pgas
	print 'nb efficiency '+str(nb)
	pcw = ((vdotcw)*rocw*(cppcw*(t8-t7)))*0.000000277386
	print 'pcw kW '+str(pcw)
	pcp = mdotst*cppcp*(100-t6)
	print 'pcp kW '+str(pcp)
	pt = pst-(pcw-pcp)
	print 'pt kW '+str(pt)
	nth = pt/pst
	pel = u*i
	nel = pel/pt
	n = pel/pgas
	print
	print 'nth '+str(nth)
	print 'nel '+str(nel)
	print 'n '+str(n)
	# fig, ax = plt.subplots()
	# Path = mpath.Path
	# path_data = [
	# 	(Path.MOVETO, (s2, t2)),
	# 	(Path.LINETO, (s3, t3)),
	# 	(Path.LINETO, (s4, t4)),
	# 	(Path.LINETO, (s5, t5)),
	# 	(Path.LINETO, (s6, t6)),
	# ]
	# codes, verts = zip(*path_data)
	# path = mpath.Path(verts, codes)
	# # plot control points and connecting lines
	# x, y = zip(*path.vertices)
	# line, = ax.plot(x, y, 'go-')
	# ax.grid()
	# ax.axis()
	# plt.show()



thermolab1()
