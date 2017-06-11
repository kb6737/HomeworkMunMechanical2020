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


def mean(numbers):
	return math.fsum(numbers)/len(numbers)


def thermolab1():
	time = [5, 5, 10]
	dt = mean(time)
	print dt
	p1 = [4, 4, 4]
	pamb = [1, 1, 1]
	t1 = [17.9, 17.8, 17.7]
	t2 = [25.7, 20.7, 21.4]
	t3 = [145.2, 145.1, 145.2]
	t4 = [116.4, 116.1, 116]
	t5 = [101.3, 101.2, 101.1]
	t6 = [32, 32.8, 33.5]
	t7 = [11.1, 11, 10.9]
	t8 = [17.8, 17.5, 17.7]
	vdotcw = [55, 55, 55]
	vdotgas = [55, 55, 55]
	vc = [72.25, 75, 135]
	u = [6.2, 6.1, 6]
	i = [158, 145, 150]
	ron = 201
	# maybe 2.01
	hu = 46354
	to = 273.15
	rogass = ron*(to/(to+t1))*(pamb/po)
	mdotgass = (vdotgas*rogass)/(3600*1000)
	pgas = mdotgass*hu
	mdotst = roc*vc/time
	pst = mdotst * (h3pp-h3p)
	nb = pst/pgas
	pcw = vdotcw*rocw*cp*(t8-t7)
	pcp = mdotst*cp*(100-t6)
	pt = pst-(pcw-pcp)
	nth = pt/pst
	pel = u*i
	nel = pel/pt
	n = pel/pgas