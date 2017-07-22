#!/usr/bi\nenv python 2.7
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

from sympy import *
import matplotlib.pyplot as plt


class Vector(object):
	vectorCount = 0

	def __init__(self, i=0, j=0, k=0):
		self.d = sqrt(i**2+j**2+k**2)
		self.i = i
		self.j = j
		self.k = k
		Vector.vectorCount += 1


def bargen(r, rp, rpp, t, tp, tpp, barnums):
	x = 0
	vectr = list()
	vectrp = list()
	vectrpp = list()
	while x < barnums:
		# r=r1(cos(t)I+sin(t)J
		vectr.append(Vector((r[x] * (cos(t[x]))), (r[x] * sin(t[x]))))
		# rp=rp(cos i+ sin j) +r(-tp sin I + tp cos J)
		vectrp.append(
			Vector(rp[x] * cos(t[x]) - r[x] * tp[x] * sin(t[x]), rp[x] * sin(t[x]) + r[x] * tp[x] * cos(t[x])))
		# rpp = rpp cos I sin J + r tpp -sin I + cos J - rtp^2 cos I + sinJ + 2rp tp -sinI + cos J
		vectrpp.append(Vector(
			rpp[x] * cos(t[x]) - r[x] * tpp[x] * sin(t[x]) - r[x] * ((tp[x]) ** 2) * cos(t[x]) - 2 * rp[x] * tp[
				x] * sin(t[x])
			, rpp[x] * sin(t[x]) + r[x] * tpp[x] * cos(t[x]) - r[x] * ((tp[x]) ** 2) * sin(t[x]) + 2 * rp[x] * tp[
				x] * cos(t[x])))
		x += 1
	return (vectr,vectrp,vectrpp)


def degtrad(t):
	rad = t*pi/180
	return float(rad)


def radtdeg(t):
	deg = t*180/pi
	return float(deg)


def q6():
	r=[]
	rp=[]
	rpp=[]
	t=[]
	tp=[]
	tpp=[]
	v,vp,vpp=bargen(r,rp,rpp,t,tp,tpp)
	eq1=Eq()
	eq2=Eq()
	eq3=Eq()
	eq4=Eq()
q6()
