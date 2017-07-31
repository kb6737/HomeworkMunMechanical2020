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


def twoangsolve1(r, t):
	A=-2*r[0]*r[3]*cos(t[0])-2*r[3]*r[2]
	B=-2*r[3]*r[0]*sin(t[0])
	C=r[3]**2+r[2]**2+r[0]**2-r[1]**2+2*r[2]*r[0]*cos(t[0])
	x = symbols('x')
	x,y=solve_linear(x,(B+(B**2-C**2+A**2)**.5)/(C-A))
	z,w=solve_linear(x,(B-(B**2-C**2+A**2)**.5)/(C-A))
	t32 = (-2 * atan(w))
	return float(t32)

def bargen(r, rp, rpp, t, tp, tpp, barnums=4):
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

#def q6(bar2,bar3,bar4,bar5,bar6,bar7,bar8):
def q6(bar2,bar3,bar4,bar5,bar6):
	master1 = list()
	master2 = list()
	master3 = list()
	master4 = list()
	master5 = list()
	master6 = list()
	master7 = list()
	master8 = list()
	master9 = list()
	master10 = list()
	master11 = list()
	master12 = list()
	master13 = list()
	master14 = list()
	master15 = list()
	x=(14**2+4**2-2*14*4*cos(degtrad(45)))**.5
	var2=x
	w = 20000 * 9.81
	n=0.0
	while x <= bar2+bar3 :
		t0 = acos((x ** 2 + bar2 ** 2 - bar3 ** 2) / (2 * x * bar2))
		t2 = acos((bar3 ** 2 + bar2 ** 2 - x ** 2) / (2 * bar3 * bar2))
		tp0,tp2,tpp0,tpp2=var('tp0 tp2 tpp0 tpp2')
		r=[x,bar2,bar3]
		rp=[0,1,0]
		rpp=[0,0,0]
		t=[(pi-t0),pi,t2]
		tp=[tp0,0,tp2]
		tpp=[tpp0,0,tpp2]
		vr,vrp,vrpp=bargen(r,rp,rpp,t,tp,tpp,3)
		eq1=Eq(vrp[0].i,vrp[1].i+vrp[2].i)
		eq2=Eq(vrp[0].j,vrp[1].j+vrp[2].j)
		a,b=solve_linear(eq1)
		c,d=solve_linear(eq2)
		e,f=solve_linear(b,d)
		tp=[tp0,0,f]
		vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp, 3)
		eq1 = Eq(vrp[0].i, vrp[1].i + vrp[2].i)
		a,b=solve_linear(eq1)
		tp=[b,0,f]
		vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp, 3)
		eq1 = Eq(vrpp[0].i, vrpp[1].i + vrpp[2].i)
		eq2 = Eq(vrpp[0].j, vrpp[1].j + vrpp[2].j)
		a, b = solve_linear(eq1)
		c, d = solve_linear(eq2)
		e, f = solve_linear(b, d)
		tpp=[f,0,tpp2]
		vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp, 3)
		eq1 = Eq(vrpp[0].i, vrpp[1].i + vrpp[2].i)
		a, b = solve_linear(eq1)
		tpp=[f,0,b]
		vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp, 3)
		r1=[((2*(bar3**2))**.5),bar4,bar5,bar6]
		r1p=[0,0,0,0]
		r1pp=[0,0,0,0]
		t1=[(t[2]-degtrad(45)),var('t11'),pi,var('t13')]
		a=twoangsolve1(r1,t1)
		t1=[t1[0],t1[1],pi,a]
		t1p=[tp[2],var('tp11'),0,var('tp13')]
		t1pp=[tpp[2],var('tpp11'),0,var('tpp13')]
		v1r,v1rp,v1rpp=bargen(r1,r1p,r1pp,t1,t1p,t1pp,4)
		eq1 = Eq(v1r[0].i+v1r[1].i,v1r[2].i+v1r[3].i)
		a, b = solve(eq1)
		t1=[t1[0],a,pi,t1[3]]
		v1r,v1rp,v1rpp=bargen(r1,r1p,r1pp,t1,t1p,t1pp,4)
		eq1 = Eq(v1rp[0].i+v1rp[1].i,v1rp[2].i+v1rp[3].i)
		eq2 = Eq(v1rp[0].j+v1rp[1].j,v1rp[2].j+v1rp[3].j)
		a, b = solve_linear(eq1)
		c, d = solve_linear(eq2)
		e, f = solve_linear(b, d)
		t1p=[tp[2],var('tp11'),0,f]
		v1r, v1rp, v1rpp = bargen(r1, r1p, r1pp, t1, t1p, t1pp, 4)
		eq1 = Eq(v1rp[0].i + v1rp[1].i, v1rp[2].i + v1rp[3].i)
		a, b = solve_linear(eq1)
		t1p=[tp[2],b,0,f]
		v1r, v1rp, v1rpp = bargen(r1, r1p, r1pp, t1, t1p, t1pp, 4)
		eq1 = Eq(v1rpp[0].i + v1rpp[1].i, v1rpp[2].i + v1rpp[3].i)
		eq2 = Eq(v1rpp[0].j + v1rpp[1].j, v1rpp[2].j + v1rpp[3].j)
		a, b = solve_linear(eq1)
		c, d = solve_linear(eq2)
		e, f = solve_linear(b, d)
		t1pp=[tpp[2],f,0,var('tpp13')]
		v1r, v1rp, v1rpp = bargen(r1, r1p, r1pp, t1, t1p, t1pp, 4)
		eq1 = Eq(v1rpp[0].i + v1rpp[1].i, v1rpp[2].i + v1rpp[3].i)
		a, b = solve_linear(eq1)
		t1pp=[tpp[2],f,0,b]
		v1r, v1rp, v1rpp = bargen(r1, r1p, r1pp, t1, t1p, t1pp, 4)
		master1.append(vr[0].i)
		master2.append(vr[0].j)
		master3.append(v1r[0].i)
		master4.append(v1r[0].j)
		master5.append(v1r[0].i+v1r[1].i)
		master6.append(v1r[0].j+v1r[0].j)
		fad,fde,ffe,fcb=var('fad fde ffe fcb')
		eq1=Eq(ffe*cos(t1[3]),fde*cos(pi-t1[1]))
		eq2=Eq(ffe*sin(t1[3])+fde*sin(pi-t1[1]),w)
		a, b = solve_linear(eq1)
		c, d = solve_linear(eq2)
		e, f = solve_linear(b, d)
		#check e,f for fde or ffe it changes based on the computer type and i check this myself every time i run the code on a new cpu
		fde=f
		eq1=Eq(ffe*cos(t1[3]),fde*cos(pi-t1[1]))
		a, b = solve_linear(eq1)
		ffe=b
		master7.append(fde)
		master8.append(ffe)
		fdei = fde * cos(pi - t1[1])
		fdej = fde * sin(pi - t1[1])
		eq1=Eq(r[2]*cos(t1[0]+degtrad(45))*sin(pi-t[0])*fcb+r[2]*sin(t1[0]+degtrad(45))*cos(pi-t[0])*fcb,r1[0]*cos(t1[0])*fdej+r1[0]*sin(t1[0])*fdei)
		a, b = solve_linear(eq1)
		fcb=b
		if fcb>0:
			master9.append(fcb)
			master10.append(t1[3])
		master11.append(r[0])
		master12.append(t1[3])
		master13.append(rp[1])
		n+=.1
		master14.append(n)
		master15.append(fcb*(x-var2)/n)
		x+=.1
	plt.figure(1)
	plt.plot(master1, master2)
	plt.title('Position analysis of B')
	plt.xlabel('point B i')
	plt.ylabel('Point B j')
	plt.savefig('1')

	plt.figure(2)
	plt.plot(master3, master4)
	plt.title('Position analysis of D')
	plt.xlabel('point D i')
	plt.ylabel('Point D j')
	plt.savefig('2')

	plt.figure(3)
	plt.plot(master5, master6)
	plt.title('Position analysis of E')
	plt.xlabel('point E i')
	plt.ylabel('Point E j')
	plt.savefig('3')

	plt.figure(4)
	plt.plot(master7,master8)
	plt.title('Force analysis of F')
	plt.xlabel('Bar Fde')
	plt.ylabel('Bar Ffe')
	plt.savefig('4')

	plt.figure(5)
	plt.plot(master10,master9)
	plt.title('Force analysis of piston')
	plt.xlabel('Theta 6')
	plt.ylabel('Bar Fcb')
	plt.savefig('5')

	plt.figure(6)
	plt.plot(master12, master11)
	plt.title('Position')
	plt.ylabel('Lenth of piston')
	plt.xlabel('Theta 6')
	plt.savefig('6')

	plt.figure(7)
	plt.plot(master13,master12)
	plt.title('Force analysis of piston')
	plt.xlabel('Velocity of piston')
	plt.ylabel('Theta 6')
	plt.savefig('7')

	plt.figure(8)
	plt.plot(master14,master15)
	plt.title('Power')
	plt.xlabel('Time')
	plt.ylabel('Force*Distance/time')
	plt.savefig('8')
	plt.show()

q6(14,4,4.5,3.2,7.75)
