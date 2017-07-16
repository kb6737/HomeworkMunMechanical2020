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


def accsolve1(r, t, rp, tp, rpp, tpp):
	t = [t[2], t[0], t[1], t[3]]
	r = [r[2], r[0], r[1], r[3]]
	rp = [rp[2], rp[0], rp[1], rp[3]] # include this in nxt update
	tp = [tp[2], tp[0], tp[1], tp[3]]
	rpp = [rpp[2], rpp[0], rpp[1], rpp[3]]
	tpp = [tpp[2], tpp[0], tpp[1], tpp[3]]
	r1tpps=r[1]*tpp[1]*sin(t[1])
	r1tppc=r[1]*tpp[1]*cos(t[1])
	r1tpc=r[1]*tp[1]*tp[1]*cos(t[1])
	r1tps=r[1]*tp[1]*tp[1]*sin(t[1])
	r2tpc=r[2]*tp[2]*tp[2]*cos(t[2])
	r2tps=r[2]*tp[2]*tp[2]*sin(t[2])
	r3tpc=r[3]*tp[3]*tp[3]*cos(t[3])
	r3tps=r[3]*tp[3]*tp[3]*sin(t[3])
	A = Matrix([[-r[2]*sin(t[2]),r[3]*sin(t[3])],[-r[2]*cos(t[2]),r[3]*cos(t[3])]])
	B = Matrix([[r1tpps+r1tpc+r2tpc-r3tpc],[r1tppc-r1tps-r2tps+r3tps]])
	C,D = A.inv()*B
	return C,D

def veosolve1(r, t, rp, tp):
	t = [t[2], t[0], t[1], t[3]]
	r = [r[2], r[0], r[1], r[3]]
	rp = [rp[2], rp[0], rp[1], rp[3]] # include this in nxt update
	tp = [tp[2], tp[0], tp[1], tp[3]]
	A = Matrix([[-r[2]*sin(t[2]),r[3]*sin(t[3])],[-r[2]*cos(t[2]),r[3]*cos(t[3])]])
	B = Matrix([[r[1]*tp[1]*sin(t[1])],[r[1]*tp[1]*cos(t[1])]])
	C,D = A.inv()*B
	return C,D


def twoangsolve2(r, t):
	t = [t[2], t[0], t[1], t[3]]
	r = [r[2], r[0], r[1], r[3]]
	A = 2 * r[0] * r[3] * cos(t[0]) - 2 * r[1] * r[3] * cos(t[1])
	B = 2 * r[0] * r[3] * sin(t[0]) - 2 * r[1] * r[3] * sin(t[1])
	C = r[0] ** 2 + r[1] ** 2 + r[3] ** 2 - r[2] ** 2 - 2 * r[0] * r[1] * (cos(t[0]) * cos(t[1]) + sin(t[0]) * sin(t[1]))
	x = symbols('x')
	z, y = solve((C - A)*x*x + 2*B*x + A + C, x)
	t32 = 2 * atan(z)
	t2 = atan2((r[0]*sin(t[0])+r[3]*sin(t32)-r[1]*sin(t[1]))/r[2], (r[0]*cos(t[0])+r[3]*cos(t32)-r[1]*cos(t[1]))/r[2])
	return float(t32), float(t2)


def twoangsolve1(r, t):
	t = [t[2], t[0], t[1], t[3]]
	r = [r[2], r[0], r[1], r[3]]
	A = 2 * r[0] * r[3] * cos(t[0]) - 2 * r[1] * r[3] * cos(t[1])
	B = 2 * r[0] * r[3] * sin(t[0]) - 2 * r[1] * r[3] * sin(t[1])
	C = r[0] ** 2 + r[1] ** 2 + r[3] ** 2 - r[2] ** 2 - 2 * r[0] * r[1] * (cos(t[0]) * cos(t[1]) + sin(t[0]) * sin(t[1]))
	x = symbols('x')
	z, y = solve((C - A)*x*x + 2*B*x + A + C, x)
	t32 = 2 * atan(y)
	t2 = atan2((r[0]*sin(t[0])+r[3]*sin(t32)-r[1]*sin(t[1]))/r[2], (r[0]*cos(t[0])+r[3]*cos(t32)-r[1]*cos(t[1]))/r[2])
	return float(t32), float(t2)


def fourbar(barnums=4):
	mstr1 = list()
	mstr2 = list()
	mstr3 = list()
	mstr4 = list()
	x = 0
	rp1 = list()
	rpp1 = list()
	tp1 = list()
	tpp1 = list()
	t1 = list()
	r1 = list()
	while x < barnums:
		rpz = symbols('rp' + str(x))
		rppz = symbols('rpp' + str(x))
		tpz = symbols('tp' + str(x))
		tppz = symbols('tpp' + str(x))
		rz = symbols('r' + str(x))
		tz = symbols('t' + str(x))
		r1.append(rz)
		t1.append(tz)
		rp1.append(rpz)
		tp1.append(tpz)
		rpp1.append(rppz)
		tpp1.append(tppz)
		x += 1
	vectr1 ,vectrp1 ,vectrpp1 = bargen(r1, rp1, rpp1, t1, tp1, tpp1, barnums)
	for x, y in enumerate(vectr1):
		print 'the r'+str(x)+' bar in I and J '+str(vectr1[x].i)+'i and '+str(vectr1[x].j)+'j'
		print 'the r'+str(x)+' dot bar in I and J '+str(vectrp1[x].i)+'i and '+str(vectrp1[x].j)+'j'
		print 'the r'+str(x)+' double dot bar in I and J '+str(vectrpp1[x].i)+'i and '+str(vectrpp1[x].j)+'j'
	x = 0
	while x <= 360:
		t1 = symbols('t1')
		t3 = symbols('t3')
		r = [1, 4, 5, 3]
		rp = [0, 0, 0, 0]
		rpp = [0, 0, 0, 0]
		t = [degtrad(x), t1, 0, t3]
		tp = [10, symbols('tp1'), 0, symbols('tp3')]
		tpp = [0, symbols('tpp1'), 0, symbols('tpp3')]
		a, s = twoangsolve1(r, t)
		t = [degtrad(x), s, 0, a]
		vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp, barnums)
		mstr1.append([vr[0].i + vr[1].i + 2 * cos(t[1])])
		mstr2.append([vr[0].j + vr[1].j + 2 * sin(t[1])])
		tp1, tp3 = veosolve1(r, t, rp, tp)
		tp = [10, tp1, 0, tp3]
		mstr3.append(tp)
		tpp1, tpp3 = accsolve1(r, t, rp, tp,rpp,tpp)
		tpp = [0, tpp1, 0, tpp3]
		mstr4.append(tpp)
		print x
		x += 6
	plt.figure(1)
	pl1 = plt.plot(mstr1, mstr2)
	plt.figure(2)
	pl2 = plt.plot(mstr3)
	plt.figure(3)
	pl3 = plt.plot(mstr4)
	plt.show()



# fourbar()

def q2():
	r = [(sin(degtrad(90+30))*12.5/sin(degtrad((180-90-30-30)))), 12.5, (sin(degtrad(30))*12.5/sin(degtrad((180-90-30-30))))]
	t = [degtrad(30), degtrad(0), degtrad(-120+180)]
	rp = [.5, 0, symbols('rp2')]
	rpp = [0,0,symbols('rpp2')]
	tp = [symbols('tp0'),0,0]
	tpp = [symbols('tpp1'),0,0]
	eq1=Eq(rp[0]*cos(t[0]) - r[0]*tp[0]*sin(t[0]), rp[2]*cos(t[2]) - r[2]*tp[2]*sin(t[2]))
	eq2=Eq(rp[0]*sin(t[0]) + r[0]*tp[0]*cos(t[0]), rp[2]*sin(t[2]) + r[2]*tp[2]*cos(t[2]))
	a,b=solve_linear(eq1)
	c,d=solve_linear(eq2)
	eq3=Eq(b,d)
	e,rp2=solve_linear(eq3)
	rp = [.5, 0, rp2]
	eq1 = Eq(rp[0] * cos(t[0]) - r[0] * tp[0] * sin(t[0]), rp[2] * cos(t[2]) - r[2] * tp[2] * sin(t[2]))
	a,tp0=solve_linear(eq1)
	tp = [tp0, 0, 0]
	eq1=Eq(rpp[0]*cos(t[0])+r[0]*tpp[0]*-sin(t[0])-r[0]*tp[0]**2*cos(t[0])+2*rp[0]*tp[0]*-sin(t[0]),rpp[2]*cos(t[2])+r[2]*tpp[2]*-sin(t[2])-r[2]*tp[2]**2*cos(t[2])+2*rp[2]*tp[0]*-sin(t[2]))
	eq2=Eq(rpp[0]*sin(t[0])+r[0]*tpp[0]*cos(t[0])-r[0]*tp[0]**2*sin(t[0])+2*rp[0]*tp[0]*cos(t[0]),rpp[2]*sin(t[2])+r[2]*tpp[2]*cos(t[2])-r[2]*tp[2]**2*sin(t[2])+2*rp[2]*tp[2]*cos(t[2]))
	a, b = solve_linear(eq1)
	c, d = solve_linear(eq2)
	eq3 = Eq(b, d)
	e, rpp2 = solve_linear(eq3)
	rpp = [0,0,rpp2]
	eq1=Eq(rpp[0]*cos(t[0])+r[0]*tpp[0]*-sin(t[0])-r[0]*tp[0]**2*cos(t[0])+2*rp[0]*tp[0]*-sin(t[0]),rpp[2]*cos(t[2])+r[2]*tpp[2]*-sin(t[2])-r[2]*tp[2]**2*cos(t[2])+2*rp[2]*tp[0]*-sin(t[2]))
	a,tpp0=solve_linear(eq1)
	tpp = [tpp0, 0, 0]
	print r
	print rp
	print rpp
	print t
	print tp
	print tpp


# q2()


def q4():
	t1=15*sin(degtrad(45))
	t2=asin((20-15*sin(degtrad(45)))/40)
	t3=40*cos(t2)+t1
	r=[15,40,t3,20]
	rp=[0,0,var('rp2'),0]
	rpp=[0,0,var('rpp2'),0]
	t=[degtrad(45),asin((20-15*sin(degtrad(45)))/40),0,pi/2]
	tp=[150,var('tp1'),0,0]
	tpp=[0,var('tpp1'),0,0]
	eq2=Eq(rp[0]*sin(t[0]) + r[0]*tp[0]*cos(t[0]) + rp[1]*sin(t[1]) + r[1]*tp[1]*cos(t[1]), rp[2]*sin(t[2]) + r[2]*tp[2]*cos(t[2]) + rp[3]*sin(t[3]) + r[3]*tp[3]*cos(t[3]))
	a,b=solve_linear(eq2)
	tp=[150,b,0,0]
	eq1=Eq(rp[0]*cos(t[0]) - r[0]*tp[0]*sin(t[0]) + rp[1]*cos(t[1]) - r[1]*tp[1]*sin(t[1]), rp[2]*cos(t[2]) - r[2]*tp[2]*sin(t[2]) + rp[3]*cos(t[3]) - r[3]*tp[3]*sin(t[3]))
	c,d=solve_linear(eq1)
	rp=[0,0,d,0]
	eq2=Eq(rpp[0]*sin(t[0])+r[0]*tpp[0]*cos(t[0])-r[0]*tp[0]**2*sin(t[0])+2*rp[0]*tp[0]*cos(t[0]) + rpp[1]*sin(t[1])+r[1]*tpp[1]*cos(t[1])-r[1]*tp[1]**2*sin(t[1])+2*rp[1]*tp[1]*cos(t[1]),rpp[2]*sin(t[2])+r[2]*tpp[2]*cos(t[2])-r[2]*tp[2]**2*sin(t[2])+2*rp[2]*tp[2]*cos(t[2]) + rpp[3]*sin(t[3])+r[3]*tpp[3]*cos(t[3])-r[3]*tp[3]**2*sin(t[3])+2*rp[3]*tp[3]*cos(t[3]))
	a,b=solve_linear(eq2)
	tpp=[0,b,0,0]
	eq1=Eq(rpp[0]*cos(t[0])+r[0]*tpp[0]*-sin(t[0])-r[0]*tp[0]**2*cos(t[0])+2*rp[0]*tp[0]*-sin(t[0]) + rpp[1]*cos(t[1])+r[1]*tpp[1]*-sin(t[1])-r[1]*tp[1]**2*cos(t[1])+2*rp[1]*tp[1]*-sin(t[1]),rpp[2]*cos(t[2])+r[2]*tpp[2]*-sin(t[2])-r[2]*tp[2]**2*cos(t[2])+2*rp[2]*tp[0]*-sin(t[2]) + rpp[3]*cos(t[3])+r[3]*tpp[3]*-sin(t[3])-r[3]*tp[3]**2*cos(t[3])+2*rp[3]*tp[3]*-sin(t[3]))
	c,d=solve_linear(eq1)
	rpp=[0,0,d,0]
	r1=sqrt(12**2+15**2)
	t1=atan(15/15)
	r=[15,r1]
	rp=[0,0]
	rpp=[0,0]
	t=[degtrad(45),t1]
	tp=[150,0]
	tpp=[0,0]
	vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,2)
	print vr[1].i
	print vr[1].j
	print vrp[1].i
	print vrp[1].j
	print vrpp[1].i
	print vrpp[1].j

# q4()

def q6():
	master1=list()
	master2=list()
	master3=list()
	master4=list()
	master5=list()
	master6=list()
	x=0
	while x<=360:
		r=[sqrt(6**2+3**2-2*6*3*cos(degtrad(90+x))),6,3]
		rp=[var('rp'),0,0]
		rpp=[var('rpp)'),0,0]
		if x>=90 and x<=270:
			t=[float(pi/2+acos((r[0]**2+6**2-3**2)/(2*r[0]*6))),pi/2,degtrad(x)+.0]
		else:
			t=[float(pi/2-acos((r[0]**2+6**2-3**2)/(2*r[0]*6))),pi/2,degtrad(x)+.0]
		tp=[var('tp0'),0,100]
		tpp=[var('tpp0'),0,0]
		# print t
		vr,vrp,vrpp=bargen(r,rp,rpp,t,tp,tpp,3)
		a,b=solve_linear(vrp[0].i,vrp[2].i)
		c,d=solve_linear(vrp[0].j,vrp[2].j)
		e,f=solve_linear(b,d)
		tp=[f,0,100]
		vr,vrp,vrpp=bargen(r,rp,rpp,t,tp,tpp,3)
		a,b=solve_linear(vrp[0].i,vrp[2].i)
		rp=[b,0,0]
		vr,vrp,vrpp=bargen(r,rp,rpp,t,tp,tpp,3)
		a,b=solve_linear(vrpp[0].i,vrpp[2].i)
		c,d=solve_linear(vrpp[0].j,vrpp[2].j)
		e,f=solve_linear(b,d)
		tpp=[f,0,0]
		vr,vrp,vrpp=bargen(r,rp,rpp,t,tp,tpp,3)
		a,b=solve_linear(vrpp[0].i,vrpp[2].i)
		rpp=[b,0,0]
		x+=5
		b=atan((12*sin(t[0])-6)/8)
		n=-cos(b)*8+cos(t[0])*12
		r1=[12,8,-n,6]
		rp1=[0,0,var('rp'),0]
		rpp1=[0,0,var('rpp'),0]
		t1=[t[0],pi+b,pi,pi/2]
		tp1=[tp[0],var('tp'),0,0]
		tpp1=[tpp[0],var('tpp'),0,0]
		vr1, vrp1, vrpp1 = bargen(r1, rp1, rpp1, t1, tp1, tpp1, 4)
		a,b=solve_linear(vrpp1[0].i+vrpp1[1].i,vrpp1[2].i+vrpp1[3].i)
		c,d=solve_linear(vrpp1[0].j+vrpp1[1].j,vrpp1[2].j+vrpp1[3].j)
		e,f=solve_linear(b,d)
		tp1=[tp[0],f,0,0]
		vr1, vrp1, vrpp1 = bargen(r1, rp1, rpp1, t1, tp1, tpp1, 4)
		a, b = solve_linear(vrpp1[0].i + vrpp1[1].i, vrpp1[2].i + vrpp1[3].i)
		rp1 = [0, 0, b, 0]
		vr1, vrp1, vrpp1 = bargen(r1, rp1, rpp1, t1, tp1, tpp1, 4)
		a,b=solve_linear(vrpp1[0].i+vrpp1[1].i,vrpp1[2].i+vrpp1[3].i)
		c,d=solve_linear(vrpp1[0].j+vrpp1[1].j,vrpp1[2].j+vrpp1[3].j)
		e,f=solve_linear(b,d)
		tpp1 = [tpp[0], f, 0, 0]
		vr1, vrp1, vrpp1 = bargen(r1, rp1, rpp1, t1, tp1, tpp1, 4)
		a, b = solve_linear(vrpp1[0].i + vrpp1[1].i, vrpp1[2].i + vrpp1[3].i)
		rpp1 = [0, 0, b, 0]
		vr1, vrp1, vrpp1 = bargen(r1, rp1, rpp1, t1, tp1, tpp1, 4)
		master1.append(vr1[0].i + vr1[1].i)
		master2.append(vr1[0].j + vr1[1].j)
		master3.append(vrp1[0].i + vrp1[1].i)
		master4.append(vrp1[0].j + vrp1[1].j)
		master5.append(vrpp1[0].i + vrpp1[1].i)
		master6.append(vrpp1[0].j + vrpp1[1].j)
	plt.figure(1)
	pl1 = plt.plot(master1, master2)
	plt.figure(2)
	pl2 = plt.plot(master3,master4)
	plt.figure(3)
	pl3 = plt.plot(master5,master6)
	plt.show()


q6()