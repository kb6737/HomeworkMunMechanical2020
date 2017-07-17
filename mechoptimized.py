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


def twoangsolve1(r, t):
	t = [t[2], t[0], t[1], t[3]]
	r = [r[2], r[0], r[1], r[3]]
	A = 2 * r[0] * r[3] * cos(t[0]) - 2 * r[1] * r[3] * cos(t[1])
	B = 2 * r[0] * r[3] * sin(t[0]) - 2 * r[1] * r[3] * sin(t[1])
	C = r[0] ** 2 + r[1] ** 2 + r[3] ** 2 - r[2] ** 2 - 2 * r[0] * r[1] * (cos(t[0]) * cos(t[1]) + sin(t[0]) * sin(t[1]))
	x = symbols('x')
	z, y = solve((C - A)*x*x + 2*B*x + A + C, x)
	t32 = 2 * atan(y)
	t1 = atan2((r[0]*sin(t[0])+r[3]*sin(t32)-r[1]*sin(t[1]))/r[2], (r[0]*cos(t[0])+r[3]*cos(t32)-r[1]*cos(t[1]))/r[2])
	return float(t32), float(t1)


def q1():
    master1 = list()
    master2 = list()
    master3 = list()
    master4 = list()
    master5 = list()
    master6 = list()
    x=0
    while x<=360:
        r = [1, 4, 5, 3]
        rp = [0, 0, 0, 0]
        rpp = [0, 0, 0, 0]
        t = [degtrad(x), var('t1'), 0, var('t3')]
        tp = [10, var('tp1'), 0, var('tp3')]
        tpp = [0, var('tpp1'), 0, var('tpp3')]
        t3, t1 = twoangsolve1(r,t)
        t = [degtrad(x),t1,0,t3]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp)
        master1.append([vr[0].i + vr[1].i + 2 * cos(t[1])])
        master2.append([vr[0].j + vr[1].j + 2 * sin(t[1])])
        eq1 = Eq(vrp[0].i+vrp[1].i,vrp[2].i+vrp[3].i)
        eq2 = Eq(vrp[0].j+vrp[1].j,vrp[2].j+vrp[3].j)
        tp3, tp1 = solve_linear(eq1)
        tp32, tp12 = solve_linear(eq2)
        V, tp1 = solve_linear(tp1, tp12)
        tp = [10, tp1, 0, var('tp3')]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp)
        eq1 = Eq(vrp[0].i+vrp[1].i,vrp[2].i+vrp[3].i)
        tp3, value = solve_linear(eq1)
        tp = [10, tp1, 0, value]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp)
        eq1 = Eq(vrpp[0].i+vrpp[1].i,vrpp[2].i+vrpp[3].i)
        eq2 = Eq(vrpp[0].j+vrpp[1].j,vrpp[2].j+vrpp[3].j)
        tp3, tp1 = solve_linear(eq1)
        tp32, tp12 = solve_linear(eq2)
        V, tp1 = solve_linear(tp1, tp12)
        tpp = [0, tp1, 0, var('tpp3')]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp)
        eq1 = Eq(vrpp[0].i+vrpp[1].i,vrpp[2].i+vrpp[3].i)
        tp3, value = solve_linear(eq1)
        tpp = [0, tp1, 0, value]
	r = [1, 6, 5, 3]
	vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp)
        master3.append(vrp[0].i+vrp[1].i)
	master4.append(vrp[0].j+vrp[1].j)
        master5.append(vrpp[0].i+vrpp[1].i)
	master6.append(vrpp[0].j+vrpp[1].j)
        x += 15
    for x, y in enumerate(vr):
        print 'the r'+str(x)+' bar in I and J '+str(vr[x].i)+'i and '+str(vr[x].j)+'j'
        print 'the r'+str(x)+' dot bar in I and J '+str(vrp[x].i)+'i and '+str(vrp[x].j)+'j'
        print 'the r'+str(x)+' double dot bar in I and J '+str(vrpp[x].i)+'i and '+str(vrpp[x].j)+'j'
        print 'a vector loop equation is an isolation of a variable knowing two sides equal zero.'
    r = [var('r0'),var('C'),var('r2'),var('r3')]
    rp = [var('rp0'),var('Cp'),var('rp2'),var('rp3')]
    rpp = [var('rpp0'),var('Cpp'),var('rpp2'),var('rpp3')]
    t = [var('t0'),var('tc'),var('t2'),var('t3')]
    tp = [var('tp0'),var('tpc'),var('tp2'),var('tp3')]
    tpp = [var('tpp0'),var('tppc'),var('tpp2'),var('tpp3')]
    vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp)
    print 'vrp[0].i + vrp[1].i = vrp[2].i + vrp[3].i'
    print 'where "vrp" stands for vector r dot, we start at the r dot becasue the r values can be solved geometricaly'
    print 'the general form for C will be'
    print str(vr[0].i+vr[1].i,vr[2].i+vr[3].i)
    print str(vr[0].j+vr[1].j,vr[2].j+vr[3].j)
    print str(vrp[0].i+vrp[1].i,vrp[2].i+vrp[3].i)
    print str(vrp[0].j+vrp[1].j,vrp[2].j+vrp[3].j)
    print str(vrpp[0].i+vrpp[1].i,vrpp[2].i+vrpp[3].i)
    print str(vrpp[0].j+vrpp[1].j,vrpp[2].j+vrpp[3].j)
    plt.figure(1)
    plt.plot(master1,master2)
    plt.figure(2)
    plt.plot(master3,master4)
    plt.figure(3)
    plt.plot(master5,master6)
    plt.show()


def q2():
    L, value = solve_linear(var('L')/sin(degtrad(90+30)),12.5/sin(degtrad(180-20-120)))
    L, value1 = solve_linear(value/sin(degtrad(90+30)),var('L')/sin(degtrad(20)))
    r = [12.5,value1,value]
    t = [degtrad(0), degtrad(60), degtrad(20)]
    rp = [0, 0.5, var('rp1')]
    rpp = [0, 0, var('rpp1')]
    tp = [0, 0, symbols('tp0')]
    tpp = [0, 0, symbols('tpp2')]
    vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
    eq1 = Eq(vrp[0].i+vrp[1].i,vrp[2].i)
    eq2 = Eq(vrp[0].j+vrp[1].j,vrp[2].j)
    a, b = solve_linear(eq1)
    c, d = solve_linear(eq2)
    e, f = solve_linear(b,d)
    rp = [0, 0.5, f]
    vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
    eq1 = Eq(vrp[0].i+vrp[1].i,vrp[2].i)
    a, b = solve_linear(eq1)
    tp = [0, 0, b]
    vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
    eq1 = Eq(vrpp[0].i+vrpp[1].i,vrpp[2].i)
    eq2 = Eq(vrpp[0].j+vrpp[1].j,vrpp[2].j)
    a, b = solve_linear(eq1)
    c, d = solve_linear(eq2)
    e, f = solve_linear(b,d)
    rpp = [0, 0, f]
    vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
    eq1 = Eq(vrpp[0].i+vrpp[1].i,vrpp[2].i)
    a, b = solve_linear(eq1)
    tpp = [0, 0, b]
    print r
    print rp
    print rpp
    print t
    print tp
    print tpp
    r = [var('r0'),var('C'),var('r2')]
    rp = [var('rp0'),var('Cp'),var('rp2')]
    rpp = [var('rpp0'),var('Cpp'),var('rpp2')]
    t = [var('t0'),var('tc'),var('t2')]
    tp = [var('tp0'),var('tpc'),var('tp2')]
    tpp = [var('tpp0'),var('tppc'),var('tpp2')]
    vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)    
    print str(vr[0].i+vr[1].i,vr[2].i)
    print str(vr[0].j+vr[1].j,vr[2].j)
    print str(vrp[0].i+vrp[1].i,vrp[2].i)
    print str(vrp[0].j+vrp[1].j,vrp[2].j)
    print str(vrpp[0].i+vrpp[1].i,vrpp[2].i)
    print str(vrpp[0].j+vrpp[1].j,vrpp[2].j)
    print 'this can be used to solve for any problem given any loop equation but to solve for the inital lengths and angles'
    print 'we use sin law for this and generate the equation ('L')/sin((90+30)) = 12.5/sin((180-x-120)) where x is any angle limited by geometry'

def q4():
    master1 = list()
    master2 = list()
    master3 = list()
    master4 = list()
    x=45
    while x==45:
        valuevar=15*sin(degtrad(x))-20
        s = sqrt(40**2-valuevar**2)
        value=15*cos(degtrad(x))+s
        r = [15, 40, value, 20]
        rp = [0, 0, var('rp2'), 0]
        rpp = [0, 0, var('rpp2'), 0]
        t = [degtrad(x), atan(valuevar/s), 0, pi / 2]
        tp = [15.71, var('tp1'), 0, 0]
        tpp = [0, var('tpp1'), 0, 0]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,4)
        eq1 = Eq(vrp[0].i+vrp[1].i,vrp[2].i+vrp[3].i)
        eq2 = Eq(vrp[0].j+vrp[1].j,vrp[2].j+vrp[3].i)
        a, b = solve_linear(eq1)
        c, d = solve_linear(eq2)
        e, f = solve_linear(b,d)
        rp = [0, 0, f, 0]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,4)
        eq1 = Eq(vrp[0].i+vrp[1].i,vrp[2].i+vrp[3].i)
        a, b = solve_linear(eq1)
        tp = [150, b, 0, 0]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,4)
        eq1 = Eq(vrpp[0].i+vrpp[1].i,vrpp[2].i+vrpp[3].i)
        eq2 = Eq(vrpp[0].j+vrpp[1].j,vrpp[2].j+vrpp[3].j)
        a, b = solve_linear(eq1)
        c, d = solve_linear(eq2)
        e, f = solve_linear(b,d)
        rpp = [0, 0, f, 0]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,4)
        eq1 = Eq(vrpp[0].i+vrpp[1].i,vrpp[2].i+vrpp[3].i)
        a, b = solve_linear(eq1)
        tpp = [0, b, 0, 0]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,4)
        master1.append(vr[0].i+(15*cos(t[1])**2+15*cos(t[1]+pi/2)))
        master2.append(vr[0].j+(15*sin(t[1])**2+15*sin(t[1]+pi/2)))
	print 'the initial mechanism is solved then the point C is solved knowing the difference in geometry from the first
        print 'then the componets of C give us a 3 bar'
	r1 = [15,15,sqrt(15**2+15**2]
        rp1 = [0,0,0]
        rpp1 = [0, 0, 0]
        t1 = [t[1],pi/2+t[1],t[1]+degtrad(15)]
        tp1 = [tp[1],0,0]
        tpp1 = [tpp[1],0,0]
        vr1, vrp1, vrpp1 = bargen(r1, rp1, rpp1, t1, tp1, tpp1,3)
	print 'C velocity in i '+str(vrp[0].i+vrp1[2].i)+' and j '+str(vrp[0].j+vrp1[2].j)
        print 'C acceleration in i '+str(vrpp[0].i+vrpp1[2].i)+' and j '+str(vrpp[0].j+vrpp1[2].j)
	x+=15
    plt.figure(4)
    pl1 = plt.plot(master1, master2)
    plt.show()
    

def q6():
    master1 = list()
    master2 = list()
    master3 = list()
    master4 = list()
    x=0
    while x<=360:
        r = [6, 3,var('r2')]
        rp = [0, 0, var('rp2')]
        rpp = [0, 0, var('rpp2')]
        if x<=90 and x>=270:
            t = [pi/2,degtrad(x),(pi/2)-(atan((cos(degtrad(x))*r[1])/r[0]))]
        else:
            t = [pi/2,degtrad(x),(pi/2)+(atan((cos(degtrad(x))*r[1])/r[0]))]
        tp = [0, 100, var('tp2')]
        tpp = [0, 0, var('tpp2')]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
        eq1 = Eq(vr[0].i+vr[1].i,vr[2].i)
        a, b = solve_linear(eq1)
        r = [6, 3, b]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
        eq1 = Eq(vrp[0].i + vrp[1].i, vrp[2].i)
        eq2 = Eq(vrp[0].j + vrp[1].j, vrp[2].j)
        a, b = solve_linear(eq1)
        c, d = solve_linear(eq2)
        e, f = solve_linear(b, d)
        tp = [0, 100, f]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
        eq1 = Eq(vrp[0].i + vrp[1].i, vrp[2].i)
        a, b = solve_linear(eq1)
        rp = [0, 0, b]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp,3)
        eq1 = Eq(vrpp[0].i + vrpp[1].i, vrpp[2].i)
        eq2 = Eq(vrpp[0].j + vrpp[1].j, vrpp[2].j)
        a, b = solve_linear(eq1)
        c, d = solve_linear(eq2)
        e, f = solve_linear(b, d)
        rpp = [0, 0, f]
        vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp, 3)
        eq1 = Eq(vrpp[0].i + vrpp[1].i, vrpp[2].i)
        a, b = solve_linear(eq1)
        tpp = [0, 0, b]
        x+=15
        h=cos(t[0])*12 - 6
        r1 = [cos(t[0])*12+8**2-h**2,6,12,8]
        t1 = [pi,pi/2,t[2],asin(h/8)+pi]
        rp1 = [var('rp0'),0,0,0]
        rpp1 = [var('rpp0'),0,0,0]
        tp1 = [0,0,tp[2],var('tp3')]
        tpp1 = [0,0,tpp[2],var('tpp3')]
        vr, vrp, vrpp = bargen(r1, rp1, rpp1, t1, tp1, tpp1)
        eq1 = Eq(vrp[0].i + vrp[1].i, vrp[2].i + vrp[3].i)
        eq2 = Eq(vrp[0].j + vrp[1].j, vrp[2].j + vrp[3].j)
        a, b = solve_linear(eq1)
        c, d = solve_linear(eq2)
        e, f = solve_linear(b, d)
        rp1 = [f, 0, 0, 0]
        vr, vrp, vrpp = bargen(r1, rp1, rpp1, t1, tp1, tpp1)
        eq1 = Eq(vrp[0].i + vrp[1].i, vrp[2].i + vrp[3].i)
        a, b = solve_linear(eq1)
        tp1 = [0,0,tp[2],b]
        vr, vrp, vrpp = bargen(r1, rp1, rpp1, t1, tp1, tpp1)
        eq1 = Eq(vrpp[0].i + vrpp[1].i, vrpp[2].i + vrpp[3].i)
        eq2 = Eq(vrpp[0].j + vrpp[1].j, vrpp[2].j + vrpp[3].j)
        a, b = solve_linear(eq1)
        c, d = solve_linear(eq2)
        e, f = solve_linear(b, d)
        tpp1 = [0,0,tpp[2],f]
        vr, vrp, vrpp = bargen(r1, rp1, rpp1, t1, tp1, tpp1)
        eq1 = Eq(vrpp[0].i + vrpp[1].i, vrpp[2].i + vrpp[3].i)
        a, b = solve_linear(eq1)
        rpp1 = [b,0,0,0]
        vr, vrp, vrpp = bargen(r1, rp1, rpp1, t1, tp1, tpp1)
        master1.append(vrp[0].i + vrp[1].i)
        master2.append(vrp[0].j + vrp[1].j)
        master3.append(vrpp[0].i + vrpp[1].i)
        master4.append(vrpp[0].j + vrpp[1].j)
	    for x, y in enumerate(vr):
                print 'the r'+str(x)+' bar in I and J '+str(vr[x].i)+'i and '+str(vr[x].j)+'j'
                print 'the r'+str(x)+' dot bar in I and J '+str(vrp[x].i)+'i and '+str(vrp[x].j)+'j'
                print 'the r'+str(x)+' double dot bar in I and J '+str(vrpp[x].i)+'i and '+str(vrpp[x].j)+'j'
                print 'a vector loop equation is an isolation of a variable knowing two sides equal zero.'
    r = [var('r0'),var('r'),var('r2'),var('r3')]
    rp = [var('rp0'),var('rp'),var('rp2'),var('rp3')]
    rpp = [var('rpp0'),var('rpp1'),var('rpp2'),var('rpp3')]
    t = [var('t0'),var('tp'),var('t2'),var('t3')]
    tp = [var('tp0'),var('tp1'),var('tp2'),var('tp3')]
    tpp = [var('tpp0'),var('tpp1'),var('tpp2'),var('tpp3')]
    vr, vrp, vrpp = bargen(r, rp, rpp, t, tp, tpp)
    print 'vrp[0].i + vrp[1].i = vrp[2].i + vrp[3].i'
    print 'where "vrp" stands for vector r dot, we start at the r dot becasue the r values can be solved geometricaly'
    print 'the general form for C will be'
    print str(vr[0].i+vr[1].i,vr[2].i+vr[3].i)
    print str(vr[0].j+vr[1].j,vr[2].j+vr[3].j)
    print str(vrp[0].i+vrp[1].i,vrp[2].i+vrp[3].i)
    print str(vrp[0].j+vrp[1].j,vrp[2].j+vrp[3].j)
    print str(vrpp[0].i+vrpp[1].i,vrpp[2].i+vrpp[3].i)
    print str(vrpp[0].j+vrpp[1].j,vrpp[2].j+vrpp[3].j)
    print 'to obtain C we need the loop QAO solved. That is done geometrically using by using a position analysis or sin law'
    print 'Qa=sqrt(6^2+(3*cos(x))^2)'
    print 'we then fill in our data and solve for what we need to find (this data is at 360 degreees)'
    r1 = [cos(t[0])*12+8**2-h**2,6,12,8]
    t1 = [pi,pi/2,t[2],asin(h/8)+pi]
    rp1 = [var('rp0'),0,0,0]
    rpp1 = [var('rpp0'),0,0,0]
    tp1 = [0,0,tp[2],var('tp3')]
    tpp1 = [0,0,tpp[2],var('tpp3')]
    print r1
    print rp1    
    print rpp1
    print t1
    print tp1
    print tpp1
    print 'the i components are equal to our rp0 and rpp0 values'
    plt.figure(5)
    pl1 = plt.plot(master1, master2)
    plt.figure(6)
    pl1 = plt.plot(master3, master4)
    plt.show()


q1()
q2()
q4()
q6()
