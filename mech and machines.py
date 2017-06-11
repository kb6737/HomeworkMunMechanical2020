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


def rev2rad(rev):
	return rev*0.1047


def interpolation(known1, known2, unknown1, known3, known4):
	unknown2 = known3+((unknown1-known1)*(known4-known3)/(known2-known1))
	return unknown2


def question1():
	wab = 2
	lab = 75
	lbc = interpolation(2, 4, 3 ,75, 100)
	lcd = 100
	lad = math.sqrt(math.pow(50, 2)+math.pow(250, 2))
	print locals()


def question4():
	woc = rev2rad(120)
	loc = 125
	lcp = 500
	lpa = 125
	laq = 250
	lqe = 125


def question6():
	woa = rev2rad(30)
	loa = 15
	loq = 10
	lqb = 12.5
	lbc = 50


def question8():
	woa = rev2rad(130)
	loa = 100
	lab = 400
	lce = 400
	lac = 150
	lef = 300


question1()