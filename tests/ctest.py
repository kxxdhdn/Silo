#!/usr/bin/env python
# -*- coding: utf-8 -*-

class parent:
	def __init__(self, x, y):
		self.x = x

		z0 = x + y
		z = self.x + y
		self.z = self.x + y + 1.
		
		print('z0, z, self.z', z0, z, self.z)
		print('>> coucou 0 <<')

	def P2(self):
		self.P1(2., 4.)
		print('>> coucou 2 <<')

	def P1(self, x1, y1):
		print('self.x', self.x)
		self.x = 100.
		print('new self.x', self.x)
		print('>> coucou 1 <<')

		return x1 * y1


class child(parent):
	def __init__(self, x, y):
		super().__init__(x, y)

		print('self.x in child', self.x)
		self.x = 101.
		print('new self.x in child', self.x)
		self.a = self.P1(20., 40.)
		print('self.P1(20., 40.) in child', self.a)
		self.P2()
		print('>> coucou 00 <<')
		
