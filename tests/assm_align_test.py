#!/usr/bin/env python

import unittest
import cStringIO
import sys
import os

direct=os.getcwd()
print direct
sys.path.append("%s/assm_align" % direct)

from assm_align import mergeLoc, getLength

class check_mergeLoc(unittest.TestCase):

	def setUp(self):
		self.loc_list1=['1:1-100', '1:50-150']
		self.merge_loc_list1=['1:1-150']

	def test_mergeLoc1(self):
		self.assertEqual(mergeLoc(self.loc_list1), self.merge_loc_list1)


if __name__ == '__main__':
	unittest.main()