import unittest
import cStringIO
import sys
import os

direct=os.getcwd()
sys.path.append(direct)

from assm_align import mergeLoc, getLength

class check_mergeLoc(unittest.TestCase):

	def setUp(self):
		self.loc_list1=['1:1-100', '1:50-150']
		self.merge_loc_list1=[('1',1,150, 149)]
		self.loc_list2=['1:1-100', '1:110-120']
		self.merge_loc_list2=[('1',1,100,99), ('1',110,120,10)]
		self.loc_list3=['1:1-100', '1:50-150', '2:1-100']
		self.merge_loc_list3=[('1',1,150,149), ('2',1,100,99)]

	def test_mergeLoc1(self):
		self.assertEqual(mergeLoc(self.loc_list1), self.merge_loc_list1)

	def test_mergeLoc2(self):
		self.assertEqual(mergeLoc(self.loc_list2), self.merge_loc_list2)

	def test_mergeLoc3(self):
		self.assertEqual(mergeLoc(self.loc_list3), self.merge_loc_list3)

class check_getLength(unittest.TestCase):

	def setUp(self):
		self.loc_list1=[('1',1,150, 149)]
		self.len1={'1': 149}
		self.loc_list2=[('1',1,150,149), ('2',1,100,99)]
		self.len2={'1': 149, '2': 99}

	def test_getLen1(self):
		self.assertEqual(getLength(self.loc_list1), self.len1)

	def test_getLen2(self):
		self.assertEqual(getLength(self.loc_list2), self.len2)

if __name__ == '__main__':
	unittest.main()
