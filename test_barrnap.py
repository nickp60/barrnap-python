# -*- coding: utf-8 -*-
"""
@author: nicholas

"""
import sys
import logging
import shutil
import os
import unittest
import time
from argparse import Namespace

from unittest.mock import Mock

from barrnap import main

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class BarrnapCase(unittest.TestCase):
    """ tests for barrnap.py
    """
    def setUp(self):
        self.startTime = time.time()

    def test_main_return_1_invalid(self):
        args = Namespace(threads = -1)
        with self.assertRaises(SystemExit):
            main(args)

    def tearDown(self):
        """
        """
        pass


if __name__ == '__main__':
    unittest.main()
