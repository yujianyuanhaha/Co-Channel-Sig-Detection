import unittest
import os
import subprocess
import sys


# To be compatible with the the standard Python unittest framework, there must be a class
# defined which derives from unittest.TestCase.
# The class name (here, ExampleTests) need not match the filename, but it is recommended
# for clarity.
class IM_Test(unittest.TestCase):
  # To define a test, simply define a class method whose name is prefixed with "test".
  # It's recommended that what follows after "test" be descriptive of the test being run.

#  def test_example_compliant_function(self):
#    # This is the preferred method for executing a test function.  The function should be designed such that it 
#    # returns as its first (or only) return value an integer representing success (0) or failure (non-zero).
#    test_script = "exit(framework_compliant_function)"
#
#    # Boilerplate.  Do not modify.
#    self.assertEqual(0, subprocess.run('matlab -r "try, ' + test_script + ', catch ME, disp(ME); exit(1), end"', 
#                        shell=True, stdout=sys.stdout, stderr=sys.stderr, cwd=os.path.dirname(__file__)).returncode)

  def test_NBI_DSS_script(self):
    
    test_script = "run './NBI - DSSS/TestDSSS'; if STATUS == 'PASSED', exit(0), end"
    # Boilerplate.  Do not modify.
    self.assertEqual(0, subprocess.run('matlab -r "try, ' + test_script + ', catch ME, disp(ME); exit(1), end"', 
                        shell=True, stdout=sys.stdout, stderr=sys.stderr, cwd=os.path.dirname(__file__)).returncode)
						
  def test_NBI_FreqDomain_script(self):
    
    test_script = "run './NBI - Frequency Domain Nonlinear/TestNbiFreqDomain'; if STATUS == 'PASSED', exit(0), end"
    # Boilerplate.  Do not modify.
    self.assertEqual(0, subprocess.run('matlab -r "try, ' + test_script + ', catch ME, disp(ME); exit(1), end"', 
                        shell=True, stdout=sys.stdout, stderr=sys.stderr, cwd=os.path.dirname(__file__)).returncode)
						
  def test_NBI_LinearCancel_script(self):
    
    test_script = "run './NBI - Time Domain Linear/TestNbiLinearCancel'; if STATUS == 'PASSED', exit(0), end"
    # Boilerplate.  Do not modify.
    self.assertEqual(0, subprocess.run('matlab -r "try, ' + test_script + ', catch ME, disp(ME); exit(1), end"', 
                        shell=True, stdout=sys.stdout, stderr=sys.stderr, cwd=os.path.dirname(__file__)).returncode)
						
  def test_NBI_NotchFilter_script(self):
    
    test_script = "run './NBI - Time Domain Linear/TestNotchFilter'; if STATUS == 'PASSED', exit(0), end"
    # Boilerplate.  Do not modify.
    self.assertEqual(0, subprocess.run('matlab -r "try, ' + test_script + ', catch ME, disp(ME); exit(1), end"', 
                        shell=True, stdout=sys.stdout, stderr=sys.stderr, cwd=os.path.dirname(__file__)).returncode)
