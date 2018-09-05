## this enables the ACSESS package to be run directly
import sys, os
sys.path.append('.')

# just import acsess doesn't run as __main__
path = os.path.dirname( os.path.realpath(__file__))
print "path:", path
os.system("python {}/acsess.py".format(path))
