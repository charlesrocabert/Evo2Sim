
#!/usr/bin/env python
# coding: utf-8

#***************************************************************************
# Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon
# Web: https://github.com/charlesrocabert/Evo2Sim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***************************************************************************

import sys
import os

#-------------#
# PRINT USAGE #
#-------------#
def print_usage():
	print ""
	print "=== CHANGE VERSION ==="
	print "Usage: python change_version.py [parameters]"
	print "Parameters are:"
	print "-h, --help:"
	print "    Print this help, then exit."
	print "-p, --path:"
	print "    Give the path of the project to modify."
	print "-vmajor, --version-major:"
	print "    Give the new major version."
	print "-vminor, --version-minor:"
	print "    Give the new minor version."
	print "-vpatch, --version-patch:"
	print "    Give the new patch version."
	print ""

#-----------------------------#
# READ COMMAND LINE ARGUMENTS #
#-----------------------------#
def read_args( args ):
	if len(args) < 2:
		print("Lack of parameters, see help (-h --help)")
		sys.exit()
	else:
		PACKAGE = ""
		MAJOR   = ""
		MINOR   = ""
		PATCH   = ""
		for i in range(len(args)):
			if args[i] == "-h" or args[i] == "--help":
				print_usage()
				sys.exit()
			elif args[i] == "-p" or args[i] == "--path":
				if i+1 >= len(args):
					print("Lack of parameters, see help (-h --help)")
					sys.exit()
				else:
					PATH = args[i+1]
			elif args[i] == "-vmajor" or args[i] == "--version-major":
				if i+1 >= len(args):
					print("Lack of parameters, see help (-h --help)")
					sys.exit()
				else:
					MAJOR = args[i+1]
			elif args[i] == "-vminor" or args[i] == "--version-minor":
				if i+1 >= len(args):
					print("Lack of parameters, see help (-h --help)")
					sys.exit()
				else:
					MINOR = args[i+1]
			elif args[i] == "-vpatch" or args[i] == "--version-patch":
				if i+1 >= len(args):
					print("Lack of parameters, see help (-h --help)")
					sys.exit()
				else:
					PATCH = args[i+1]
		if PATH == "" or MAJOR == "" or MINOR == "" or PATCH == "":
			print("Lack of parameters, see help (-h --help)")
			sys.exit()
		else:
			return PATH, MAJOR, MINOR, PATCH

#-----------------------------#
# MODIFY CMakeLists.txt       #
#-----------------------------#
def modify_CMakeLists( path, major, minor, patch ):
	f = open(path+"CMakeLists.txt", "r")
	data = ""
	l = f.readline()
	while l:

		if l.startswith("set(PACKAGE"):
			data += "set(PACKAGE \"\\\""+"Evo2Sim"+"\\\"\")\n"

		elif l.startswith("set(VERSION_MAJOR"):
			data += "set(VERSION_MAJOR "+str(major)+")\n"

		elif l.startswith("set(VERSION_MINOR"):
			data += "set(VERSION_MINOR "+str(minor)+")\n"

		elif l.startswith("set(VERSION_PATCH"):
			data += "set(VERSION_PATCH "+str(patch)+")\n"

		else:
			data += l

		l = f.readline()
	f.close()
	f = open(path+"CMakeLists.txt", "w")
	f.write(data)
	f.close()

#-----------------------------#
# MODIFY LICENSE.md           #
#-----------------------------#
def modify_license( path, major, minor, patch ):
	version = major+"."+minor+"."+patch
	f = open(path+"LICENSE.md", "w")
	f.write("\n# Evo<sup>2</sup>Sim "+version+"\n")
	f.write("\n")
	f.write("Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon.\n")
	f.write("All rights reserved\n")
	f.write("\n")
	f.write("This program is free software: you can redistribute it and/or modify\n")
	f.write("it under the terms of the GNU General Public License as published by\n")
	f.write("the Free Software Foundation, either version 3 of the License, or\n")
	f.write("(at your option) any later version.\n")
	f.write("\n")
	f.write("This program is distributed in the hope that it will be useful,\n")
	f.write("but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
	f.write("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
	f.write("GNU General Public License for more details.\n")
	f.write("\n")
	f.write("You should have received a copy of the GNU General Public License\n")
	f.write("along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\n")
	f.close()

#-----------------------------#
# MODIFY README.md            #
#-----------------------------#
def modify_readme( path, major, minor, patch ):
	f = open(path+"README.md", "r")
	data = "\n# Evo<sup>2</sup>Sim "+str(major)+"."+str(minor)+"."+str(patch)+" <br />&nbsp;&nbsp;<img src=\"logo/logo_evo2sim_small.png\"><br />\n"
	l = f.readline()
	l = f.readline()
	while l:
		data += l
		l = f.readline()
	f.close()
	f = open(path+"README.md", "w")
	f.write(data)
	f.close()

##################
#      MAIN      #
##################
if __name__ == '__main__':
	PATH, MAJOR, MINOR, PATCH = read_args(sys.argv)
	modify_CMakeLists(PATH, MAJOR, MINOR, PATCH)
	modify_license(PATH, MAJOR, MINOR, PATCH)
	modify_readme(PATH, MAJOR, MINOR, PATCH)



