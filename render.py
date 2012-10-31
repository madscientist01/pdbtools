#!/usr/bin/python
#
# Batch Pymol Renderer
# Require Python 2.7 (argparse) and Pymol
#
#

import sys, subprocess, re, os, glob
import argparse

def pdbparsing(filename, regex):
#
# Parse pdb content and extract using regex
#
	if os.path.exists(filename):
		pdb = open(filename)
		pdbcontent = pdb.readlines()
		matched = ''
		matchedcount=0
		reg = re.compile(regex)
		for line in pdbcontent:
			match = reg.match(line)
			if match:
				if matchedcount==0:
					gap = " "
				else:
					gap = ""
				matched = matched+gap+match.group(1).strip()
		return (matched)
	else:
		print "{0} is not found!".format(filename)
		sys.exit()


def htmlout(description):
#
# Generate render.html which display rendered png files
#
#
	htmlheader = """
<!DOCTYPE html>
<html>
	<head>
    	<title>List</title>
 		<meta charset='utf-8'>
	</head>
	<style>
	</style>
	<body>
"""
	imgtag = """
		<img src="./{0}">
		<p><a href="http://www.rcsb.org/pdb/explore/explore.do?structureId={1}">{2}</a></p>
"""

	htmlfooter = """
	</body>
</html>
"""
	htmloutput = open('render.html','w')
	htmloutput.write(htmlheader)
	for pdb,desc in description.items():
		htmloutput.write(imgtag.format(pdb[:-4]+".png", pdb[:len(pdb)-4], pdb[:len(pdb)-4]+":"+desc))
	htmloutput.write(htmlfooter)
	htmloutput.close()


def main(arglist):
#
# PyMol executable path
# In cases of MacPymol, It assumed that pymol is installed in 'Application' folder.
# Other cases, plase set up the exact location of PyMOL executable
#
	PyMOLPath = '/Applications/MacPyMOL.app/Contents/MacOS/MacPyMol'
	if len(arglist.file)>0:
		subjectList = arglist.file
	else:
		subjectList = glob.glob('*.pdb')
	if len(subjectList)>0 :
		pymolLoadingFile = open('render.pml','w')
		pymolLoadingFile.write('bg_color white\n')
		description = {}

		style = "cartoon" #Defaults
		if arglist.ribbon :
			style = "ribbon"
		if arglist.line :
			style = "line"
		if arglist.stick :
			style = "stick"

		if arglist.view:
			if os.path.exists(arglist.view):
				f = open(arglist.view)
				view = f.readlines()
				if re.match("^set_view",view[0]):
					for line in view:
						pymolLoadingFile.write(line)

		for pdb in subjectList:	
			if os.path.exists(pdb):			
				print pdb
				pymolLoadingFile.write("load {0}/{1}\n".format(os.getcwd(),pdb))
				pymolLoadingFile.write("hide all\n")
				pymolLoadingFile.write("show {0}, {1}\n".format(style, pdb[:-4]))
				if arglist.cbc:
					pymolLoadingFile.write("util.cbc\n")
				else:
					pymolLoadingFile.write("spectrum count,rainbow,{0}\n".format(pdb[:-4]))
				if arglist.surface:
					pymolLoadingFile.write("show surface,{0}\n".format(pdb[:-4]))
				if not arglist.view:	
					pymolLoadingFile.write("orient {0}\n".format(pdb[:-4]))
				if arglist.ray:
					pymolLoadingFile.write("ray\n")
				pymolLoadingFile.write("png {0}.png\n".format(pdb[:-4]))
				description[pdb]=pdbparsing(pdb, '^TITLE    [ \d](.*)$')
		pymolLoadingFile.close()
		print "Running PyMol..."
		p = subprocess.Popen([PyMOLPath,'-c','render.pml'],
							  stdout=subprocess.PIPE)
		p_stdout = p.stdout.read()
		htmlout(description)
	else:
		print "Usage : render.py"
		sys.exit()		


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--ribbon', action='store_true', dest='ribbon',default=False,
	                    help='Draw as ribbon')
	parser.add_argument('-l', '--line', action='store_true', dest='line',default=False,
	                    help='Draw as Line')
	parser.add_argument('-t', '--stick', action='store_true', dest='stick',default=False,
	                    help='Draw as Stick')
	parser.add_argument('-s', '--surface', action='store_true', dest='surface',default=False,
	                    help='Draw as Surface')
	parser.add_argument('-c', '--color_by_chain', action='store_true', dest='cbc',default=False,
	                    help='Color by Chain')
	parser.add_argument('-rt', '--ray_trace', action='store_true', dest='ray',default=False,
	                    help='Ray tracing')
	parser.add_argument('-f', '--files', nargs='*', dest='file',default=[],
	                    help='File to process')
	parser.add_argument('-v', '--fixed_view', action='store', dest='view')
	results = parser.parse_args()
	main(results)
