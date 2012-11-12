#!/usr/bin/python
#
# Batch Pymol Renderer
# Require Python 2.7 (argparse) and Pymol
#
#

import sys, subprocess, re, os, glob
import argparse

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
		div.head {
		font-family: Sans-Serif;
		font-size: 14px;
		border:3px solid #CCCCCC;
		border-radius: 10px;
		padding: 10px;
	    align :center;
		background-color: #CCEEFF;
		}
		
	</style>
	<body>
	<table>
		<tr>
"""
	imgtag = """
		<tr>
			<td>
				<img src="./{0}">		
			</td>
			<td class="content">
			<div class="head">
			<p><a href="http://www.rcsb.org/pdb/explore/explore.do?structureId={1}">{1}</a></p>
			<p>{2}</p>
			<p>Resolution:{3}</p>
			<p>Authors:{4}</p>
			</div>
			</td>
		</tr>
"""

	htmlfooter = """

	</table>
	</body>
</html>
"""
	htmloutput = open('render.html','w')
	htmloutput.write(htmlheader)
	for pdb,desc in description.items():
		htmloutput.write(imgtag.format(desc.filename[:len(desc.filename)-4]+".png", desc.pdb, desc.title, desc.resolution, desc.author))
	htmloutput.write(htmlfooter)
	htmloutput.close()

class PDBdescription():

	def __init__(self):
		self.title = ""
		self.resolution = ""
		self.filename = ""
		self.pdb = ""
		self.author = ""

	def parse(self,filename):

		regexs=[]
		regexs.append('^TITLE    [ \d](.*)$')
		regexs.append('^AUTHOR    (.*)$')
		regexs.append('^REMARK   2 RESOLUTION.\s+([.\d]+) ANGSTROMS')
		regexs.append('^DBREF  (\S{4}) A')
		parsed = self.pdbparsing(filename,regexs)
		self.title = parsed['^TITLE    [ \d](.*)$']
		self.resolution = parsed['^REMARK   2 RESOLUTION.\s+([.\d]+) ANGSTROMS']
		self.pdb = parsed['^DBREF  (\S{4}) A']
		self.filename = filename
		self.author = parsed['^AUTHOR    (.*)$']

	def pdbparsing(self, filename,regexs):
#
# Parse pdb content and extract using regex
#
		if os.path.exists(filename):
			pdb = open(filename)
			pdbcontent = pdb.readlines()
			pdb.close()
			matched = {}
			matchedcounts={}
			compiledRegex = {}
			for regex in regexs:
				compiledRegex[regex]=re.compile(regex)
				matchedcounts[regex]=0
				matched[regex]=''

			for line in pdbcontent:

				for regex in compiledRegex.iterkeys():

					match = compiledRegex[regex].match(line)
					if match:
						if matchedcounts[regex]==0:
							gap = ""
						else:
							gap = " "

						matched[regex] = matched[regex]+gap+match.group(1).strip()
						matchedcounts[regex]+=1
			
			return (matched)
	
class Render():
	def __init__(self, **kwargs):
		self.files=kwargs.get('files')
		self.fof = kwargs.get('fof')
		self.ribbon=kwargs.get('ribbon')
		self.line = kwargs.get('line')
		self.stick=kwargs.get('stick')
		self.surface=kwargs.get('surface')
		self.cbc=kwargs.get('cbc')
		self.ray = kwargs.get('ray')
		self.view = kwargs.get('view')
		self.PyMOLPath = kwargs.get('pymolpath','/Applications/MacPyMOL.app/Contents/MacOS/MacPyMol')
		
	def go(self):
	#
	# PyMol executable path
	# In cases of MacPymol, It assumed that pymol is installed in 'Application' folder.
	# Other cases, plase set up the exact location of PyMOL executable
	#
		if self.files:
			pdbList=self.files
			
		elif self.fof and os.path.exists(self.fof):
			f = open (self.fof)
			pdbList=f.readlines()
			f.close()
			pdbList = [x.strip() for x in pdbList]
		else:
			pdbList=glob.glob('*.pdb')


		if len(pdbList)>0 :
			pymolLoadingFile = open('render.pml','w')
			pymolLoadingFile.write('bg_color white\n')
			description = {}

			style = "cartoon" #Defaults
			if self.ribbon :
				style = "ribbon"
			if self.line :
				style = "line"
			if self.stick :
				style = "stick"

			if self.view:
				if os.path.exists(self.view):
					f = open(self.view)
					view = f.readlines()
					if re.match("^set_view",view[0]):
						for line in view:
							pymolLoadingFile.write(line)
					f.close()

			for pdb in pdbList:	
				if os.path.exists(pdb):			
					print pdb
					pymolLoadingFile.write("load {0}/{1}\n".format(os.getcwd(),pdb))
					pymolLoadingFile.write("hide all\n")
					pymolLoadingFile.write("show {0}, {1}\n".format(style, pdb[:-4]))
					if self.cbc:
						pymolLoadingFile.write("util.cbc\n")
					else:
						pymolLoadingFile.write("spectrum count,rainbow,{0}\n".format(pdb[:-4]))
					if self.surface:
						pymolLoadingFile.write("show surface,{0}\n".format(pdb[:-4]))
					if not self.view:	
						pymolLoadingFile.write("orient {0}\n".format(pdb[:-4]))
					if self.ray:
						pymolLoadingFile.write("ray\n")
					pymolLoadingFile.write("png {0}.png\n".format(pdb[:-4]))
					pdbparse = PDBdescription()
					pdbparse.parse(pdb)
					description[pdb]=pdbparse
			pymolLoadingFile.close()
			print "Running PyMol..."
			p = subprocess.Popen([self.PyMOLPath,'-c','render.pml'],
								  stdout=subprocess.PIPE)
			p_stdout = p.stdout.read()
			htmlout(description)
		else:
			print "Usage : render.py"
				


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--ribbon', action='store_true', dest='ribbon',default=False,
	                    help='Draw as ribbon')
	filegroup = parser.add_mutually_exclusive_group()
	filegroup.add_argument('-T', '--files-from',  dest='fof',
	                    help='Files from file of names')
	filegroup.add_argument('-f', '--files', nargs='*', dest='files',default=[],
	                    help='File to process')
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
	parser.add_argument('-v', '--fixed_view', action='store', dest='view')
	parser.add_argument('-p', '--pymolpath', action='store', dest='pymolpath', default='/Applications/MacPyMOL.app/Contents/MacOS/MacPyMol')
	
	results = parser.parse_args()
	rendered = Render(files=results.files, fof=results.fof, ribbon=results.ribbon, line=results.line,\
						stick=results.stick, surface=results.surface, cbc=results.cbc, ray=results.ray,\
						view=results.view)	 	
	fileList = rendered.go()
	
