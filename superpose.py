#!/usr/bin/python
#
# Wrapper of Superpose (CCP4)
#
#

import sys, subprocess, re, os, glob, argparse


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



def htmlout(queryid, queryDescription, renderingSubjectDescriptions,tableSubjectDescriptions,rmsd,qscore,align,filename):
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
		table{
			    font-family: "Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
			    font-size: 12px;
			    #margin: 45px;
			    width:100%;
			    text-align: left;
			    border-collapse: collapse;  
				}
			tr:nth-child(odd) {
				background-color: #999999;
			}
	</style>
	<body>
"""
	imgtag = """
		<img src="./{0}">
		<p><a href="http://www.rcsb.org/pdb/explore/explore.do?structureId={1}">{1}:{2}</a></p>
		<p><a href="http://www.rcsb.org/pdb/explore/explore.do?structureId={3}">{3}:{4}</a></p>
		<p>RMSD: {5} Angstrom</p>
		<p>Q Score: {6}</p>
		<p>Aligned: {7}</p>
"""

	tableHeader = """
		<table>
		<tr>
			<td>PDB Reference id </td>
			<td>Reference </td>
			<td>PDB Query id </td>
			<td>Query </td>
			<td>RMSD (A) </td>
			<td>Q Score (0-1) </td>
			<td>Aligned AA (#) </td>
		</tr>
"""
	tableContent = """
		<tr>
			<td> <a href="http://www.rcsb.org/pdb/explore/explore.do?structureId={0}">{0} </td>
			<td> {1} </td>
			<td> <a href="http://www.rcsb.org/pdb/explore/explore.do?structureId={2}">{2} </td>
			<td> {3} </td>
			<td> {4} </td>
			<td> {5} </td>
			<td> {6} </td>
		</tr>
"""


	htmlfooter = """
	</body>
</html>
"""
	htmloutput = open(filename,'w')
	htmloutput.write(htmlheader)
	table = tableHeader
	for pdb,desc in renderingSubjectDescriptions.items():

		match = re.match("(\S+)_(\S+).pdb",pdb)
		if match:
			query = match.group(1)
			subject = match.group(2)
		htmloutput.write(imgtag.format(pdb[:-4]+".png", query,queryDescription,subject,desc,rmsd[pdb],qscore[pdb],align[pdb]))

	for pdb,desc in tableSubjectDescriptions.items():
		match = re.match("(\S+)_(\S+).pdb",pdb)
		if match:
			query = match.group(1)
			subject = match.group(2)
		if desc == 'Alignment Failed':
			table = table+tableContent.format(query,queryDescription,subject,desc,'None','None','None')
		else:
			table = table+tableContent.format(query,queryDescription,subject,desc,rmsd[pdb],qscore[pdb],align[pdb])
	htmloutput.write(table+"</table>\n")
	htmloutput.write(htmlfooter)
	htmloutput.close()


def main(argument):

	pdbList=argument.files
	PyMOLPath = '/Applications/MacPyMOL.app/Contents/MacOS/MacPyMol'
	if len(pdbList)>0 :
		for pdb in pdbList:	
			if os.path.exists(pdb):

				subjectList = glob.glob('*.pdb')
				pdbLoadList = []
				RMSDDic = {}
				QDic = {}
				alignedDic = {}
				
				informationFile = open(pdb+'.txt','w')
				pymolloadingFile = open(pdb+'.pml','w')
				if argument.view:
					if os.path.exists(argument.view):
						f = open(argument.view)
						view = f.readlines()
						if re.match("^set_view",view[0]):
							for line in view:
								pymolloadingFile.write(line)
						f.close()

				if len(subjectList)>0:
					tableSubjectDescriptions = {}
					for onePdb in subjectList:
						if (onePdb != pdb) :
							print "Processing {0}".format(onePdb)
							superposedFilename = pdb[:len(pdb)-4]+'_'+onePdb
							#
							# run superpose using subprocess. superpose should be in your PATH.
							#
							p = subprocess.Popen(['superpose',onePdb,pdb,superposedFilename],
												  stdout=subprocess.PIPE)
							#
							# stdout from superpose will be saved into p_stdout
							#
							p_stdout = p.stdout.read()
							
							#
							# regex for the parsing output from superpose
							# RMSD, Q, Transformation matrix and vector will be parsed out.
							
							regex = re.compile(" at RMSD =\s+([\d.]+),\s+Q=\s+([\d.]+)\s+and alignment length\s+([\d.]+)")
							regex2 = re.compile("Rx         Ry         Rz           T\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+")
							match = regex.search(p_stdout)

							if match:
								RMSD= float(match.group(1))
								Q=float(match.group(2))
								aligned=float(match.group(3))
								RMSDDic[superposedFilename] = RMSD
								QDic[superposedFilename] = Q
								alignedDic[superposedFilename] = aligned
								tableSubjectDescriptions[superposedFilename]=pdbparsing(superposedFilename, '^TITLE    [ \d](.*)$')		
							else:
								RMSD= 99
								Q=0
								aligned=0
								RMSDDic[superposedFilename] = RMSD
								QDic[superposedFilename] = Q
								alignedDic[superposedFilename] = aligned
								tableSubjectDescriptions[superposedFilename]='Alignment Failed'
								
							tablelist = {}
							match2 = regex2.search(p_stdout)
							if match2:
								matrix=[[float(match2.group(1)),float(match2.group(2)),float(match2.group(3))],
										[float(match2.group(5)),float(match2.group(6)),float(match2.group(7))],
										[float(match2.group(9)),float(match2.group(10)),float(match2.group(11))]]
								vector=[float(match2.group(4)),float(match2.group(5)),float(match2.group(6))]
							
							informationFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(pdb, onePdb, superposedFilename,RMSD,Q,aligned))
							if match:
								pdbLoadList.append(superposedFilename)
				
				informationFile.close() 
				pymolloadingFile.write("bg_color white\n")
				pymolloadingFile.write("hide all\n")
				pymolloadingFile.write("load {0}/{1}\n".format(os.getcwd(),pdb))
				renderingSubjectDescriptions={}

				for singlepdb in pdbLoadList:
					write = True
					
					if (argument.score and (QDic[singlepdb]>=float(argument.score))):
						print "{0} is skipped. Q Score:{1}".format(singlepdb, RMSDDic[singlepdb])
						write = False	

					if  (argument.rmsd and (RMSDDic[singlepdb]>=float(argument.rmsd))):
						print "{0} is skipped. RMSD:{1}".format(singlepdb, RMSDDic[singlepdb])
						write = False

					if (argument.align and (alignDic[singlepdb]>=float(argument.align))):
						print "{0} is skipped. RMSD:{1}".format(singlepdb, RMSDDic[singlepdb])
						write = False

					if write:

						if argument.stepwise:
							renderingSubjectDescriptions[singlepdb]=pdbparsing(singlepdb, '^TITLE    [ \d](.*)$')							
							pymolloadingFile.write("load {0}/{1}\n".format(os.getcwd(),singlepdb))
							pymolloadingFile.write("hide all\n")			
							pymolloadingFile.write("show cartoon, {0}\n".format(singlepdb[:-4]))
							pymolloadingFile.write("show cartoon, {0}\n".format(pdb[:-4]))
							pymolloadingFile.write("color red, {0}\n".format(singlepdb[:-4]))
							pymolloadingFile.write("color yellow, {0}\n".format(pdb[:-4]))
							if not argument.view:
								pymolloadingFile.write("orient {0}\n".format(pdb[:-4]))
							pymolloadingFile.write("show cartoon, {0}\n".format(pdb[:-4]))
							pymolloadingFile.write("png {0}.png\n".format(singlepdb[:-4]))
						else:
							pymolloadingFile.write("load {0}/{1}\n".format(os.getcwd(),singlepdb))
				
				pymolloadingFile.write("hide all\n")
				pymolloadingFile.write("show cartoon\n")
				pymolloadingFile.close()
				
				if argument.render:
					print "Running PyMol..."
					p = subprocess.Popen([PyMOLPath,'-c',pdb+'.pml'],
										  stdout=subprocess.PIPE)
					p_stdout = p.stdout.read()
				if argument.stepwise:
					htmlout(pdb,pdbparsing(pdb, '^TITLE    [ \d](.*)$'),renderingSubjectDescriptions,tableSubjectDescriptions,RMSDDic,QDic,alignedDic,pdb+'.html')
		sys.exit()		


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--reference_pdb', nargs='+', dest='files',default=[],
	                    help='Files to process')
	parser.add_argument('-d', '--rmsd_cutoff', action='store', dest='rmsd',
						help='Set RMSD cutoff')
	parser.add_argument('-s', '--score_cutoff', action='store', dest='score',
						help='Set Q Score cutoff')	
	parser.add_argument('-a', '--align_cutoff', action='store', dest='align',
						help='Set alignment length cutoff')		
	parser.add_argument('-t', '--stepwise_compare_rendering', action='store_true', dest='stepwise', default=False,
						help='Render 1:1 comparison')
	parser.add_argument('-p', '--pymol', action='store_true', dest='render', default=False,
						help='Execute PyMOL')
	parser.add_argument('-v', '--fixed_view', action='store', dest='view',
						help='Filename which contains set_view command of PyMOL')

	results = parser.parse_args()
	main(results)