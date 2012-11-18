#!/usr/bin/python
#
# Wrapper of Superpose (CCP4)
# Written by MadScientist (http://madscientist.wordpress.com)
#
#

import sys, subprocess, re, os, glob, argparse
from pdbextract import *


class Aligned():
	#
	# Storage class for the structural alignment in Superpose output file
	#
	def __init__(self, **kwargs):
		self.pdb = kwargs.get('pdb')
		self.chain = kwargs.get('chain')
		self.number = kwargs.get('number')
		self.aa = kwargs.get('aa')
		self.distance = kwargs.get('distance')
		self.aligned = kwargs.get('aligned')
		if self.chain and self.number and self.aa:
			self.id = self.chain+"_"+str(self.number)+"_"+self.aa
		else:
			self.id = None
		self.alignid = None
		self.gap = kwargs.get('gap')
		self.outfile = ""

class Superpose():

	def __init__(self, **kwargs):
		self.queryPDB=kwargs.get('queryPDB')
		self.subjectPDB = kwargs.get('subjectPDB')
		self.superposedPDB=kwargs.get('superposedPDB')
		self.outfile = ""
		self.rmsd = 0
		self.aligned = 0
		self.qscore = 0	
		self.queryAlign = []
		self.subjectAlign = []
		self.querySeq = ""
		self.subjectSeq = ""
		self.upper = 0
		self.lower = 0
	

	def fastaSplit(self,fasta,width):
		#
		# Split sequence in fasta string into predefined width of string lists
		#
		cursor = 0
		end = 0
		buffer = []
		while (cursor <len(fasta)):
			if (len(fasta)-cursor) > width:
				end+=width
			else:
				end=len(fasta)		
			buffer.append(fasta[cursor:end]+"\n")
			cursor+=width
		return(buffer)	

	def FASTAoutput(self):
		#
		# Generate structural alignment as FASTA formats
		#
		fasta = []
		if len(self.queryAlign)>0:
			
			queryHeader = (self.queryAlign[self.upper].chain+":"+
						   str(self.queryAlign[self.upper].number)+"-"+
						   self.queryAlign[self.lower].chain+":"+
						   str(self.queryAlign[self.lower].number))

			fasta.append(">{0}:{1}\n".format(self.queryPDB[:len(self.queryPDB)-4],queryHeader))
			fasta=fasta+self.fastaSplit(self.querySeq,60)
			subjectHeader=	(self.subjectAlign[self.upper].chain+":"+
							 str(self.subjectAlign[self.upper].number)+"-"+
							 self.subjectAlign[self.lower].chain+":"+
							 str(self.subjectAlign[self.lower].number))

			fasta.append(">{0}:{1}\n".format(self.subjectPDB[:len(self.subjectPDB)-4],subjectHeader))
			fasta=fasta+self.fastaSplit(self.subjectSeq,60)
		return(fasta)
		

	def generateAlignment(self):
		
		oneLetter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
					'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
					'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
					'GLY':'G', 'PRO':'P', 'CYS':'C', 'MSE':'M', '-':'-'}
		sequence = ""
		aligned = ""

		for i in range(0,len(self.queryAlign)-1):
			residue = self.queryAlign[i]
		
			if residue.aa in oneLetter:
				sequence = sequence+oneLetter[residue.aa]
			else:
				sequence = sequence+"X"
			if residue.alignid.aa in oneLetter:
				aligned = aligned+oneLetter[residue.alignid.aa]
			else:
				aligned = aligned+"X"

		uTrim = 0
		lTrim = len(sequence)-1
		for i in range(len(sequence)-1):
			if self.queryAlign[i].distance:
				uTrim = i
				break

		for i in reversed(range(len(sequence)-1)):
			if self.queryAlign[i].distance:	 
				lTrim = i
				break
		self.querySeq=sequence[uTrim:lTrim+1]
		self.subjectSeq=aligned[uTrim:lTrim+1]
		self.upper = uTrim
		self.lower = lTrim

		return 

	def extractSuperposedPDB(self,filename):

		if os.path.exists(self.superposedPDB) and len(self.queryAlign)>0:
			startchain = self.queryAlign[self.upper].chain;
			start = self.queryAlign[self.upper].number;
			endchain = self.queryAlign[self.lower].chain;
			end = self.queryAlign[self.lower].number;
			[header,chains] = PDBParse(self.superposedPDB,None,None)

			if startchain ==endchain:
				extractRegion = "{0}:{1}-{2}".format(startchain,start,end)
			else:
				extractRegion = "{0}:{1}-{2} {3}:{4}-{5}".format(startchain,start,9999,endchain,1,end)
			pdbExtract = PDBExtract(extractregion = extractRegion)
			newFileName=pdbExtract.extractRegions(header,chains,filename)
			return(newFileName)
		else:
			print "problem"
			return(filename)
	


	def run(self):

		#
		# Run superpose with two PDB file (self.queryPDB and self.subjectPDB)
		#
		if os.path.exists(self.queryPDB) and os.path.exists(self.subjectPDB):
			p = subprocess.Popen(['superpose',self.queryPDB,self.subjectPDB,self.superposedPDB],
											  stdout=subprocess.PIPE)
		#
		# stdout from superpose will be saved into p_stdout
		#
			p_stdout = p.stdout.read()
			self.outfile = self.superposedPDB[:len(self.superposedPDB)-4]+'.out'
			fw = open(self.outfile,'w')
			fw.write(p_stdout)
			fw.close()
		#
		# Parsing output file. RMSD, Q-score, # of aligned residues are parsed using regex
		#
			regex = re.compile(" at RMSD =\s+([\d.]+),\s+Q=\s+([\d.]+)\s+and alignment length\s+([\d.]+)")
			match = regex.search(p_stdout)
			
			if match:
				self.rmsd= float(match.group(1))
				self.qscore=float(match.group(2))
				self.aligned=aligned=float(match.group(3))
				success=True
			else:
				success=False
				return(success)
		#
		# Parsing sequence alignment 
		#			
			data = p_stdout.split("\n")	
			capture = False

			for line in data:
				if line == "|-------------+------------+-------------|":
					capture = True
			
				if line == "`-------------\'------------\'-------------\'":
					capture = False

				if capture and line !="|-------------+------------+-------------|":

					querySecondary = line[1:2].strip()
					queryChain = line[4:5].strip()
					queryAmino = line[6:9].strip()
					queryAminoNum = line[9:13].strip()
					queryDistance = line[19:23].strip()
					if line[32:33]=="|":
						a = 5
					else:
						a = 0
					subjectSecondary = line[28+a:29+a].strip()
					subjectChain=line[31+a:32+a].strip()
					subjectAmino=line[33+a:36+a].strip()
					subjectAminoNum=line[36+a:40+a].strip()
					
					query = False
					subject = False
					if len(querySecondary)+len(queryChain)+len(queryAmino)+len(queryAminoNum)>0:
						queryresidue = Aligned(pdb=self.queryPDB, chain=queryChain, aa=queryAmino,
											   number=int(queryAminoNum),gap=False)
						query = True
					else:
						queryresidue = Aligned(pdb=self.queryPDB,aa="-",gap=True)

					if len(subjectSecondary)+len(subjectChain)+len(subjectAmino)+len(subjectAminoNum)>0:
						subjectresidue = Aligned(pdb=self.subjectPDB, chain=subjectChain, 
												 aa=subjectAmino,number=int(subjectAminoNum),gap=False)
						subject = True
					else:
						subjectresidue = Aligned(pdb=self.subjectPDB,aa="-", gap=True)
					
					if len(queryDistance):
						queryresidue.distance = float(queryDistance)
						subjectresidue.distance = float(queryDistance)

					queryresidue.alignid = subjectresidue
					subjectresidue.alignid = queryresidue
					self.queryAlign.append(queryresidue)
					self.subjectAlign.append(subjectresidue)

			self.generateAlignment()
			return(success)				

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



def htmlout(table,argument,pdb):
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
		htmloutput.write(imgtag.format(pdb[:-4]+".png", subject,subjectDescription,subject,desc,
			             rmsd[pdb],qscore[pdb],align[pdb]))

	
	for pdb,desc in tableSubjectDescriptions.items():
		match = re.match("(\S+)_(\S+).pdb",pdb)
		if match:
			query = match.group(1)
			subject = match.group(2)
		if desc == 'Alignment Failed':
			table = table+tableContent.format(query,queryDescription,subject,desc,'None','None','None')
		else:
			table = table+tableContent.format(query,queryDescription,subject,desc,
											  rmsd[pdb],qscore[pdb],align[pdb])
	htmloutput.write(table+"</table>\n")
	htmloutput.write(htmlfooter)
	htmloutput.close()



def pymolRendering(table,argument,pdb):

	PyMOLPath = '/Applications/MacPyMOL.app/Contents/MacOS/MacPyMol'
	pymolFile = open(pdb[:len(pdb)-4]+".pml","w")
	if argument.view:
		if os.path.exists(argument.view):
			f = open(argument.view)
			view = f.readlines()
			if re.match("^set_view",view[0]):
				for line in view:
					pymolFile.write(line)
			f.close()

	pymolFile.write("bg_color white\n")
	pymolFile.write("hide all\n")
	pymolFile.write("load {0}/{1}\n".format(os.getcwd(),pdb))

	for superposed,align in table.items():
		write = True
		if (argument.score and (align.qscore>=float(argument.score))):
			print "{0} is skipped. Q Score:{1}".format(align.queryPDB, align.qscore)
			write = False	

		if (argument.rmsd and (align.rmsd>=float(argument.rmsd))):
			print "{0} is skipped. Q Score:{1}".format(align.queryPDB, align.rmsd)
			write = False
		
		if (argument.align and (align.aligned>=float(argument.align))):
			print "{0} is skipped. Q Score:{1}".format(align.queryPDB, align.aligned)
			write = False

		if write:

			if argument.stepwise:			
				pymolFile.write("load {0}/{1}\n".format(os.getcwd(),align.superposedPDB))
				pymolFile.write("hide all\n")			
				pymolFile.write("show cartoon, {0}\n".format(align.superposedPDB[:-4]))
				pymolFile.write("show cartoon, {0}\n".format(pdb[:-4]))
				pymolFile.write("util.cbc {0}\n".format(align.superposedPDB[:-4]))

				if not argument.view:
					pymolFile.write("orient {0}\n".format(pdb[:-4]))
				pymolFile.write("show cartoon, {0}\n".format(pdb[:-4]))
				pymolFile.write("png {0}.png\n".format(align.superposedPDB[:-4]))
			else:
				pymolFile.write("load {0}/{1}\n".format(os.getcwd(),singlepdb))
	
	pymolFile.write("hide all\n")
	pymolFile.write("show cartoon\n")
	pymolFile.close()
	
	if argument.render:
		print "Running PyMol..."
		p = subprocess.Popen([PyMOLPath,'-c',pdb+'.pml'],
							  stdout=subprocess.PIPE)
		p_stdout = p.stdout.read()
	# if argument.stepwise:
		# htmlout(pdb,pdbparsing(pdb, '^TITLE    [ \d](.*)$'),renderingSubjectDescriptions,
		# 					   tableSubjectDescriptions,RMSDDic,QDic,alignedDic,pdb+'.html')
	


def main(argument):

	pdbList=argument.files
	
	if len(pdbList)>0 :
		for pdb in pdbList:	
			if os.path.exists(pdb):

				subjectList = glob.glob('*.pdb')
				pdbLoadList = []
				RMSDDic = {}
				QDic = {}
				alignedDic = {}
	
				if len(subjectList)>0:
					table = {}
					for onePdb in subjectList:
						if (onePdb != pdb) :
							print "Processing {0}".format(onePdb)
							superposedFilename = pdb[:len(pdb)-4]+'_'+onePdb

							sup = Superpose(queryPDB=onePdb, subjectPDB=pdb, 
											superposedPDB=superposedFilename)
							if sup.run():
								if argument.extract:
									superposedFilename = sup.extractSuperposedPDB(superposedFilename)
									sup.superposedPDB = superposedFilename
								
								table[superposedFilename]=sup
								pdbLoadList.append(superposedFilename)
								if argument.stalign:
									aligned = sup.FASTAoutput()
									if aligned:
										f = open(superposedFilename[:len(superposedFilename)-4]+".fasta","w")
										f.writelines(aligned)
										f.close()		
							
					pymolRendering(table, argument,pdb)
				else:
					print "File is not exist!"
					sys.exit()		

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--reference_pdb', nargs='+', dest='files',default=[],
	                    help='Files to process')
	parser.add_argument('-e', '--extract_superimpose', action='store_true', dest='extract',default=False,
						help='Extract superposed region')
	parser.add_argument('-d', '--rmsd_cutoff', action='store', dest='rmsd',type=float,
						help='Set RMSD cutoff')
	parser.add_argument('-s', '--score_cutoff', action='store', dest='score',type=float,
						help='Set Q Score cutoff')	
	parser.add_argument('-a', '--align_cutoff', action='store', dest='align',type=float,
						help='Set alignment length cutoff')		
	parser.add_argument('-t', '--stepwise_compare_rendering', action='store_true', dest='stepwise', default=False,
						help='Render 1:1 comparison')
	parser.add_argument('-p', '--pymol', action='store_true', dest='render', default=False,
						help='Execute PyMOL')
	parser.add_argument('-n', '--align', action='store_true', dest='stalign', default=False,
						help='Save structure based sequence alignment')
	parser.add_argument('-v', '--fixed_view', action='store', dest='view',
						help='Filename which contains set_view command of PyMOL')

	results = parser.parse_args()
	main(results)