#!/usr/bin/python
#
# pdb blast search using RESTful web service form www.rcsb.org
#
#

import urllib2, urllib, os,sys,re
import xml.etree.ElementTree as ET
import argparse

def readfasta(filename) :

	sequence = ''
	sequenceName = ''
	if os.path.exists(filename):
		
		f = open(filename)
		while True:
			line = f.readline().strip()
			if not line: 
				break 
			match = re.match('^>(.*)',line)
			if match:
				sequenceName = match.group(1)
			else:
				sequence = sequence+line			
		return(sequenceName,sequence)	    	
	else:
		print "Fasta file is not exists!"
		sys.exit()
		return()


def generateQuery(parsing_results) :
#
# Generate XML query for PDB RESTful service based on input arguments
#
#

	keywordQuery= """
				<orgPdbQuery>
					<queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
					<keywords>{0}</keywords>
				</orgPdbQuery>
	"""

	descriptionQuery= """
				<orgPdbQuery>
					<queryType>org.pdb.query.simple.StructDescQuery</queryType>
					<entity.pdbx_description.comparator>{1}</entity.pdbx_description.comparator>
					<entity.pdbx_description.value>{0}</entity.pdbx_description.value>
				</orgPdbQuery>
	"""

	resolutionQuery = """
				<orgPdbQuery>
					<queryType>org.pdb.query.simple.ResolutionQuery</queryType>
					<refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
					<refine.ls_d_res_high.max>{0}</refine.ls_d_res_high.max>
				</orgPdbQuery>
	"""

	homologReductionQuery="""

				<orgPdbQuery>
					<queryType>org.pdb.query.simple.HomologueReductionQuery</queryType>
				    <identityCutoff>{0}</identityCutoff>
				</orgPdbQuery>
	"""

	pfamIDQuery ="""
				<orgPdbQuery>
				    <queryType>org.pdb.query.simple.PfamIdQuery</queryType>
				    <pfamID>{0}</pfamID>
				</orgPdbQuery>
	"""

	BLASTQuery = """
				<orgPdbQuery>
					<queryType>org.pdb.query.simple.SequenceQuery</queryType>
					<sequence>{0}</sequence>
					<eCutOff>{1}</eCutOff>
					<searchTool>blast</searchTool>
				</orgPdbQuery>
	"""

	molecularWeightQuery = """
				<orgPdbQuery>
					<queryType>org.pdb.query.simple.MolecularWeightQuery</queryType>
					<mvStructure.structureMolecularWeight.min>{0}</mvStructure.structureMolecularWeight.min>
    				<mvStructure.structureMolecularWeight.max>{1}</mvStructure.structureMolecularWeight.max>
					<searchTool>blast</searchTool>
				</orgPdbQuery>
	"""
	techniqueQuery = """
				<orgPdbQuery>
					<queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
				    <mvStructure.expMethod.value>{0}</mvStructure.expMethod.value>
				</orgPdbQuery>
	"""

	queryHead="""
<?xml version="1.0" encoding="UTF-8"?>
<orgPdbCompositeQuery version=\"1.0\">
"""

	queryEnd = """
	</orgPdbCompositeQuery>
	"""
	queryRefinementHead="""
		<queryRefinement>
			<queryRefinementLevel>{0}</queryRefinementLevel>
	"""
	queryRefinementEnd ="""
		</queryRefinement>
	"""

	query = queryHead
	refinementLevel = 0
	if results.fasta_filename:
		(sequenceName, sequence) = readfasta(results.fasta_filename)
		query = query+queryRefinementHead.format(refinementLevel)+BLASTQuery.format(sequence, results.evalue)+queryRefinementEnd
		refinementLevel+=1
	if results.keyword:
		query = query+queryRefinementHead.format(refinementLevel)+keywordQuery.format(results.keyword)+queryRefinementEnd
		refinementLevel+=1
	if results.description:
		query = query+queryRefinementHead.format(refinementLevel)+descriptionQuery.format(results.description, results.description_comparator)+queryRefinementEnd
		refinementLevel+=1
	if results.resolution:
		query = query+queryRefinementHead.format(refinementLevel)+resolutionQuery.format(results.resolution)+queryRefinementEnd
		refinementLevel+=1
	if results.homologyReduction:
		query = query+queryRefinementHead.format(refinementLevel)+homologReductionQuery.format(results.homologyReduction)+queryRefinementEnd
		refinementLevel+=1
	if results.pfamID:
		query = query+queryRefinementHead.format(refinementLevel)+pfamIDQuery.format(results.pfamID)+queryRefinementEnd
		refinementLevel+=1
	techlist = ['X-RAY','SOLUTION NMR', 'SOLID-STATE NMR', 'ELECTRON MICROSCOPY', 'ELECTRON CRYSTALLOGRAPHY', 'FIBER DIFFRACTION',
				'NEUTRON DIFFRACTION', 'SOLUTION SCATTERING', 'OTHER', 'HYBRID' ]
	if results.technique and results.technique in techlist:
		query = query+queryRefinementHead.format(refinementLevel)+techniqueQuery.format(results.technique)+queryRefinementEnd
		refinementLevel+=1
			
	if results.molecularWeight:
		if re.match("\d+-\d+",results.molecularWeight):
			(low,high) = results.molecularWeight.split('-')
			query = query+queryRefinementHead.format(refinementLevel)+molecularWeightQuery.format(low,high)+queryRefinementEnd
			refinementLevel+=1
		else:
			print "Format error!"
			sys.exit()

	query = query+queryEnd

	return(query)

def executeQuery(queryText,options):

	url = 'http://www.rcsb.org/pdb/rest/search'

	biological_unit = True
	download = True
	reprot = """
	http://www.rcsb.org/pdb/rest/customReport?pdbids={0}&customReportColumns=structureId,structureTitle,resolution,experimentalTechnique,depositionDate,structureAuthor,classification,structureMolecularWeight&service=wsdisplay&format=xml
	"""

	print "querying PDB...\n"

	#
	# Fetch pdb list using XML query in queryText
	# and pdb list was in pdblist
	#
	print queryText
	req = urllib2.Request(url, data=queryText)
	f = urllib2.urlopen(req)

	result = f.read()

	if result:
		print result	
		pdblist = []
		pdblist = result.split()
		print len (pdblist)

		if options.biological_unit:
			pdburl = "ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/{0}.pdb1.gz"
		else:
			pdburl = "http://www.rcsb.org/pdb/files/{0}.pdb.gz "
		

		if options.download:
			f=open('load.pml','w')

			for pdb in pdblist:
				print "download {0}".format(pdb)
				urllib.urlretrieve(pdburl.format(pdb.lower()), pdb+".pdb.gz")
				os.system("gunzip {0}.pdb.gz".format(pdb))
				f.write("load {0}/{1}.pdb\n".format(os.getcwd(),pdb))
				
			f.write("hide all\n".format(pdb))
			f.write("show cartoon\n".format(pdb))

			for i in range(1,len(pdblist)):
				f.write("cealign {0}, {1}\n".format(pdblist[0],pdblist[i]))
			f.close()
		
		req = urllib2.Request(reprot.format(",".join(pdblist)))
		f = urllib2.urlopen(req)
		result = f.read()

		if result:
			f=open('list.txt','w')	
			root = ET.fromstring(result)
			for record in root.iter('record'):
				recordline = []
				for field in record:
					recordline.append(field.text)
				f.write("\t".join(recordline)+"\n")
			f.close()
		else:
		    print "Failed to retrieve results" 
	else:
	    print "Failed to retrieve results" 
	sys.exit()


if __name__ == "__main__":


	parser = argparse.ArgumentParser()
	parser.add_argument('-b', action='store', dest='fasta_filename',
	                    help='BLAST search of given fasta file')
	parser.add_argument('-be', action='store', dest='evalue',default=1e-10,
	                    help='e-value cutoff of BLAST search')
	parser.add_argument('-k', action='store', dest='keyword',
	                    help='search pdb based on keywords ')
	parser.add_argument('-de', action='store', dest='description',
	                    help='search pdb based on keywords ')
	parser.add_argument('-c', action='store', dest='keyword_comparator',default='contains',
	                    help='search pdb does not contain keywords')
	parser.add_argument('-r', action='store', dest='resolution',
	                    help='resolution limit')
	parser.add_argument('-m', action='store', dest='molecularWeight',
	                    help='molecular weight range: min-max Dalton (example : "20000-40000")')

	parser.add_argument('-t', action='store', dest='technique',
	                    help='technique in X-RAY,SOLUTION NMR,SOLID-STATE NMR,ELECTRON MICROSCOPY,ELECTRON CRYSTALLOGRAPHY,FIBER DIFFRACTION,NEUTRON DIFFRACTION,SOLUTION SCATTERING,OTHER,HYBRID')
	parser.add_argument('-i', action='store', dest='homologyReduction',default=90,
	                    help='Reduce Homolog (Exclude identity)')
	parser.add_argument('-p', action='store', dest='pfamID',
	                    help='search based on Pfam ID')
	parser.add_argument('-d', action='store_true', dest='download', default=False,
	                    help='Download pdb files')
	parser.add_argument('-asu', action='store_true', dest='biological_unit', default=False,
	                    help='Download Whole content of asymetric Unit')


	results = parser.parse_args()
	query = generateQuery(results)
	executeQuery(query,results)
	sys.exit()

