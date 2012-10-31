#
# pdb blast search using RESTful web service form www.rcsb.org
#
#

import urllib2, urllib, os
import xml.etree.ElementTree as ET

url = 'http://www.rcsb.org/pdb/rest/search'
#
# 

q="""
<?xml version="1.0" encoding="UTF-8"?>
	<orgPdbCompositeQuery version="1.0">
	
		<queryRefinement>
			<queryRefinementLevel>0</queryRefinementLevel>
	
				<orgPdbQuery>
					<queryType>org.pdb.query.simple.SequenceQuery</queryType>
					<sequence>MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYALPHAILRLDLAGRDLTDYLMKILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDGQVITIGNERFRCPEALFQPSFLGMESCGIHETTFNSIMKCDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWISKQEYDESGPSIVHRKCF</sequence>
					<eCutOff>1e-30</eCutOff>
					<searchTool>blast</searchTool>
				</orgPdbQuery>
	
		</queryRefinement>
	
	</orgPdbCompositeQuery>
"""


queryXML = """
<?xml version="1.0" encoding="UTF-8"?>
<orgPdbCompositeQuery version=\"1.0\">
	<queryRefinement>
		<queryRefinementLevel>0</queryRefinementLevel>
		
			<orgPdbQuery>
				<queryType>org.pdb.query.simple.SequenceQuery</queryType>
				<sequence>{0}</sequence>
				<eCutOff>{1}</eCutOff>
				<searchTool>blast</searchTool>
			</orgPdbQuery>
</queryRefinement>
<queryRefinement>
	<queryRefinementLevel>1</queryRefinementLevel>
		<orgPdbQuery>
			<queryType>org.pdb.query.simple.ResolutionQuery</queryType>
			<description>ResolutionQuery </description>
			<refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
			<refine.ls_d_res_high.max>{2}</refine.ls_d_res_high.max>
		</orgPdbQuery>
</queryRefinement>
<queryRefinement>
	<queryRefinementLevel>2</queryRefinementLevel>
		<orgPdbQuery>
		 <queryType>org.pdb.query.simple.HomologueReductionQuery</queryType>
    	 <identityCutoff>90</identityCutoff>
		</orgPdbQuery>
</queryRefinement>
</orgPdbCompositeQuery>
"""
#
#
sequence = """
ETTDSFWEVGNYKRTVKRIDDGHRLCNDLMNCVQERAKIEKAYGQQLTDWAKRWRQLIEKGPQYGSLERAWGAIMTEADK
VSELHQEVKNNLLNEDLEKVKNWQKDAYHKQIMGGFKETKEAEDGFRKAQKPWAKKMKELEAAKKAYHLACKEEKLAMTR
EMNSKTEQSVTPEQQKKLQDKVDKCKQDVQKTQEKYEKVLEDVGKTTPQYMENMEQVFEQCQQFEEKRLVFLKEVLLDIK
RHLNLAENSSYIHVYRELEQAIRGADAQEDLRWFRSTSGPGMPMNWPQFEEWNPD
"""
resolution="2.8"
eCutOff="1e-30"
biological_unit = True
reprot = """
http://www.rcsb.org/pdb/rest/customReport?pdbids={0}&customReportColumns=structureId,structureTitle,resolution,experimentalTechnique,depositionDate,structureAuthor,classification,structureMolecularWeight&service=wsdisplay&format=xml
"""

queryText = queryXML.format(sequence, eCutOff, resolution)
print "querying PDB...\n"

#
# Fetch pdb list using XML query in queryText
# and pdb list was in pdblist
#

#req = urllib2.Request(url, data=queryText)
req = urllib2.Request(url, data=q)
f = urllib2.urlopen(req)

result = f.read()

if result:	
	pdblist = []
	pdblist = result.split()


	if biological_unit:
		pdburl = "ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/{0}.pdb1.gz"
	else:
		pdburl = "http://www.rcsb.org/pdb/files/{0}.pdb.gz "
	
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
