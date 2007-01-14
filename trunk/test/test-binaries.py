#!/usr/bin/python
import os, re

"""Test framework to test executable programs that are part of the evoli project.

Run as "python ./test/test-binaries.py" from the top directory. The script creates a
temporary directory in which it executes the binary with a given set of parameters,
and compares the output to prerecorded output. If no problem is found, the
temporary directory is deleted. Otherwise, the framework reports an error, and the
temporary directory is left in place, so that the error can be investigated.

To create a new test, follow the steps below:
	1. Add a suitable control file to the directory "test/bin_tests". (Control files have
	   names of the form "<testname>.control", and contain the name of the binary
	   to be tested and the argument list. See tr-driver_tr for an example.) Let's
	   assume your new control file is called "newtest.control".
	2. Run the test script with "python ./test/test-binaries.py". The new test will fail.
	3. Copy the entire contents of the temporary directory into a directory "newtest"
	   in "test/bin_tests":
		cp -r newtest.DKhje8bS3.tmp_dir/ test/bin_tests/newtest
	4. That's it.
"""


testsdir = './test/bin_tests/' # directory that holds all the test data


def printTestFailed():
	print '\n   **********   Test failed!   **********'

def printTestPassed():
	print '      test passed.'

def runProgram( testdir_name, binary, arguments ):
	print '   creating temporary directory:', testdir_name
	try:
		os.mkdir( testdir_name )
	except:
		print 'cannot create directory:', testdir_name
		return False
	
	os.chdir( testdir_name )
	try:
		command = os.path.join( "../", binary ) + " " + arguments + " > testrun_log.txt"
		print '   executing:', command
		os.system( command )
	except:
		print 'cannot run executable'
		return False
	os.chdir( ".." )
	return True

def cleanup( testdir_name ):
	if not os.path.exists( testdir_name ):
		return
	for root, dirs, files in os.walk(testdir_name, topdown=False):
		for name in files:
			#print "   removing file:", os.path.join(root, name)
			os.remove( os.path.join(root, name) )
		for name in dirs:
			#print "   removing directory:", os.path.join(root, name)
			os.rmdir( os.path.join(root, name) )
	print "   removing directory:", testdir_name
	os.rmdir( testdir_name )

def parseControlFile( filename ):
	regbinary = re.compile( r'^\s*Binary:\s*(.*)' )
	regarguments = re.compile( r'^\s*Arguments:\s*(.*)' )
	
	lines = open( filename ).readlines()	
	for line in lines:
		m = regbinary.match( line.strip() )
		if m:
			binary = m.groups()[0]
			continue
		m = regarguments.match( line.strip() )
		if m:
			arguments = m.groups()[0]
			continue

	return ( binary, arguments )

def compareFiles( template_file, test_file ):
	template_raw_lines = open( template_file ).readlines()
	template_stripped_lines = []
	for line in template_raw_lines:
		sline = re.sub("\s+", " ", line)
		if sline == ' ': # skip empty lines
			continue
		template_stripped_lines.append( sline )
	test_raw_lines = open( test_file ).readlines()
	test_stripped_lines = []
	for line in test_raw_lines:
		sline = re.sub("\s+", " ", line)
		if sline == ' ': # skip empty lines
			continue
		test_stripped_lines.append( sline )
	if not len( template_stripped_lines ) == len( test_stripped_lines ):
		printTestFailed()
		print "   files"
		print "     ", template_file
		print "   and"
		print "     ", test_file
		print "   differ in length."
		return False
	a = zip( template_stripped_lines, test_stripped_lines )
	for i in a:
		if not i[0] == i[1]:
			printTestFailed()
			print "   files"
			print "     ", template_file
			print "   and"
			print "     ", test_file
			print "   differ in content."
			return False
	return True


def compareDirectories( template_dir, test_dir ):
	"""Compare the file contents of two directories. Does not work for subdirectories!"""
	if not os.path.exists( template_dir ):
		print "   ***********************************************"
		print "   Cannot complete test. No comparison data exists."
		print "   Leaving test data in temporary directory."
		print "   ***********************************************"
		return False
	template_dir_files = os.listdir( template_dir )
	template_dir_files.sort()
	test_dir_files = os.listdir( test_dir )
	test_dir_files.sort()
	result = True
	if not template_dir_files == test_dir_files:
		printTestFailed()
		print '   directories "' + template_dir + '" and "' + test_dir + '" do not contain the same set of files.'
		print '   files in "' + template_dir + '":'
		for i in template_dir_files:
			print '\t' + i
		print '   files in "' + test_dir + '":'
		for i in test_dir_files:
			print '\t' + i
		result = False
	if result:
		for file in template_dir_files:
			if os.path.isdir( os.path.join( template_dir, file ) ):
				print 'Warning: subdirectory found in directory containing test data.'
				print 'Currently the function "compareDirectories()" does not work in this case.'
				print 'Test could not be completed.'
				result = False
			if not compareFiles( os.path.join( template_dir, file ), os.path.join( test_dir, file ) ):
				result = False
	if result:
		cleanup( test_dir )
	if result:
		printTestPassed()
	return result

def testProgram( testsdir, control_file ):
	m = re.compile( r'^(.*).control$' ).match( control_file )
	testname = m.groups()[0]
	( binary, arguments ) = parseControlFile( os.path.join( testsdir, control_file ) )
	print 'running test:', testname
	testdir_name = testname + '.DKhje8bS3.tmp_dir'
	cleanup( testdir_name ) # cleanup from previous failed run, just in case
	runProgram( testdir_name, binary, arguments )
	compareDirectories( os.path.join( testsdir, testname ), testdir_name )

def runTests():
	files = os.listdir( testsdir )
	for file in files:
		if file.endswith( '.control' ):
			testProgram( testsdir, file )


runTests()