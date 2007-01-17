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
	   names of the form "<testname>.control", and contain information about the
	   test to run. See details below, or the existing tests for a examples.) Let's
	   assume your new control file is called "newtest.control".
	2. Run the test script with "python ./test/test-binaries.py". The new test will fail.
	3. Copy the entire contents of the temporary directory into a directory "newtest"
	   in "test/bin_tests":
		cp -r newtest.DKhje8bS3.tmp_dir/ test/bin_tests/newtest
	4. That's it.

Structure of the control file:
The control file consists of lines of the form <variable name>: <value>. Comments
are *not* allowed in the control file.

The following variables exist:
Binary: The name of the program to execute.
Arguments: The set of arguments to give to the binary
Exclude: A regular expression to exclude lines from the file comparison that determines
	whether the test passed or not. The regular expression has to be enclosed in 
	slashes (/).
	For example, the directive:
		Exclude: /excluded line/
	will exclude all lines containing the string "excluded line" from the file comparison.
	There can be several Exclude directives in a control file. All will be applied.
	
You can also add a file called 'DONT_RUN' into the directory holding all the tests (i.e.,
"test/bin_tests"). All tests whose name is listed in this file will not be run. Several
tests can be listed on a single line, separated by commas, or on separate lines. You can
use the symbol "#" to comment out parts of the file. The 'DONT_RUN' file is useful
when your are debugging specific tests.
"""


testsdir = './test/bin_tests/' # directory that holds all the test data
dontrun = 'DONT_RUN' # name of the don't run file.


class ControlRecord:
        def __init__( self ):
                self.binary = ""
                self.argument = ""
                self.excludes = []


def printTestFailed():
	print '\n   **********   Test failed!   **********'

def printTestPassed():
	print '      test passed.'

def readDontRunFile():
	result = []
	filename = os.path.join( testsdir, dontrun )
	try:
		lines = file( filename ).readlines()
	except:
		return result
	for line in lines:
		# get rid of comments:
		contents = line.strip().split( '#', 1 )[0]
		if re.match( r'^\s*$', contents ): # ignore empty lines
			continue
		names = contents.split( ',' )
		for name in names:
			result.append( name.strip() )
	if not result == []:
		print "Skipping tests:"
		for name in result:
			print "\t", name
	return result


def runProgram( testdir_name, cr ):
	'''
	'testdir_name' is the name of the temporary directory in which the test is run.
	'cr' is a ControlRecord object holding all the information necessary to run the test.
	'''
	print '   creating temporary directory:', testdir_name
	try:
		os.mkdir( testdir_name )
	except:
		print 'cannot create directory:', testdir_name
		return False
	
	os.chdir( testdir_name )
	try:
		command = os.path.join( "../", cr.binary ) + " " + cr.arguments + " > testrun_log.txt"
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
	regexclude = re.compile( r'^\s*Exclude:\s*/(.*)/$' )
	
	lines = open( filename ).readlines()
	record = ControlRecord()
	for line in lines:
		m = regbinary.match( line.strip() )
		if m:
			record.binary = m.groups()[0]
			continue
		m = regarguments.match( line.strip() )
		if m:
			record.arguments = m.groups()[0]
			continue
		m = regexclude.match( line.strip() )
		if m:
			record.excludes.append( m.groups()[0] )
	return record

def prepareFileContents( lines, cr ):
	'''
	This function removes whitespace and other information from the files that are
	being compared. It is given an array of lines (i.e., strings) to process, and returns
	a new array of lines.
	
	'lines' is the array of strings to process
	'cr' is a ControlRecord object containing information about the test
	'''
	result = []
	for line in lines:
		exclude = False
		for s in cr.excludes: # test whether we should exclude the line
			if re.search( s, line ): # does line match an exclude?
				exclude = True # exclude this line
				break
		if exclude:
			continue
		sline = re.sub("\s+", " ", line) # if we are still in the game, strip whitespace
		if sline == ' ': # skip empty lines
			continue
		result.append( sline )
	return result


def compareFiles( template_file, test_file, cr ):
	'''
	'template_file' is the name of the file against which we are comparing
	'test_file' is the name of the file we compare against the template
	'cr' is a ControlRecord object holding information about the test
	'''
	template_lines = prepareFileContents( open( template_file ).readlines(), cr )
	test_lines = prepareFileContents( open( test_file ).readlines(), cr )
	if not len( template_lines ) == len( test_lines ):
		printTestFailed()
		print "   files"
		print "     ", template_file
		print "   and"
		print "     ", test_file
		print "   differ in length."
		return False
	a = zip( template_lines, test_lines )
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


def compareDirectories( template_dir, test_dir, cr ):
	"""Compare the file contents of two directories. Does not work for subdirectories!
	'cr' is a ControlRecord object
	"""
	if not os.path.exists( template_dir ):
		print "   ***********************************************"
		print "   Cannot complete test. No comparison data exists."
		print "   Leaving test data in temporary directory."
		print "   ***********************************************"
		return False
	exclude_list = ['.svn'] # files to exclude from comparison
	template_dir_files = os.listdir( template_dir )
	for i in exclude_list:
		if i in template_dir_files:
			template_dir_files.remove( i )
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
			if not compareFiles( os.path.join( template_dir, file ), os.path.join( test_dir, file ), cr ):
				result = False
	if result:
		cleanup( test_dir )
	if result:
		printTestPassed()
	return result

def testProgram( control_file, dontrun ):
	m = re.compile( r'^(.*).control$' ).match( control_file )
	testname = m.groups()[0]
	if testname in dontrun: # exclude tests listed in dontrun
		return
	cr = parseControlFile( os.path.join( testsdir, control_file ) )
	print 'running test:', testname
	testdir_name = testname + '.DKhje8bS3.tmp_dir'
	cleanup( testdir_name ) # cleanup from previous failed run, just in case
	runProgram( testdir_name, cr )
	compareDirectories( os.path.join( testsdir, testname ), testdir_name, cr )

def runTests():
	files = os.listdir( testsdir )
	dontrun = readDontRunFile()
	for file in files:
		if file.endswith( '.control' ):
			testProgram( file, dontrun )


runTests()