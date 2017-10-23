#include <iostream>
#include <string>
#include <cerrno>
#include <vector>

#include "util/err.h"

#include "boost/program_options.hpp"

#include "phase/PhaseTable.h"
#include "phase/ReadPhaser.h"

// Namespaces
using namespace phase;
using namespace std;

namespace po = boost::program_options;

// Declarations
char *progName;


/**
 * Partition reads with phased SNVs.
 *
 * @param argc Number of arguments.
 * @param argv Array of argument strings.
 *
 * @return Integer return status. 0 for no error or a constant in ``cerrno'' on failure.
 */
int main(int argc, char *argv[]) {

	//
	// Declare options
	//

	vector<string> alignFileList;
	string snvTableFileName;
	string outFileName;
	bool verbose;

	progName = argv[0];


	//
	// Command-line options
	//

	// Set default options
	outFileName = "";

	// Setup options
	po::options_description prog_opts("Phase PacBio reads in a BAM based with SNVs from short-read data");

	prog_opts.add_options()
			("help,h", "Print help")
			("snvtable,s", po::value<string>(&snvTableFileName), "SNV Table")
			("verbose,v", po::bool_switch(&verbose)->default_value(false), "Print verbose information")
			("out,o", po::value<string>(&outFileName), "Table of phased reads and HET-SNV counts")
			;

	// Positional argument: infile
	po::positional_options_description p;
	p.add("infile", -1);

	po::options_description hidden_opts("Hidden options");
	hidden_opts.add_options()
			("infile", po::value<vector<string>>(&alignFileList), "Input SAM/BAM")
			;

	// Merge options
	po::options_description cmdline_opts;
	cmdline_opts.add(prog_opts).add(hidden_opts);

	po::options_description visible("Allowed options");
	visible.add(prog_opts);

	// Parse options
	po::variables_map vm;
	try {
		po::store(po::command_line_parser(argc, argv).options(cmdline_opts).positional(p).run(), vm);

		po::notify(vm);

	} catch (const std::exception &ex) {
		err(ex);
		return EINVAL;
	}

	// Print help
	if (vm.count("help")) {
		cout << prog_opts << endl;
		return 0;
	}

	// Check arguments
	if (! vm.count("infile")) {
		err("Missing input file (see -h)");
		return EINVAL;
	}

	if (! vm.count("snvtable")) {
		err("Missing input SNV table (--snvtable, see -h)");
		return EINVAL;
	}

	if (outFileName == "") {
		err("Missing output file name (--out, see -h)");
		return EINVAL;
	}

	// Report if verbose
	if (verbose) {
		cout << "SNV Table: " << vm["snvtable"].as<string>() << endl;

		cout << "Input SAM/BAM: ";

		for (vector<string>::const_iterator iter = alignFileList.begin(); iter != alignFileList.end(); ++iter) {

			if (iter != alignFileList.begin())
				cout << ", ";

			cout << *iter;
		}

		cout << endl;
	}


	//
	// Read SNV table
	//

	PhaseTable phaseTable;

	if (verbose)
		cout << "Loading SNV table..." << endl;

	try {
		phaseTable.load(snvTableFileName);

	} catch (exception &ex) {
		err(ex);
		return EINVAL;
	}

	if (verbose)
		cout << "Loaded " << phaseTable.size() << " SNV entries." << endl;


	//
	// Phase reads
	//

	ReadPhaser readPhaser(&phaseTable);
	readPhaser.setVerbose(verbose);
	readPhaser.phase(&alignFileList, outFileName);

	return 0;
}

