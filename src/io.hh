// some input/output routines
//

#ifndef INCLUDED_io_HH
#define INCLUDED_io_HH

#include "misc.hh"


void
check_file(
	ifstream const & data,
	string const & filename
)
{
	if ( !data.good() ) {
		cerr << "\n\n[ERROR]  Unable to open " << filename << "\n\n" << endl;
		exit(1);
	}
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// just looks for "tcr: " in line and takes the next thing...
// or if not "tcr:" in the line it takes the whole line
//
void
read_tcrs_from_file(
	string const filename,
	strings & tcrs
)
{
	ifstream data( filename.c_str() );
	check_file( data, filename );

	string line;
	while (getline( data, line ) ) {
		Size pos( line.find("tcr: ") );
		if ( pos != string::npos ) {
			istringstream l( line.substr(pos+5) );
			string tcr;
			l >> tcr;
			tcrs.push_back( tcr );
		} else {
			tcrs.push_back( line );
		}
	}
	data.close();

}

#endif
