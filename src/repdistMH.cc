#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "tcrdist.hh"




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this is templated because there are different flavors of TCR depending on whether we know the V gene at the family
// or allele level
template< typename T >
void
calc_MH_score(
	Real const sigma,
	Reals const & f1_probs,
	strings const & f1_tcrs,
	vector< T > const & f1_dtcrs,
	Reals const & f2_probs,
	strings const & f2_tcrs,
	vector< T > const & f2_dtcrs,
	TCRdistCalculator const & tcrdist,
	string const & outprefix
)
{
	Real const neginvar( -1.0/(sigma*sigma) );

	Reals all_simpsons;

	for ( Size r=1; r<= 3; ++r ) { // 1-1 2-2 1-2
		vector< T > const & i_dtcrs( r==2 ? f2_dtcrs : f1_dtcrs );
		vector< T > const & j_dtcrs( r==1 ? f1_dtcrs : f2_dtcrs );
		Reals const & i_probs( r==2 ? f2_probs : f1_probs );
		Reals const & j_probs( r==1 ? f1_probs : f2_probs );

		Size const i_end( i_dtcrs.size() ), j_end( j_dtcrs.size() );

		cerr << "compute MH: " << i_end << ' ' << j_end << endl;
		Real dist, simpson(0);
		for ( Size i=0; i< i_end; ++i ) {
			for ( Size j=0; j< j_end; ++j ) {
				dist = tcrdist( i_dtcrs[i], j_dtcrs[j] );
				simpson += i_probs[i] * j_probs[j] * exp( neginvar*dist*dist );
			}
		}
		all_simpsons.push_back( simpson );
		cerr << "Simpsons " << i_end << ' ' << j_end << ' ' << simpson << endl;
	}

	string const outfile(outprefix + "repdist_scores.tsv" );
	cerr << "writing MH scores to file: " << outfile << endl;

	Real const mh_index( 2*all_simpsons[2]/(all_simpsons[0]+all_simpsons[1]) );

	ofstream out( outfile.c_str() );
	out << "score_name\tvalue\n";
	out << "sigma\t" << sigma << '\n';
	out << "MH\t" << mh_index << '\n';
 	out << "Simpsons_file1\t" << all_simpsons[0] << '\n';
 	out << "Simpsons_file2\t" << all_simpsons[1] << '\n';
 	out << "Simpsons_file1_vs_file2\t" << all_simpsons[2] << '\n';
	out.close();

}





int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "Calculate version of Morisita-Horn (MH) overlap index that "
			"counts sharing of similar as well as identical sequences. Degree of similarity "
			"is controlled by a sigma parameter which acts like a standard deviation: the "
			"contribution for a pair of TCRs falls off exponentially as a function of "
			"exp( -1 * ( distance/sigma ) * (distance/sigma) ).  Input TCR format looks like "
			"34,V19,CASSIRSSYEQYF or 34,TRBV19*01,CASSIRSSYEQYF or 34,TRAV3*01,CAVPPDSWGKLQF "
			"where 34 is the clone size (weight) for the given TCR in the MH calculation. "
			"Output is written to .tsv files whose filenames are prepended with the string "
			"that is passed in with the -o or --outfile_prefix argument. See the test/ "
			"directory for examples of input files and usage.",
			' ', "0.1" );

		// path to database files
 		TCLAP::ValueArg<string> database_arg("d","database","Path to database directory",true,
			"unk","path",cmd);

 		TCLAP::ValueArg<string> organism_arg("g","organism","TCR source organism; either mouse or human",
			true,"unk","string",cmd);

 		TCLAP::ValueArg<string> tcrs_file2_arg("j","tcrs_file2","File containing the second set of TCRs.",
			true,"unk","string",cmd);

 		TCLAP::ValueArg<string> tcrs_file1_arg("i","tcrs_file1","File containing TCRs for distance computation. "
			"Will compute matrix of distances between TCRs in this file and the TCRs in tcrs_file2.",true,"unk","string",cmd);

 		TCLAP::ValueArg<string> outprefix_arg("o","outprefix","String that will be prepended to any output files created.",
			true,"unk","string",cmd);

 		TCLAP::ValueArg<Real> sigma_arg("s","sigma","Standard-deviation or sigma parameter that controls how quickly"
			"similarity decays with increasing TCRdist.",
		 true,12.5,"float",cmd);

		// TCLAP::SwitchArg with_counts_arg("w","with_counts","tcrs include count as first element", cmd, false);

		cmd.parse( argc, argv );

		string const outprefix( outprefix_arg.getValue() );
		string const organism( organism_arg.getValue() );
		if ( organism != "mouse" && organism != "human" ) {
			cout << "--organism must be either mouse or human" << endl;
			exit(1);
		}
		Real const sigma( sigma_arg.getValue() );

		string db( database_arg.getValue() );
		if ( db.back() != '/' ) db += "/";
		string const tcrdist_dbfile( db+"tcrdist_info_both_chains_"+organism+".txt" );

		string const tcrs_file1( tcrs_file1_arg.getValue() );
		string const tcrs_file2( tcrs_file2_arg.getValue() );
		//bool const with_counts( with_counts_arg.getValue() );

		strings f1_tcrs, f2_tcrs;
		Sizes f1_counts, f2_counts;

		read_tcrs_and_counts_from_file( tcrs_file1, f1_tcrs, f1_counts );
		read_tcrs_and_counts_from_file( tcrs_file2, f2_tcrs, f2_counts );

		if ( f1_tcrs.empty() || f2_tcrs.empty() ) {
			cout << "No TCRs read successfully from tcrs_file:  num1= " << f1_tcrs.size() << ' ' <<
				" num2= " << f2_tcrs.size() << endl;
			exit(1);
		}


		// which kind of input do we have? V-family or exact V-gene?
		bool by_family;
		char tcr_chain( 'B' ); // the default

		{ // hacky stuff here
			string const random_tcr( f1_tcrs[ f1_tcrs.size()/2 ] ); // the first could be a csv header, for example
			if ( random_tcr[0] == 'V' ) {
				by_family = true;
			} else {
				by_family = false; // full allele-level gene information
				runtime_assert( random_tcr.substr(0,2) == "TR" );
				tcr_chain = random_tcr[2];
			}
		} // scope for io checking

		TCRdistCalculator tcrdist( tcr_chain, tcrdist_dbfile );

		for ( int i=f1_tcrs.size()-1; i>=0; --i ) {
			if ( !tcrdist.check_tcr_string_ok( f1_tcrs[i] ) ) {
				cerr << "[WARNING] bad file1 tcr: " << i << ' ' << f1_tcrs[i] << endl;
				f1_tcrs.erase( f1_tcrs.begin()+i );
				f1_counts.erase( f1_counts.begin()+i );
			}
		}
		for ( int i=f2_tcrs.size()-1; i>=0; --i ) {
			if ( !tcrdist.check_tcr_string_ok( f2_tcrs[i] ) ) {
				cerr << "[WARNING] bad file2 tcr: " << i << ' ' << f2_tcrs[i] << endl;
				f2_tcrs.erase( f2_tcrs.begin()+i );
				f2_counts.erase( f2_counts.begin()+i );
			}
		}

		// convert counts to probabilities
		runtime_assert( f1_counts.size() == f1_tcrs.size() );
		runtime_assert( f2_counts.size() == f2_tcrs.size() );
		Reals f1_probs, f2_probs;
		Size total(0);
		for ( Size c: f1_counts ) total += c;
		for ( Size c: f1_counts ) f1_probs.push_back( Real(c)/total );
		total = 0;
		for ( Size c: f2_counts ) total += c;
		for ( Size c: f2_counts ) f2_probs.push_back( Real(c)/total );

		if ( by_family ) {
			vector< DistanceTCR_f > f1_dtcrs, f2_dtcrs;
			for ( string tcr : f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			for ( string tcr : f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			calc_MH_score( sigma, f1_probs, f1_tcrs, f1_dtcrs, f2_probs, f2_tcrs, f2_dtcrs, tcrdist, outprefix );
		} else {

			vector< DistanceTCR_g > f1_dtcrs, f2_dtcrs;
			for ( string tcr : f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}
			for ( string tcr : f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}

			calc_MH_score( sigma, f1_probs, f1_tcrs, f1_dtcrs, f2_probs, f2_tcrs, f2_dtcrs, tcrdist, outprefix );
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
