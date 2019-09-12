#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "tcrdist.hh"



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
Real
compute_NNdistance_nosorting(
	Size const num,
	Reals const & dists // should already be sorted in increasing order
)
{
	runtime_assert( dists.size() >= num );

	Real nbrdist(0), totalwt(0);
	for ( Size i=0; i<num; ++i ) {
		// Real const wt( 1.0 );
		Real const wt( 1.0 - Real(i)/num );
		totalwt += wt;
		nbrdist += wt * dists[i];
	}
	return nbrdist / totalwt;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// uses the closest decile
//
inline
Real
compute_NNdistance_nosorting(
	Reals const & dists // should already be sorted in increasing order
)
{
	if ( dists.empty() ) return 0.0;
	//return compute_NNdistance_nosorting( max( Size(1), dists.size()/2 ), dists );
	return compute_NNdistance_nosorting( max( Size(1), dists.size()/10 ), dists );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// rank_score = 0 if x < min(vals)-epsilon
// rank_score = 1 if x > max(vals)+epsilon
// otherwise interpolate
//
//
Real
compute_rank_score(
	Real const x,
	Reals const & vals,
	Real const epsilon = 1e-6
)
{
	if ( vals.empty() ) return 0.0;

	Size num_below(0), num_above(0);

	for ( Real v: vals ) {
		if ( v < x-epsilon ) {
			++num_below;
		} else if ( v > x+epsilon ) {
			++num_above;
		}
	}

	return 0.5 * ( Real(num_below)/vals.size() + 1.0 - Real(num_above)/vals.size() );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this is templated because there are different flavors of TCR depending on whether we know the V gene at the family
// or allele level
template< typename T >
void
calc_repdist_scores(
	strings const & f1_tcrs,
	vector< T > const & f1_dtcrs,
	strings const & f2_tcrs,
	vector< T > const & f2_dtcrs,
	TCRdistCalculator const & tcrdist,
	string const & outprefix
)
{

	// make a big distance matrix including all the tcrs
	Size const num_f1_tcrs( f1_tcrs.size() ), num_f2_tcrs( f2_tcrs.size() ), num_tcrs( num_f1_tcrs + num_f2_tcrs );
	runtime_assert( num_f1_tcrs>0 ); // we may not have f2 tcrs, though
	vector<Reals> D( num_tcrs);

	for ( Reals & ds : D ) { ds.resize( num_tcrs ); }

	vector< T > all_dtcrs( f1_dtcrs );
	for ( T const & dtcr : f2_dtcrs ) all_dtcrs.push_back( dtcr );
	strings all_tcrs( f1_tcrs );
	for ( string const & tcr : f2_tcrs ) all_tcrs.push_back(tcr);

	// fill distance matrix
	cerr << "computing distances: num_tcrs= " << num_tcrs << endl;
	for ( Size i=0; i< num_tcrs; ++i ) {
		D[i][i] = 0.0;
		for ( Size j=i+1; j< num_tcrs; ++j ) {
			D[i][j] = D[j][i] = tcrdist( all_dtcrs[i], all_dtcrs[j] );
		}
	}

	// get average distances for simple repertoire comparison score
	Sizes counts{0,0,0};
	Reals avgdists{0,0,0};
	for ( Size i=0; i< num_tcrs; ++i ) {
		Size const i_ind( i>=num_f1_tcrs );
		for ( Size j=0; j< num_tcrs; ++j ) {
			Size const j_ind( j>=num_f1_tcrs ), counts_ind( i_ind==j_ind ? i_ind : 2 );
			++counts[ counts_ind ];
			avgdists[ counts_ind ] += D[i][j];
		}
	}
	// take the average TCRdist between the two repertoires, subtract off the average intra-repertoire distance,
	//  so that the distance between a repertoire and itself will be 0
	Real const avgdist_score( avgdists[2]/counts[2] - 0.5*( avgdists[0]/counts[0] + avgdists[1]/counts[1] ) );


	// compute NNdists
	cerr << "computing nndists: num_tcrs= " << num_tcrs << endl;
	vector< Reals > all_nndists;
	for ( Size ii=0; ii<2; ++ii ) { // loop over the two files
		// create a tsv file with the NNdist and TStar scores
		string const outfile( outprefix + "file" + to_string(ii+1) + "_NNdist_and_TStar_scores.tsv" );
		cerr << "Writing NNdist and TStar scores to file: " << outfile << endl;
		ofstream out( outfile.c_str() );
		// header
		out << "index\ttcr\tNNdist_source_rep\tNNdist_other_rep\tTStar\n";

		Size const ib( ii == 0 ? 0 : num_f1_tcrs ), ie( ii == 0 ? num_f1_tcrs : num_tcrs );
		Reals i_nndists, dists;
		i_nndists.reserve(num_tcrs);
		dists.reserve(num_tcrs);
		for ( Size i=ib; i< ie; ++i ) {
			i_nndists.clear();
			for ( Size jj=0; jj<2; ++jj ) { // which file
				dists.clear();
				Size const jb( jj == 0 ? 0 : num_f1_tcrs ), je( jj == 0 ? num_f1_tcrs : num_tcrs );
				for ( Size j=jb; j< je; ++j ) {
					if ( i!=j ) dists.push_back( D[i][j] );
				}
				sort( dists.begin(), dists.end() ); // SLOW STEP -- could use a heap instead to get the bottom decile...
				i_nndists.push_back( compute_NNdistance_nosorting( dists ) );
			}
			Real const NNdist_src( i_nndists[ii==1] ), NNdist_other( i_nndists[ii==0] ), TStar( NNdist_other-NNdist_src );
			out << i-ib << '\t' <<
				(ii==0 ? f1_tcrs[i-ib] : f2_tcrs[i-ib] ) << '\t' <<
				NNdist_src << '\t' <<
				NNdist_other << '\t' <<
				TStar << '\n';
			all_nndists.push_back( i_nndists );
		}
		out.close();
	}


	/// now convert nndists to rank scores, get avg mutual rank score
	vector< Reals > file_self_nndists(2);
	for ( Size i=0; i<num_tcrs; ++i ) {
		Size const ind( i >= num_f1_tcrs);
		file_self_nndists[ind].push_back( all_nndists[i][ind] );
	}


	cerr << "computing NNrank scores: num_tcrs= " << num_tcrs << endl;
	Reals avg_other_nnranks{0.,0.};

	for ( Size i=0; i<num_tcrs; ++i ) {
		Size const ind( i >= num_f1_tcrs), other_ind( i < num_f1_tcrs );
		Real const other_nnrank( compute_rank_score( all_nndists[i][ other_ind ], file_self_nndists[other_ind] ) );
		avg_other_nnranks[ ind ] += other_nnrank;
	}
	avg_other_nnranks[0] /= num_f1_tcrs;
	avg_other_nnranks[1] /= num_f2_tcrs;

	// compute Repdist /////////////////////////
	//
	Real const raw_repdist_score( 0.5 * ( avg_other_nnranks[0] + avg_other_nnranks[1] ) - 0.5 );

	//  fabs here handles rare case where NNrank scores are less than 0.5, which can come up especially with
	//  small and fairly similar repertoires
	Real repdist_score( 0.5 * ( fabs( avg_other_nnranks[0] - 0.5 ) + fabs( avg_other_nnranks[1] - 0.5 ) ) );

	// self-repdists are a little wonky since i-i distances are not included within repertoires but are counted
	// between repertoires. So hack this here:
	if ( fabs(avgdist_score) < 1e-6 ) repdist_score=0.;

	string const outfile(outprefix + "repdist_scores.tsv" );
	cerr << "writing Repdist scores to file: " << outfile << endl;

	ofstream out( outfile.c_str() );
	out << "score_name\tvalue\n";
	out << "num_file1_tcrs\t" << num_f1_tcrs << '\n';
	out << "num_file2_tcrs\t" << num_f2_tcrs << '\n';
	out << "repdist\t" << repdist_score << '\n';
	out << "avgdist\t" << avgdist_score << '\n';
	out << "raw_repdist\t" << raw_repdist_score << '\n';

	// compute repertoire diversity measures
	cerr << "Computing TCRdiv repertoire diversity measures" << endl;

	Reals const thresholds{ 19.26, 38.52, 57.78 }; // 19.26 is the tcrdist value for single-chain analyses
	vector< Reals > file_tcrdivs(3);
	for ( Size ind=0; ind<2; ++ind ) { // which file
		Size const ib( ind == 0 ? 0 : num_f1_tcrs ), ie( ind == 0 ? num_f1_tcrs : num_tcrs );
		for ( Real sdev : thresholds ) {
			Real overlap_sum(0), total(0), norm;
			for ( Size i=ib; i< ie; ++i ) {
				for ( Size j=i+1; j< ie; ++j ) {
					norm = D[i][j]/sdev;
					overlap_sum += exp( -1* norm*norm );
					total += 1.;
				}
			}
			file_tcrdivs[ind].push_back( total/overlap_sum );
		}
	}
	for ( Size ind=0; ind<2; ++ind ) {
		for ( Size i=0; i<thresholds.size(); ++i ) {
			out << "file" << ind+1 <<"_TCRdiv_using_sdev" << i+1 << '\t' << file_tcrdivs[ind][i] << '\n';
		}
	}
	out.close();

}



int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "Calculate Repdist and TStar scores for two repertoires. "
			"The TCRs can be defined at the V-beta family level (e.g. V19,CASSIRSSYEQYF) "
			"or the V-allele level (e.g. TRBV19*01,CASSIRSSYEQYF or TRAV3*01,CAVPPDSWGKLQF). "
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

		cmd.parse( argc, argv );

		string const outprefix( outprefix_arg.getValue() );
		string const organism( organism_arg.getValue() );
		if ( organism != "mouse" && organism != "human" ) {
			cout << "--organism must be either mouse or human" << endl;
			exit(1);
		}

		string db( database_arg.getValue() );
		if ( db.back() != '/' ) db += "/";
		string const tcrdist_dbfile( db+"tcrdist_info_both_chains_"+organism+".txt" );

		string const tcrs_file1( tcrs_file1_arg.getValue() );
		string const tcrs_file2( tcrs_file2_arg.getValue() );

		strings f1_tcrs, f2_tcrs;

		read_tcrs_from_file( tcrs_file1, f1_tcrs );
		if ( !tcrs_file2.empty() ) read_tcrs_from_file( tcrs_file2, f2_tcrs );

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
			}
		}
		for ( int i=f2_tcrs.size()-1; i>=0; --i ) {
			if ( !tcrdist.check_tcr_string_ok( f2_tcrs[i] ) ) {
				cerr << "[WARNING] bad file2 tcr: " << i << ' ' << f2_tcrs[i] << endl;
				f2_tcrs.erase( f2_tcrs.begin()+i );
			}
		}

		if ( by_family ) {
			vector< DistanceTCR_f > f1_dtcrs, f2_dtcrs;
			foreach_( string tcr, f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			foreach_( string tcr, f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			calc_repdist_scores( f1_tcrs, f1_dtcrs, f2_tcrs, f2_dtcrs, tcrdist, outprefix );
		} else {

			vector< DistanceTCR_g > f1_dtcrs, f2_dtcrs;
			foreach_( string tcr, f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}
			foreach_( string tcr, f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}
			calc_repdist_scores( f1_tcrs, f1_dtcrs, f2_tcrs, f2_dtcrs, tcrdist, outprefix );
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
