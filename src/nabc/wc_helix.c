// wc_helix() - create Watson/Crick duplex

#include "nabc.h"

//   see user's manual for explanation of the code.
//   If sreslib is empty, use nab.lib libraries by default;
//     (same for aresib; sreslib and areslib should both be either
//     empty or non-empty).

MOLECULE_T *wc_helix(
	char *seq,  char *snatype, char *aseq, char *anatype,
	double xoff, double incl, double twist, double rise,
	char *opts )
{
	MOLECULE_T  *m1, *m2, *m3;
	MATRIX_T xomat, inmat, mat;
	char *arname, *srname;
	RESIDUE_T *sres, *ares;
	int	has_s, has_a;
	int i, slen;
	double	ttwist, trise;
	// string loup[ hashed ];

	// loup["g"] = "G"; loup["a"] = "A";
	// loup["t"] = "T"; loup["c"] = "C"; loup["u"] = "U";

	has_s = 1; has_a = 1;

	if( seq == NULL && aseq == NULL ){
		fprintf( stderr, "wc_helix: no sequence\n" );
		return( NULL );
	}else if( seq == NULL ){
		seq = wc_complement( aseq, snatype );
		has_s = 0;
	}else if( aseq == NULL ){
		aseq = wc_complement( seq, anatype );
		has_a = 0;
	}
#if 0

	slen = strlen( seq );

	setreslibkind( sreslib_use, snatype );
	srname = loup[ substr( seq, 1, 1 ) ];
	if( !strcmp(snatype, "dna"))  srname = "D" + loup[ substr( seq, 1, 1 ) ];
	if( opts =~ "s5" )
		sres = getresidue( srname + "5", sreslib_use );
	else if( opts =~ "s3" && slen == 1 )
		sres = getresidue( srname + "3", sreslib_use );
	else sres = getresidue( srname, sreslib_use );

	setreslibkind( areslib_use, anatype );
	arname = loup[ substr( aseq, 1, 1 ) ];
	if( anatype == "dna") arname = "D" + loup[ substr( aseq, 1, 1 ) ];
	if( opts =~ "a3" )
		ares = getresidue( arname + "3", areslib_use );
	else if( opts =~ "a5" && slen == 1 )
		ares = getresidue( arname + "5", areslib_use );
	else ares = getresidue( arname, areslib_use );
	m1 = wc_basepair( sres, ares );
	freeresidue( sres );
	freeresidue( ares );

	xomat = newtransform(xoff, 0., 0., 0., 0., 0. );
	transformmol( xomat, m1, NULL );
	inmat = newtransform( 0., 0., 0., incl, 0., 0.);
	transformmol( inmat, m1, NULL );

	trise = rise; ttwist = twist;
	for( i = 2; i <= slen-1; i = i + 1 ){
		srname = loup[ substr( seq, i, 1 ) ];
		if( snatype == "dna" ) srname = "D" + loup[ substr( seq, i, 1 ) ];
		setreslibkind( sreslib_use, snatype );
		sres = getresidue( srname, sreslib_use );
		arname = loup[ substr( aseq, i, 1 ) ];
		if( anatype == "dna" ) arname = "D" + loup[ substr( aseq, i, 1 ) ];
		setreslibkind( areslib_use, anatype );
		ares = getresidue( arname, areslib_use );
		m2 = wc_basepair( sres, ares );
		freeresidue( sres );
		freeresidue( ares );
		transformmol( xomat, m2, NULL );
		transformmol( inmat, m2, NULL );
		mat = newtransform( 0., 0., trise,
			0., 0., ttwist );
		transformmol( mat, m2, NULL );
		mergestr( m1, "sense", "last",
			m2, "sense", "first" );
		connectres( m1, "sense",
			i-1, "O3'", i, "P" );
		mergestr( m1, "anti", "first",
			m2, "anti", "last" );
		connectres( m1, "anti",
			1, "O3'", 2, "P" );
		trise = trise + rise;
		ttwist = ttwist + twist;
		freemolecule( m2 );
	}

	i = slen;         // add in final residue pair
	if( i > 1 ){

		srname = loup[ substr( seq, i, 1 ) ];
		if( snatype == "dna" ) srname = "D" + loup[ substr( seq, i, 1 ) ];
		setreslibkind( sreslib_use, snatype );
		if( opts =~ "s3"  )
			sres = getresidue( srname + "3", sreslib_use );
		else 
			sres = getresidue( srname, sreslib_use );
		arname = loup[ substr( aseq, i, 1 ) ];
		if( anatype == "dna" ) arname = "D" + loup[ substr( aseq, i, 1 ) ];
		setreslibkind( areslib_use, anatype );
		if( opts =~ "a5" )
			ares = getresidue( arname + "5", areslib_use );
		else 
			ares = getresidue( arname, areslib_use );

		m2 = wc_basepair( sres, ares );
		freeresidue( sres );
		freeresidue( ares );
		transformmol( xomat, m2, NULL );
		transformmol( inmat, m2, NULL );
		mat = newtransform( 0., 0., trise,
			0., 0., ttwist );
		transformmol( mat, m2, NULL );
		mergestr( m1, "sense", "last",
			m2, "sense", "first" );
		connectres( m1, "sense",
			i-1, "O3'", i, "P" );
		mergestr( m1, "anti", "first",
			m2, "anti", "last" );
		connectres( m1, "anti",
			1, "O3'", 2, "P" );
		trise = trise + rise;
		ttwist = ttwist + twist;
		freemolecule( m2 );
	}

	m3 = newmolecule();
	addstrand( m3, "sense" );
	addstrand( m3, "anti" );
	if( has_s )
		mergestr( m3, "sense", "last", m1, "sense", "first" );
	if( has_a )
		mergestr( m3, "anti", "last", m1, "anti", "first" );

	freemolecule( m1 );

#endif
	return( m3 );
};
