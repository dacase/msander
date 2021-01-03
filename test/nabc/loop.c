#include    "nabc.h"

FILE *nabout;    /* can this go into nabc.h?  */

int main( int arg, char *argv[] )
{

    int rn, an;
    MOLECULE_T *m;
    STRAND_T *s;
    RESIDUE_T *r;
    ATOM_T *a;

    m = getpdb("2.pdb", NULL);

    //  strands are stored in a simple linked-list:
    for (s = m->m_strands; s; s = s->s_next) {

        printf( "strand %d:  %s\n", s->s_strandnum, s->s_strandname );

        // but residues are stored in a residue array inside each strand:
        for (int rn = 0; rn < s->s_nresidues; rn++) {
            r = s->s_residues[rn];

            printf( "    residue %d: %s\n", rn, r->r_resname );

            // atoms are stored in an atom array inside each residue:
            for (int an = 0; an < r->r_natoms; an++) {
                a = &r->r_atoms[an];

                NAB_arc( a, "fullname" );  // needed to update a_fullname
                printf( "        atom %d: %s\n", an, a->a_fullname );

            }                   /* end loop over atoms in this residue  */

        }                       /* end loop over residues in this strand  */

    }                           /* end loop over strands  */

}
