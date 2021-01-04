#include    "nabc.h"

FILE *nabout;    /* can this go into nabc.h?  */

int main( int arg, char *argv[] )
{

    MOLECULE_T *m;
    STRAND_T *sp;    
    RESIDUE_T *rp;
    ATOM_T *ap;

    m = getpdb("2.pdb", NULL);
    upd_molnumbers( m );

    //  strands are stored in a simple linked-list:
    for (sp = m->m_strands; sp; sp = sp->s_next) {

        printf( "strand %d:  %s\n", sp->s_strandnum, sp->s_strandname );

        // but residues are stored in a residue array inside each strand:
        for (int rn = 0; rn < sp->s_nresidues; rn++) {
            rp = sp->s_residues[rn];

            printf( "    residue %d: %s\n", rp->r_resnum, rp->r_resname );

            // atoms are stored in an atom array inside each residue:
            for (int an = 0; an < rp->r_natoms; an++) {
                ap = &rp->r_atoms[an];

                NAB_arc( ap, "fullname" );  // needed to update a_fullname
                // printf( "        atom %d: %s\n", an, ap->a_atomname );

            }                   /* end loop over atoms in this residue  */

        }                       /* end loop over residues in this strand  */

    }                           /* end loop over strands  */

}
