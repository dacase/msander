/* gaussian esp generated by iop(6/50=1) */

#include "eprintf.h"


int rgesp(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo, MOLINFO minfo)
{
    FILE *fpin;
    int rflag = 0;
    int numatom;
    int overflow_flag = 0;
    char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR];
    char tmpchar4[MAXCHAR], tmpchar5[MAXCHAR];
    char line[MAXCHAR];


    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, minfo.resname);
    numatom = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
/*     printf("\nFinished reading %s file.", cinfo.ifilename); */
            break;
        }
        sscanf(line, "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5);

        if (rflag == 1 && strcmp(tmpchar1, "DIPOLE") == 0)
            break;
        if (rflag == 0 && strncmp(line, " ATOMIC COORDINATES", 19) == 0) {
            rflag = 1;
            continue;
        }

        if (rflag == 1) {
            if (overflow_flag == 0) {
                strcpy(atom[numatom].name, tmpchar1);
                strcpy(atom[numatom].element, tmpchar1);
                atom[numatom].x = translate(tmpchar2);
                atom[numatom].y = translate(tmpchar3);
                atom[numatom].z = translate(tmpchar4);
                atom[numatom].x *= Bohr;
                atom[numatom].y *= Bohr;
                atom[numatom].z *= Bohr;
                atom[numatom].charge = translate(tmpchar5);
            }
            numatom++;
            if (numatom >= cinfo.maxatom && overflow_flag == 0) {
                printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", numatom);
                overflow_flag = 1;
            }
        }
    }
    fclose(fpin);
    *atomnum = numatom;
    return overflow_flag;
}

void wgesp()
{
    printf("Warning: Cannot write a Gaussian ESP file.\n"
           "         You must run Gaussian 09 with iop(6/50=1).\n");
}