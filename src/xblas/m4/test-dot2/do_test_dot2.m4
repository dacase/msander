dnl Generates test code for dot2
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_DOT2_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_DOT2_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on DOT2  
 *
 * Arguments
 * =========
 *
 * n         (input) int
 *           The size of vector being tested
 *
 * ntests    (input) int
 *           The number of tests to run for each set of attributes.
 *
 * seed      (input/output) int         
 *           The seed for the random number generator used in testgen().
 *
 * thresh    (input) double
 *           When the ratio returned from test() exceeds the specified
 *           threshold, the current size, r_true, r_comp, and ratio will be
 *           printed.  (Since ratio is supposed to be O(1), we can set thresh
 *           to ~10.)
 *
 * debug     (input) int
 *           If debug=3, print summary 
 *           If debug=2, print summary only if the number of bad ratios > 0
 *           If debug=1, print complete info if tests fail
 *           If debug=0, return max ratio
 *
 * min_ratio (output) double
 *           The minimum ratio
 * 
 * num_bad_ratio (output) int
 *               The number of tests fail; they are above the threshold.
 *
 * num_tests (output) int
 *           The number of tests is being performed.
 *
 * Return value
 * ============
 *
 * The maximum ratio if run successfully, otherwise return -1 
 *
 */')dnl
dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_DOT2(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abr_typeltr : the type and precision of alpha, beta and r
dnl        x_typeltr   : the type and precision of x
dnl        y_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
define(`DO_TEST_DOT2', 
  `double DO_TEST_DOT2_NAME($1, $2, $3, $4)(int n, int ntests, dnl 
       int *seed, double thresh, int debug, float test_prob, dnl
       double *min_ratio, int *num_bad_ratio, int *num_tests)
   DO_TEST_DOT2_COMMENT($4)
   DO_TEST_DOT2_BODY($1, $2, $3, $4)')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_DOT2_BODY(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abr_typeltr : the type and precision of alpha, beta and r
dnl        x_typeltr   : the type and precision of x
dnl        y_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_DOT2_BODY',
`{                                        
  /* function name */
  const char fname[] = "DOT2_NAME($1, $2, $3, $4)";

  /* max number of debug lines to print */
  const int max_print = 32;

  /* Variables in the "x_val" form are loop vars for corresponding
     variables */
  int i;            /* iterate through the repeating tests */
  int incx_val, incy_val; /* for testing different inc values */
  int incx_gen, incy_gen; /* for complex case inc=2, for real case inc=1 */
  int d_count;      /* counter for debug */
  int find_max_ratio; /* find_max_ratio = 1 only if debug = 3 */
  int p_count;      /* counter for the number of debug lines printed*/
  int tot_tests;    /* total number of tests to be done */
  int norm;         /* input values of near underflow/one/overflow */
  double ratio_max; /* the current maximum ratio */
  double ratio_min; /* the current minimum ratio */
  double ratio;     /* the per-use test ratio from test() */
  int bad_ratios = 0;   /* the number of ratios over the threshold */
  double eps_int;   /* the internal epsilon expected--2^(-24) for float */
  double un_int;    /* the internal underflow threshold */
  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(x, $2_type)
  DECLARE_VECTOR(head_y, $3_type)
  DECLARE_VECTOR(tail_y, $3_type)
  
  /* x_gen and y_gen are used to store vectors generated by testgen.
     They are eventually copied into x and y for real runs */
  DECLARE_VECTOR(x_gen, $2_type)
  DECLARE_VECTOR(head_y_gen, $3_type) 
  DECLARE_VECTOR(tail_y_gen, $3_type) 

  /* the true r calculated by testgen(), in double-double */
  DECLARE(r_true, EXTRA_TYPE($1_type, $2_type, $3_type))

  DECLARE(r, $1_type)   /* the generated r */
  DECLARE(r_comp, $1_type)  /* the r computed  by DOT2_NAME($1, $2, $3, $4) */
  int alpha_val;
  int alpha_flag;   /* input flag for DOT2_TESTGEN_NAME($1, $2, $3, $4) */
  int beta_val;
  int beta_flag;    /* input flag for DOT2_TESTGEN_NAME($1, $2, $3, $4) */
  int conj_val;
  enum blas_conj_type conj_type;
  ifelse(`$4', `_x', `int prec_val;')
  enum blas_prec_type prec;
  int saved_seed;   /* for saving the original seed */
  int count = 0, old_count = 0;  /* used for counting the number of
                                    testgen calls * 2 */ 

  FPU_FIX_DECL;

  /* test for bad arguments */
  if (n < 0 )
    BLAS_error(fname,  -1,  n, NULL);   
  else if ( ntests < 0 )
    BLAS_error(fname,  -2,  ntests, NULL);   

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0){
    *min_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return 0.0;
  }

  FPU_FIX_START;

  /* initialization */
  saved_seed = *seed;
  ratio_min = 1e308;
  ratio_max = 0.0;
  tot_tests = 0;
  p_count = 0;
  count = 0;
  find_max_ratio = 0;
  if (debug == 3) find_max_ratio = 1;
  incx_gen = incy_gen = 1;
  INC_ADJUST(incx_gen, $2_type)
  INC_ADJUST(incy_gen, $3_type)

  /* get space for calculation */
  MALLOC_VECTOR(x, $2_type, n*2*incx_gen)
dnl  MALLOC_VECTOR(y, $3_type, n*2*incy_gen)
  MALLOC_VECTOR(head_y, $3_type, n*2*incy_gen)
  MALLOC_VECTOR(tail_y, $3_type, n*2*incy_gen)
  MALLOC_VECTOR(x_gen, $2_type, n*incx_gen)
  MALLOC_VECTOR(head_y_gen, $3_type, n*incy_gen)
  MALLOC_VECTOR(tail_y_gen, $3_type, n*incy_gen)

  /* The debug iteration:
     If debug=1, then will execute the iteration twice. First, compute the
     max ratio. Second, print info if ratio > (50% * ratio_max). */
  for (d_count=0; d_count<= find_max_ratio; d_count++) {
    bad_ratios = 0; /* set to zero */ 

    if ((debug == 3) && (d_count == find_max_ratio))
      *seed = saved_seed; /* restore the original seed */

    /* varying alpha */
    for (alpha_val=1; alpha_val<3; alpha_val++) { 
      SET_ALPHA($1_type)
        
      /* varying beta */
      for (beta_val=0; beta_val<3; beta_val++) {
        SET_BETA($1_type)

        ifelse($4, _x, `
        /* varying extra precs */
        for (prec_val = 0; prec_val <= 2; prec_val++) {')
          SET_INTERNAL_PARAMS($1_type, $4)
           
          /* values near underflow, 1, or overflow */
          for (norm = -1; norm <= 1; norm++) {

            /* number of tests */
            for (i=0; i<ntests; i++){           

              /* conj or not */
              for (conj_val=0; conj_val<2; conj_val++) {
                switch(conj_val){
                case 0: conj_type = blas_no_conj; break;
                case 1: 
                default: conj_type = blas_conj; break;
                }

                /* For the sake of speed, we throw out this case at random */
                if ( xrand(seed) >= test_prob ) continue;

                DOT2_TESTGEN_NAME($1, $2, $3)(n, 0, 0, norm, conj_type, dnl
                             &alpha, alpha_flag, &beta, beta_flag, dnl
                             head_y_gen, tail_y_gen, x_gen, seed, &r, dnl
                             PASS_BY_REF(r_true, EXTRA_TYPE($1_type)));
                count++;

                /* varying incx */
                for (incx_val = -2; incx_val <= 2; incx_val++){
                  if (incx_val == 0) continue;         

                  $2copy_vector(x_gen, n, 1, x, incx_val);

                  /* varying incy */
                  for (incy_val = -2; incy_val <= 2; incy_val++){
                    if (incy_val == 0) continue;

                    $3copy_vector(head_y_gen, n, 1, head_y, incy_val);
                    $3copy_vector(tail_y_gen, n, 1, tail_y, incy_val);

                    /* call DOT2_NAME($1, $2, $3, $4) to get r_comp */
                    ASSIGN(r_comp, $1_type, r, $1_type)

                    FPU_FIX_STOP;
                    DOT2_NAME($1, $2, $3, $4)(conj_type, n, alpha, x, dnl
                        incx_val, beta, head_y, tail_y, incy_val, dnl
                        &r_comp ifelse(`$4', `_x', `, prec'));
                    FPU_FIX_START;

                    /* computing the ratio */
                    TEST_DOT2_NAME($1, $2, $3, $4)(n, conj_type, alpha, dnl
                        beta, r, r_comp, HEAD(r_true), TAIL(r_true), x, dnl
                        incx_val, head_y, tail_y, incy_val, eps_int, dnl
                        un_int, &ratio);

                    /* Increase the number of bad ratios, if the ratio
                       is bigger than the threshold.
                       The !<= below causes NaN error to be detected.
                       Note that (NaN > thresh) is always false. */
                    if ( !(ratio <= thresh) ) {
                      bad_ratios++;
                  
                      if ((debug == 3) &&        /* print only when debug is on */
                         (count != old_count) && /* print if old vector is different 
                                                    from the current one */ 
                         (d_count == find_max_ratio) &&  
                         (p_count <= max_print) &&
                         (ratio > 0.5*ratio_max))
                      {
                        old_count = count;      

                        printf("FAIL> %s: n = %d, ntests = %d, threshold = %4.2f,\n", 
                                fname, n, ntests, thresh);
        
                        /* Print test info */
                        PRINT_PREC(prec)
                        PRINT_NORM(norm)
                        PRINT_CONJ(conj_type)

                        printf("incx=%d, incy=%d:\n", incx_val, incy_val);
                        
                        $2print_vector(x, n, incx_val, "x");
                        $3print_vector(head_y, n, incx_val, "head_y");
                        $3print_vector(tail_y, n, incx_val, "tail_y");

                        printf("      "); PRINT_VAR(alpha, $1_type) 
                        printf("\n      "); 
                        PRINT_VAR(beta, $1_type) printf("\n");
                        printf("      "); PRINT_VAR(r, $1_type) 
                        printf("\n      "); 
                        PRINT_VAR(r_comp, $1_type) printf("\n");
                        printf("      "); PRINT_VAR(r_true, EXTRA_TYPE($1_type))
                        printf("      ratio=%.4e\n", ratio);
                        p_count++;
                      }
                    }

                    if (d_count == 0) {

                      if (ratio > ratio_max)
                        ratio_max = ratio;

                      if (ratio != 0.0 && ratio < ratio_min)
                        ratio_min = ratio;

                      tot_tests++;
                    }
                  } /* incy */
                } /* incx */
              } /* conj */
            } /* tests */
          } /* norm */
ifelse(`$4', `_x', `         } /* prec */')
      } /* beta */
    } /* alpha */
  } /* debug */

  if ((debug == 2) ||
     ((debug == 1) && (bad_ratios > 0))){
    printf("      %s:  n = %d, ntests = %d, thresh = %4.2f\n", 
            fname, n, ntests, thresh);
    printf("      bad/total = %d/%d=%3.2f, min_ratio = %.4e, max_ratio = %.4e\n\n",
            bad_ratios, tot_tests,
           ((double)bad_ratios)/((double)tot_tests), ratio_min, ratio_max);
  }

  FREE_VECTOR(x, $2_type)
  FREE_VECTOR(head_y, $3_type)
  FREE_VECTOR(tail_y, $3_type)
  FREE_VECTOR(x_gen, $2_type)
  FREE_VECTOR(head_y_gen, $3_type)
  FREE_VECTOR(tail_y_gen, $3_type)

  *min_ratio = ratio_min;
  *num_bad_ratio = bad_ratios;
  *num_tests = tot_tests;

  FPU_FIX_STOP;
  return ratio_max;
}' )dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DOT2_NAME(abr_typeltr, x_typeltr, y_typeltr, extended) 
dnl        create a dot2 name 
dnl --------------------------------------------------------------------
define(`DOT2_NAME', `ifelse(
        `$2&&$3', `$1&&$1',
        `BLAS_$1dot2$4',
        `BLAS_$1dot2_$2_$3$4')') dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TEST_DOT2_NAME(abr_typeltr, x_typeltr, y_typeltr, extended) 
dnl        create a test_dot2 name 
dnl --------------------------------------------------------------------
define(`TEST_DOT2_NAME', `ifelse(
        `$2&&$3', `$1&&$1',
        `test_BLAS_$1dot2',
        `test_BLAS_$1dot2_$2_$3')')dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: DO_TEST_DOT2_NAME(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl        create do_test_dot2 name
dnl -------------------------------------------------------------------
define(`DO_TEST_DOT2_NAME', `ifelse(
        `$2&&$3', `$1&&$1',
        `do_test_$1dot2$4',
        `do_test_$1dot2_$2_$3$4')') dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: CALL_DO_TEST_DOT2(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abr_type : the type and precision of alpha, beta and r
dnl        x_type   : the type and precision of x
dnl        y_type   : the type and precision of y
dnl        extended : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl -------------------------------------------------------------------
define(`CALL_DO_TEST_DOT2',           
  `min_ratio = 1e308; max_ratio = 0.0;
   total_bad_ratios = 0; total_tests = 0;
   fname = "DOT2_NAME($1, $2, $3, $4)";
   printf("Testing %s...\n", fname);
   for(n=0; n<=nsizes; n++){
    
     total_max_ratio = DO_TEST_DOT2_NAME($1, $2, $3, $4)(n, ntests, dnl
         &seed, thresh, debug, test_prob, &total_min_ratio, dnl
         &num_bad_ratio, &num_tests);
     if (total_max_ratio > max_ratio)
       max_ratio = total_max_ratio;

     if (total_min_ratio != 0.0 && total_min_ratio < min_ratio)
       min_ratio = total_min_ratio;

     total_bad_ratios += num_bad_ratio;
     total_tests += num_tests; 
   }

   if (total_bad_ratios == 0)
     printf("PASS> ");
   else {
     printf("FAIL> ");
     nr_failed_routines++;
   }
   nr_routines++;

   if (min_ratio == 1e308) min_ratio = 0.0;

   printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n", dnl
       fname, total_bad_ratios, total_tests, max_ratio);
')
dnl
dnl
DO_TEST_DOT2(s, s, s, _x)
DO_TEST_DOT2(d, d, d, _x)
DO_TEST_DOT2(c, c, c, _x)
DO_TEST_DOT2(z, z, z, _x)

MAIN(`', `
CALL_DO_TEST_DOT2(s, s, s, _x)
CALL_DO_TEST_DOT2(d, d, d, _x)
CALL_DO_TEST_DOT2(c, c, c, _x)
CALL_DO_TEST_DOT2(z, z, z, _x)
')dnl
dnl
dnl