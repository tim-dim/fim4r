/*----------------------------------------------------------------------
  File    : fim4r.c
  Contents: frequent item set mining & association rule induction for R
  Author  : Christian Borgelt
  History : 2016.03.30 file created
            2016.03.31 functions genpsp(), estpsp(), patred() completed
            2016.04.05 function patred() redesigned to work on pairs
            2016.04.06 first version of all mining functions completed
            2016.11.05 adapted to modified apriori interface
            2016.11.11 adapted to modified eclat interface
            2016.11.15 adapted to modified accretion interface
            2016.11.21 adapted to modified fpgrowth interface
            2017.03.25 adapted to modified carpenter/ista interfaces
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include "R.h"
#include "Rinternals.h"
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <pthread.h>
#endif
#ifndef TATREEFN
#define TATREEFN
#endif
#ifndef TA_SURR
#define TA_SURR
#endif
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#ifndef PSP_ESTIM
#define PSP_ESTIM
#endif
/*-#include "sigint.h"--*/
#include "report.h"
#include "apriori.h"
#include "eclat.h"
#include "fpgrowth.h"
#include "sam.h"
#include "relim.h"
#include "carpenter.h"
#include "ista.h"
#include "accretion.h"
#include "fpgpsp.h"
#include "patred.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define BLKSIZE         1024    /* block size for arrays */
#define IS_NA(x)        (isnan(x) || isinf(x) || ((x) < 0))

/* --- error handling --- */
#define MYERROR(msg)    do { sig_remove(); error(msg); } while (0)
#define ERR_MEM()       MYERROR("out of memory")
#define ERR_ABORT()     MYERROR("user abort")

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- item set report data --- */
  SEXP   res;                   /* constructed result object */
  size_t size;                  /* size     of result array */
  size_t cnt;                   /* elements in result array */
  int    istr;                  /* item type (integer or string) */
  int    len;                   /* number     of values to report */
  CCHAR  *rep;                  /* indicators of values to report */
  int    err;                   /* error flag */
} REPDATA;                      /* (item set report data) */

/*----------------------------------------------------------------------
  Parameter Functions
----------------------------------------------------------------------*/

static int get_lgl (SEXP p, int dflt)
{                               /* --- get a logical parameter */
  assert(p);                    /* check the function argument */
  if (length(p) <  1)       return dflt;
  if (TYPEOF(p) == LGLSXP)  return (LOGICAL(p)[0])      ? -1 : 0;
  if (TYPEOF(p) == INTSXP)  return (INTEGER(p)[0] != 0) ? -1 : 0;
  if (TYPEOF(p) == REALSXP) return (REAL(p)[0]    != 0) ? -1 : 0;
  return dflt;                  /* extract and convert first element */
}  /* get_lgl() */

/*--------------------------------------------------------------------*/

static int get_int (SEXP p, int dflt)
{                               /* --- get an integer parameter */
  assert(p);                    /* check the function argument */
  if (length(p) <  1)       return dflt;
  if (TYPEOF(p) == LGLSXP)  return (LOGICAL(p)[0]) ? 1 : 0;
  if (TYPEOF(p) == INTSXP)  return (int)INTEGER(p)[0];
  if (TYPEOF(p) == REALSXP) return (int)REAL(p)[0];
  return dflt;                  /* extract and convert first element */
}  /* get_int() */

/*--------------------------------------------------------------------*/

static long int get_lng (SEXP p, long int dflt)
{                               /* --- get an integer parameter */
  assert(p);                    /* check the function argument */
  if (length(p) <  1)       return dflt;
  if (TYPEOF(p) == LGLSXP)  return (LOGICAL(p)[0]) ? 1 : 0;
  if (TYPEOF(p) == INTSXP)  return (long int)INTEGER(p)[0];
  if (TYPEOF(p) == REALSXP) return (long int)REAL(p)[0];
  return dflt;                  /* extract and convert first element */
}  /* get_lng() */

/*--------------------------------------------------------------------*/

static double get_dbl (SEXP p, double dflt)
{                               /* --- get a double parameter */
  assert(p);                    /* check the function argument */
  if (length(p) <  1)       return dflt;
  if (TYPEOF(p) == LGLSXP)  return (LOGICAL(p)[0]) ? 1.0 : 0.0;
  if (TYPEOF(p) == INTSXP)  return (double)INTEGER(p)[0];
  if (TYPEOF(p) == REALSXP) return (double)REAL(p)[0];
  return dflt;                  /* extract and convert first element */
}  /* get_dbl() */

/*--------------------------------------------------------------------*/

static CCHAR* get_str (SEXP p, CCHAR *dflt)
{                               /* --- get a string parameter */
  assert(p);                    /* check the function argument */
  if (length(p) <  1)       return dflt;
  if (TYPEOF(p) == STRSXP)  return CHAR(STRING_ELT(p, 0));
  return dflt;                  /* extract and convert first element */
}  /* get_str() */

/*--------------------------------------------------------------------*/

static int get_target (SEXP p, const char *targets)
{                               /* --- get target */
  CCHAR *s;                     /* target string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != STRSXP) return -1;
  if (length(p) < 1) s = "s";   /* get the target string */
  else s = CHAR(STRING_ELT(p, 0));
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "set")        == 0) s = "s";
    else if (strcmp(s, "sets")       == 0) s = "s";
    else if (strcmp(s, "all")        == 0) s = "s";
    else if (strcmp(s, "allset")     == 0) s = "s";
    else if (strcmp(s, "allsets")    == 0) s = "s";
    else if (strcmp(s, "frq")        == 0) s = "s";
    else if (strcmp(s, "freq")       == 0) s = "s";
    else if (strcmp(s, "frequent")   == 0) s = "s";
    else if (strcmp(s, "frqset")     == 0) s = "s";
    else if (strcmp(s, "frqsets")    == 0) s = "s";
    else if (strcmp(s, "freqset")    == 0) s = "s";
    else if (strcmp(s, "freqsets")   == 0) s = "s";
    else if (strcmp(s, "cls")        == 0) s = "c";
    else if (strcmp(s, "clsd")       == 0) s = "c";
    else if (strcmp(s, "closed")     == 0) s = "c";
    else if (strcmp(s, "max")        == 0) s = "m";
    else if (strcmp(s, "maxi")       == 0) s = "m";
    else if (strcmp(s, "maximal")    == 0) s = "m";
    else if (strcmp(s, "gen")        == 0) s = "g";
    else if (strcmp(s, "gens")       == 0) s = "g";
    else if (strcmp(s, "generas")    == 0) s = "g";
    else if (strcmp(s, "generators") == 0) s = "g";
    else if (strcmp(s, "rule")       == 0) s = "r";
    else if (strcmp(s, "rules")      == 0) s = "r";
    else if (strcmp(s, "arule")      == 0) s = "r";
    else if (strcmp(s, "arules")     == 0) s = "r";
  }
  if (s[0] && !s[1]             /* check for a valid string */
  && (strchr(targets, s[0]) != NULL)) {
    switch (s[0]) {             /* evaluate the target code */
      case 's': return ISR_SETS;
      case 'a': return ISR_SETS;
      case 'f': return ISR_SETS;
      case 'c': return ISR_CLOSED;
      case 'm': return ISR_MAXIMAL;
      case 'g': return ISR_GENERAS;
      case 'r': return ISR_RULES;
    }
  }
  return -1;                    /* return an error code */
}  /* get_target() */

/*--------------------------------------------------------------------*/

static int get_stat (SEXP p)
{                               /* --- get statistic code */
  CCHAR *s;                     /* surrogate string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != STRSXP) return -1;
  if (length(p) < 1) s = "s";   /* get the target string */
  else s = CHAR(STRING_ELT(p, 0));
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "none")      == 0) s = "x";
    else if (strcmp(s, "X2")        == 0) s = "p";
    else if (strcmp(s, "chi2")      == 0) s = "p";
    else if (strcmp(s, "X2pval")    == 0) s = "p";
    else if (strcmp(s, "chi2pval")  == 0) s = "p";
    else if (strcmp(s, "yates")     == 0) s = "t";
    else if (strcmp(s, "yatespval") == 0) s = "t";
    else if (strcmp(s, "info")      == 0) s = "g";
    else if (strcmp(s, "infopval")  == 0) s = "g";
    else if (strcmp(s, "gpval")     == 0) s = "g";
    else if (strcmp(s, "fetprob")   == 0) s = "f";
    else if (strcmp(s, "fetX2")     == 0) s = "h";
    else if (strcmp(s, "fetchi2")   == 0) s = "h";
    else if (strcmp(s, "fetinfo")   == 0) s = "m";
    else if (strcmp(s, "fetsupp")   == 0) s = "s";
  }
  if (s[0] && !s[1]) {          /* check for a valid string */
    switch (s[0]) {             /* evaluate the statistic code */
      case 'x': return RE_NONE;
      case 'c': return RE_CHI2PVAL;
      case 'p': return RE_CHI2PVAL;
      case 'n': return RE_CHI2PVAL;
      case 'y': return RE_YATESPVAL;
      case 't': return RE_YATESPVAL;
      case 'i': return RE_INFOPVAL;
      case 'g': return RE_INFOPVAL;
      case 'f': return RE_FETPROB;
      case 'h': return RE_FETCHI2;
      case 'm': return RE_FETINFO;
      case 's': return RE_FETSUPP;
    }
  }
  return -1;                    /* return an error code */
}  /* get_stat() */

/*--------------------------------------------------------------------*/

static int get_eval (SEXP p)
{                               /* --- get evaluation measure code */
  CCHAR *s;                     /* surrogate string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != STRSXP) return -1;
  if (length(p) < 1) s = "s";   /* get the target string */
  else s = CHAR(STRING_ELT(p, 0));
  if (s[0] && s[1]) {
    if (strcmp(s, "none")    == 0) return 'x';
    if (strcmp(s, "ldratio") == 0) return 'b';
  }
  if (strchr("xb", s[0]) != NULL) return s[0];
  return -1;                    /* return an error code */
}  /* get_eval() */

/*--------------------------------------------------------------------*/

static int get_evalx (SEXP p)
{                               /* --- get evaluation measure code */
  CCHAR *s;                     /* surrogate string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != STRSXP) return -1;
  if (length(p) < 1) s = "s";   /* get the target string */
  else s = CHAR(STRING_ELT(p, 0));
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if (strcmp(s, "none")       == 0) s = "x";
    if (strcmp(s, "supp")       == 0) s = "o";
    if (strcmp(s, "support")    == 0) s = "o";
    if (strcmp(s, "conf")       == 0) s = "c";
    if (strcmp(s, "confidence") == 0) s = "c";
    if (strcmp(s, "confdiff")   == 0) s = "d";
    if (strcmp(s, "lift")       == 0) s = "l";
    if (strcmp(s, "liftdiff")   == 0) s = "a";
    if (strcmp(s, "liftquot")   == 0) s = "q";
    if (strcmp(s, "cvct")       == 0) s = "v";
    if (strcmp(s, "conviction") == 0) s = "v";
    if (strcmp(s, "cvctdiff")   == 0) s = "e";
    if (strcmp(s, "cvctquot")   == 0) s = "r";
    if (strcmp(s, "cprob")      == 0) s = "k";
    if (strcmp(s, "import")     == 0) s = "j";
    if (strcmp(s, "importance") == 0) s = "j";
    if (strcmp(s, "cert")       == 0) s = "z";
    if (strcmp(s, "chi2")       == 0) s = "n";
    if (strcmp(s, "X2")         == 0) s = "n";
    if (strcmp(s, "chi2pval")   == 0) s = "p";
    if (strcmp(s, "X2pval")     == 0) s = "p";
    if (strcmp(s, "yates")      == 0) s = "y";
    if (strcmp(s, "yatespval")  == 0) s = "t";
    if (strcmp(s, "info")       == 0) s = "i";
    if (strcmp(s, "infopval")   == 0) s = "g";
    if (strcmp(s, "gpval")      == 0) s = "g";
    if (strcmp(s, "fetprob")    == 0) s = "f";
    if (strcmp(s, "fetchi2")    == 0) s = "h";
    if (strcmp(s, "fetX2")      == 0) s = "h";
    if (strcmp(s, "fetinfo")    == 0) s = "m";
    if (strcmp(s, "fetsupp")    == 0) s = "s";
    if (strcmp(s, "ldratio")    == 0) s = "b";
  }
  if (s[0] && !s[1]) {          /* check for a valid string */
    switch (s[0]) {             /* evaluate the measure code */
      case 'x': return RE_NONE;
      case 'o': return RE_SUPP;
      case 'c': return RE_CONF;
      case 'd': return RE_CONFDIFF;
      case 'l': return RE_LIFT;
      case 'a': return RE_LIFTDIFF;
      case 'q': return RE_LIFTQUOT;
      case 'v': return RE_CVCT;
      case 'e': return RE_CVCTDIFF;
      case 'r': return RE_CVCTQUOT;
      case 'k': return RE_CPROB;
      case 'j': return RE_IMPORT;
      case 'z': return RE_CERT;
      case 'n': return RE_CHI2;
      case 'p': return RE_CHI2PVAL;
      case 'y': return RE_YATES;
      case 't': return RE_YATESPVAL;
      case 'i': return RE_INFO;
      case 'g': return RE_INFOPVAL;
      case 'f': return RE_FETPROB;
      case 'h': return RE_FETCHI2;
      case 'm': return RE_FETINFO;
      case 's': return RE_FETSUPP;
      case 'b': return RE_FNCNT;
    }
  }
  return -1;                    /* return an error code */
}  /* get_evalx() */

/*--------------------------------------------------------------------*/

static int get_agg (SEXP p)
{                               /* --- get aggregation mode */
  CCHAR *s;                     /* surrogate string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != STRSXP) return -1;
  if (length(p) < 1) s = "p";   /* get the surrogate function string */
  else s = CHAR(STRING_ELT(p, 0));
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "none")    == 0) s = "x";
    else if (strcmp(s, "min")     == 0) s = "m";
    else if (strcmp(s, "minimum") == 0) s = "m";
    else if (strcmp(s, "max")     == 0) s = "n";
    else if (strcmp(s, "maximum") == 0) s = "n";
    else if (strcmp(s, "avg")     == 0) s = "a";
    else if (strcmp(s, "average") == 0) s = "a";
  }
  if (s[0] && !s[1]) {          /* check for a valid string */
    switch (s[0]) {             /* evaluate the aggregation code */
      case 'x': return IST_NONE;
      case 'm': return IST_MIN;
      case 'n': return IST_MAX;
      case 'a': return IST_AVG;
    }
  }
  return -1;                    /* return an error code */
}  /* get_agg() */

/*--------------------------------------------------------------------*/

static int get_surr (SEXP p)
{                               /* --- get surrogate function code */
  CCHAR *s;                     /* surrogate string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != STRSXP) return -1;
  if (length(p) < 1) s = "p";   /* get the surrogate function string */
  else s = CHAR(STRING_ELT(p, 0));
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "ident")     == 0) s = "i";
    else if (strcmp(s, "identity")  == 0) s = "i";
    else if (strcmp(s, "random")    == 0) s = "r";
    else if (strcmp(s, "randomize") == 0) s = "r";
    else if (strcmp(s, "swap")      == 0) s = "s";
    else if (strcmp(s, "perm")      == 0) s = "s";
    else if (strcmp(s, "permute")   == 0) s = "s";
    else if (strcmp(s, "shuffle")   == 0) s = "h";
  }
  if (s[0] && !s[1]) {          /* translate surrogate method string */
    switch (s[0]) {             /* evaluate the surrogate method code */
      case 'i': return FPG_IDENTITY;
      case 'r': return FPG_RANDOM;
      case 's': return FPG_SWAP;
      case 'h': return FPG_SHUFFLE;
    }
  }
  return -1;                    /* return an error code */
}  /* get_surr() */

/*--------------------------------------------------------------------*/

static int get_red (SEXP p)
{                               /* --- get random function code */
  CCHAR *s;                     /* density function string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != STRSXP) return -1;
  if (length(p) < 1) s = "S";   /* get the random function string */
  else s = CHAR(STRING_ELT(p, 0));
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "none")        == 0) s = "x";
    else if (strcmp(s, "coins")       == 0) s = "c";
    else if (strcmp(s, "coins0")      == 0) s = "c";
    else if (strcmp(s, "coins1")      == 0) s = "C";
    else if (strcmp(s, "coins+1")     == 0) s = "C";
    else if (strcmp(s, "items")       == 0) s = "i";
    else if (strcmp(s, "items2")      == 0) s = "i";
    else if (strcmp(s, "neurons")     == 0) s = "i";
    else if (strcmp(s, "cover")       == 0) s = "s";
    else if (strcmp(s, "cover0")      == 0) s = "s";
    else if (strcmp(s, "covered")     == 0) s = "s";
    else if (strcmp(s, "covered0")    == 0) s = "s";
    else if (strcmp(s, "cover1")      == 0) s = "S";
    else if (strcmp(s, "covered1")    == 0) s = "S";
    else if (strcmp(s, "leni")        == 0) s = "l";
    else if (strcmp(s, "leni0")       == 0) s = "l";
    else if (strcmp(s, "lenient")     == 0) s = "l";
    else if (strcmp(s, "lenient0")    == 0) s = "l";
    else if (strcmp(s, "leni1")       == 0) s = "L";
    else if (strcmp(s, "lenient1")    == 0) s = "L";
    else if (strcmp(s, "strict")      == 0) s = "t";
    else if (strcmp(s, "strict0")     == 0) s = "t";
    else if (strcmp(s, "strict1")     == 0) s = "T";
  }
  if (s[0] && !s[1]) {          /* translate surrogate method string */
    switch (s[0]) {             /* evaluate the surrogate method code */
      case 'x': return PSR_NONE;
      case 'c': return PSR_COINS0;
      case 'C': return PSR_COINS1;
      case 'i': return PSR_ITEMS2;
      case 's': return PSR_COVER0;
      case 'S': return PSR_COVER1;
      case 'l': return PSR_LENIENT0;
      case 'L': return PSR_LENIENT1;
      case 't': return PSR_STRICT0;
      case 'T': return PSR_STRICT1;
    }
  }
  return -1;                    /* return an error code */
}  /* get_red() */

/*--------------------------------------------------------------------*/

static int get_app (SEXP p)
{                               /* --- get item appearance indicator */
  CCHAR *s;                     /* target string */

  assert(p);                    /* check the function argument */
  if (TYPEOF(p) != CHARSXP) return -1;
  s = CHAR(p);                  /* get the target string */
  if (s[0] && !s[1]) {          /* translate single character */
    switch (s[0]) {             /* evaluate single character */
      case 'n': s = "-"; break;
      case 'i': s = "a"; break;
      case 'b': s = "a"; break;
      case 'o': s = "c"; break;
      case 'h': s = "c"; break;
    } }
  else if (s[0] && s[1]) {      /* evaluate textual identifier */
    if      (strcmp(s, "none")       == 0) s = "-";
    else if (strcmp(s, "neither")    == 0) s = "-";
    else if (strcmp(s, "ign")        == 0) s = "-";
    else if (strcmp(s, "ignore")     == 0) s = "-";
    else if (strcmp(s, "in")         == 0) s = "a";
    else if (strcmp(s, "inp")        == 0) s = "a";
    else if (strcmp(s, "input")      == 0) s = "a";
    else if (strcmp(s, "out")        == 0) s = "c";
    else if (strcmp(s, "output")     == 0) s = "c";
    else if (strcmp(s, "ante")       == 0) s = "a";
    else if (strcmp(s, "antecedent") == 0) s = "a";
    else if (strcmp(s, "cons")       == 0) s = "c";
    else if (strcmp(s, "consequent") == 0) s = "c";
    else if (strcmp(s, "body")       == 0) s = "a";
    else if (strcmp(s, "head")       == 0) s = "c";
    else if (strcmp(s, "io")         == 0) s = "x";
    else if (strcmp(s, "i&o")        == 0) s = "x";
    else if (strcmp(s, "o&i")        == 0) s = "x";
    else if (strcmp(s, "inout")      == 0) s = "x";
    else if (strcmp(s, "in&out")     == 0) s = "x";
    else if (strcmp(s, "ac")         == 0) s = "x";
    else if (strcmp(s, "a&c")        == 0) s = "x";
    else if (strcmp(s, "c&a")        == 0) s = "x";
    else if (strcmp(s, "canda")      == 0) s = "x";
    else if (strcmp(s, "bh")         == 0) s = "x";
    else if (strcmp(s, "b&h")        == 0) s = "x";
    else if (strcmp(s, "h&b")        == 0) s = "x";
    else if (strcmp(s, "both")       == 0) s = "x";
  }
  if (s[0] && !s[1]) {          /* check for a valid string */
    switch (s[0]) {             /* evaluate the appearance code */
      case '-': return APP_NONE;
      case 'a': return APP_BODY;
      case 'c': return APP_HEAD;
      case 'x': return APP_BOTH;
    }
  }
  return -1;                    /* return an error code */
}  /* get_app() */

/*----------------------------------------------------------------------
  Auxiliary Functions
----------------------------------------------------------------------*/

static int chktracts (SEXP tracts, SEXP wgts, SEXP appear)
{                               /* --- check transaction arguments */
  R_xlen_t i, n;                /* loop variables */
  int      t;                   /* item type (integer or string) */

  assert(tracts);               /* check the function arguments */
  if (TYPEOF(tracts) != VECSXP) return -1;
  n = length(tracts);           /* check the data type */
  if (n <= 0) return -1;        /* and the list length */
  t = TYPEOF(VECTOR_ELT(tracts, 0));
  if ((t != INTSXP) && (t != STRSXP))
    return -2;                  /* check the first transaction */
  for (i = 1; i < n; i++)       /* check type of transactions */
    if (TYPEOF(VECTOR_ELT(tracts, i)) != t) return -2;
  if ((wgts != R_NilValue)      /* if transaction weights are given */
  && ((TYPEOF(wgts) != INTSXP) || (length(wgts) != n)))
    return -3;                  /* check type/length of weights array */
  if (appear == R_NilValue)     /* if no item appearances are given, */
    return 0;                   /* simply abort the function */
  if ((TYPEOF(appear) != VECSXP) || (length(appear) != 2)
  ||  ( TYPEOF(VECTOR_ELT(appear, 0)) != t)
  ||  ((TYPEOF(VECTOR_ELT(appear, 1)) != INTSXP)
  &&   (TYPEOF(VECTOR_ELT(appear, 1)) != STRSXP))
  ||  ( length(VECTOR_ELT(appear, 0)) != length(VECTOR_ELT(appear, 1))))
    return -4;                  /* check for proper types and lengths */
  return 0;                     /* return 'ok' */
}  /* chktracts() */

/*--------------------------------------------------------------------*/

static int ib_appRObj (ITEMBASE *ibase, SEXP appear)
{                               /* --- get item appearances */
  int      m;                   /* integer item (from R object) */
  CCHAR    *s;                  /* string  item (from R object) */
  ITEM     i;                   /* item identifier */
  int      app;                 /* item appearance indicator */
  R_xlen_t k, n;                /* loop variables */
  int      t, a;                /* item and app. indicator type */
  SEXP     items;               /* items (first array in appear) */

  assert(ibase);                /* check the function arguments */
  if (appear == R_NilValue)     /* if no item appearances are given, */
    return 0;                   /* simply abort the function */
  items  = VECTOR_ELT(appear,0);/* get the items and their type */
  t = TYPEOF(items);            /* (items are integer or string) */
  appear = VECTOR_ELT(appear,1);/* get the appearance indicators */
  a = TYPEOF(appear);           /* (app. inds. are integer or string) */
  n = length(items);            /* get the array lengths and */
  for (k = 0; k < n; k++) {     /* traverse the items/indicators */
    if (t == INTSXP) {          /* if items are integers */
      m = INTEGER(items)[k];    /* get -1 for NA entries */
      if (m <= INT_MIN) i = -1; /* or add item to item base */
      else if ((i = ib_add(ibase, &m)) < 0) return -1; }
    else {                      /* if items are strings */
      s = CHAR(STRING_ELT(items, k));
      if (s[0] == 0)    i = -1; /* get -1 for empty strings */
      else if ((i = ib_add(ibase, s))  < 0) return -1;
    }                           /* or add item to item base */
    if (a == STRSXP) {          /* if app. indicators are strings */
      app = get_app(STRING_ELT(appear, k));
      if (app < 0) return -1; } /* decode the appearance indicator */
    else {                      /* if app. indicators are integer */
      app = 0; m = INTEGER(appear)[k];
      if (m & 1) app |= APP_BODY;
      if (m & 2) app |= APP_HEAD;
    }                           /* get app. indicators from bits */
    ib_setapp(ibase, i, app);   /* set appearance of item */
  }                             /* (or default appearance if i < 0) */
  return 0;                     /* return 'ok' */
}  /* ib_appRObj() */

/*--------------------------------------------------------------------*/

static TABAG* tbg_fromRObj (SEXP tracts, SEXP wgts, SEXP appear)
{                               /* --- create a transaction bag */
  int      e = 0;               /* error flag */
  TID      k, m;                /* trans. identifier, loop variable */
  ITEM     i, n;                /* item   identifier, loop variable */
  int      t;                   /* item type (integer or string) */
  SEXP     p;                   /* to traverse the transactions */
  ITEMBASE *ibase;              /* underlying item base */
  TABAG    *tabag;              /* created transaction bag */

  assert(tracts);               /* check the function argument */
  t = TYPEOF(VECTOR_ELT(tracts, 0));
  ibase = (t == INTSXP)         /* according to the item type */
        ? ib_create(IB_OBJNAMES, 0, ST_INTFN, (OBJFN*)0)
        : ib_create(0, 0);      /* create an item base */
  if (!ibase) return NULL;      /* for integers or strings */
  if (ib_appRObj(ibase, appear) != 0) { /* add item appearances */
    ib_delete(ibase); return NULL; }    /* to the item base */
  tabag = tbg_create(ibase);    /* create a transaction bag */
  if (!tabag) { ib_delete(ibase); return NULL; }
  m = length(tracts);           /* get the number of transactions */
  for (k = 0; k < m; k++) {     /* and traverse the transactions */
    ib_clear(ibase);            /* clear the internal transaction */
    p = VECTOR_ELT(tracts, k);  /* get the next R transaction */
    n = (ITEM)length(p);        /* and its length */
    if (t == INTSXP) {          /* if items are integer numbers */
      for (i = 0; i < n; i++) { /* traverse the items (integers) */
        if (ib_add2ta(ibase, INTEGER(p) +i) < 0) {
          e = -1; break; }      /* add items to internal transaction */
      } }                       /* and check for success */
    else {                      /* if items are character strings */
      for (i = 0; i < n; i++) { /* traverse the items (strings) */
        if (ib_add2ta(ibase, CHAR(STRING_ELT(p, i))) < 0) {
          e = -1; break; }      /* add items to internal transaction */
      }                         /* and check for success */
    }                           /* afterward set transaction weight */
    if (e) break;               /* check for success */
    ib_finta(ibase, (wgts != R_NilValue) ? INTEGER(wgts)[k] : 1);
    if (tbg_addib(tabag) < 0) { e = -1; break; }
  }                             /* add the transaction to the bag */
  if (e) { tbg_delete(tabag, 1); return NULL; }
  return tabag;                 /* return the created transaction bag */
}  /* tbg_fromRObj() */

/*--------------------------------------------------------------------*/

static int isr_Rborder (ISREPORT *rep, SEXP border)
{                               /* --- set reporter filtering border */
  ITEM   n;                     /* loop variables, border length */
  RSUPP  s;                     /* minimum support */
  int    *b;                    /* elements of the border array */
  double *d;                    /* elements of the border array */

  assert(rep && border);        /* check the function arguments */
  if (border == R_NilValue) return 0;
  n = (ITEM)length(border);     /* if no border is given or */
  if (n <= 0) return 0;         /* if it has length 0, abort */
  if      (TYPEOF(border) == INTSXP) {
    b = INTEGER(border);        /* if integer border */
    while (--n >= 0) {          /* traverse the border elements */
      s = (b[n] < 0)    ? RSUPP_MAX : (RSUPP)b[n];
      if (isr_setbdr(rep, (ITEM)n, s) < 0) return -1;
    } }                         /* set the support values */
  else if (TYPEOF(border) == REALSXP) {
    d = REAL(border);           /* if real-valued border */
    while (--n >= 0) {          /* traverse the border elements */
      s = (IS_NA(d[n])) ? RSUPP_MAX : (RSUPP)d[n];
      if (isr_setbdr(rep, (ITEM)n, s) < 0) return -1;
    }                           /* set the support values */
  }                             /* in the item set reporter */
  return 0;                     /* return 'ok' */
}  /* isr_Rborder() */

/*--------------------------------------------------------------------*/

static void isr_iset2RObj (ISREPORT *rep, void *data)
{                               /* --- report an item set */
  REPDATA *rd = data;           /* type the data pointer */
  size_t  k, n;                 /* size of result buffer, loop var. */
  ITEM    i, m;                 /* loop variable for items */
  int     v;                    /* loop variable for info. values */
  SEXP    p;                    /* resized result buffer, item buffer */
  SEXP    rset;                 /* new R object for item set */
  int     *iset;                /* item set elements (if integer) */
  SEXP    info;                 /* information for item set */
  double  *r;                   /* real array for information */
  RSUPP   supp, base;           /* item set support and base support */
  SEXP    relt;                 /* result element */

  assert(rep && data);          /* check the function arguments */
  if (rd->err) return;          /* if there was an error, do nothing */
  n = rd->size;                 /* get the current array size */
  if (rd->cnt >= n) {           /* if the result array is full */
    n += (n > BLKSIZE) ? n >> 1 : BLKSIZE;
    p = PROTECT(allocVector(VECSXP, (R_xlen_t)n));
    if (rd->res) {              /* if there is an old vector/list */
      for (k = 0; k < rd->cnt; k++)
        SET_VECTOR_ELT(p, (R_xlen_t)k, VECTOR_ELT(rd->res,(R_xlen_t)k));
      UNPROTECT(2); PROTECT(p); /* copy the existing item sets and */
    }                           /* transfer protection to new array */
    rd->res = p; rd->size = n;  /* set the (new) array/list */
  }                             /* and the new array size */
  m = isr_cnt(rep);             /* get the number of items */
  if (rd->istr) {               /* if items are strings */
    rset = PROTECT(allocVector(STRSXP, (R_xlen_t)m));
    for (i = 0; i < m; i++) {   /* create object for item set */
      p = mkChar(isr_basename(rep, isr_itemx(rep, i)));
      SET_STRING_ELT(rset,i,p); /* map identifiers to items and */
    } }                         /* store items in string array */
  else {                        /* if items are integer numbers */
    rset = PROTECT(allocVector(INTSXP, (R_xlen_t)m));
    iset = INTEGER(rset);       /* create object for item set */
    for (i = 0; i < m; i++)     /* map identifiers to items */
      iset[i] = (int)(ptrdiff_t)isr_itemobj(rep, isr_itemx(rep, i));
  }                             /* store items in integer array */
  info = PROTECT(allocVector(REALSXP, (R_xlen_t)rd->len));
  r    = REAL(info);            /* create an information array */
  supp = isr_supp(rep);         /* get the item set support */
  base = isr_suppx(rep, 0);     /* and the total transaction weight */
  for (v = 0; v < rd->len; v++){/* traverse the values to store */
    switch (rd->rep[v]) {       /* evaluate the value indicator */
      case 'a': r[v] = (double)supp;                    break;
      case 's': r[v] = (double)supp /(double)base;      break;
      case 'S': r[v] = (double)supp /(double)base *100; break;
      case 'p': r[v] = isr_eval(rep);                   break;
      case 'P': r[v] = isr_eval(rep) *100;              break;
      case 'e': r[v] = isr_eval(rep);                   break;
      case 'E': r[v] = isr_eval(rep) *100;              break;
      case 'Q': r[v] = (double)base;                    break;
      default : r[v] = 0.0;                             break;
    }                           /* get the requested value and */
  }                             /* store it in the information array */
  relt = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(relt,0,rset);  /* build element for each item set */
  SET_VECTOR_ELT(relt,1,info);  /* and store set and information */
  SET_VECTOR_ELT(rd->res, (R_xlen_t)rd->cnt, relt);
  rd->cnt += 1;                 /* store and count created item set */
  UNPROTECT(3);                 /* release the sub-objects */
}  /* isr_iset2RObj() */

/*--------------------------------------------------------------------*/

static double lift (RSUPP supp, RSUPP body, RSUPP head, RSUPP base)
{                               /* --- compute lift value of a rule */
  return ((body <= 0) || (head <= 0)) ? 0
       : ((double)supp*(double)base) /((double)body*(double)head);
}  /* lift() */

/*--------------------------------------------------------------------*/

static void isr_rule2RObj (ISREPORT *rep, void *data,
                           ITEM item, RSUPP body, RSUPP head)
{                               /* --- report an association rule */
  REPDATA *rd = data;           /* type the data pointer */
  ITEM    i, m, o, z;           /* loop variable, array size */
  int     v;                    /* loop variable, index offset */
  size_t  k, n;                 /* size of result buffer */
  SEXP    p;                    /* resized result array */
  SEXP    cons;                 /* new R object for rule head */
  SEXP    ante;                 /* new R object for rule body */
  int     *iset;                /* item set elements (if integer) */
  SEXP    info;                 /* information for item set */
  double  *r;                   /* real array for information */
  RSUPP   supp, base;           /* item set support and base support */
  SEXP    relt;                 /* result element */

  assert(rep && data            /* check the function arguments */
  &&    (body > 0) && (head > 0));
  assert(isr_uses(rep, item));  /* head item must be in item set */
  if (rd->err) return;          /* if there was an error, do nothing */
  n = rd->size;                 /* get the current array size */
  if (rd->cnt >= n) {           /* if the result array is full */
    n += (n > BLKSIZE) ? n >> 1 : BLKSIZE;
    p = PROTECT(allocVector(VECSXP, (R_xlen_t)n));
    if (rd->res) {              /* if there is an old vector/list */
      for (k = 0; k < rd->cnt; k++)
        SET_VECTOR_ELT(p, (R_xlen_t)k, VECTOR_ELT(rd->res,(R_xlen_t)k));
      UNPROTECT(2); PROTECT(p); /* copy the existing item sets and */
    }                           /* transfer protection to new array */
    rd->res = p; rd->size = n;  /* set (new) array/list and its size */
  }
  m = isr_cnt(rep);             /* get the number or items */
  if (rd->istr) {               /* if items are strings */
    cons = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(cons, 0, mkChar(isr_basename(rep, item)));
    ante = PROTECT(allocVector(STRSXP, (R_xlen_t)(m-1)));
    for (i = o = 0; i < m; i++){/* map identifiers to items */
      z = isr_itemx(rep, i);    /* get the next item and skip it */
      if (z == item) continue;  /* if it is the head of the rule */
      SET_STRING_ELT(ante, o, mkChar(isr_basename(rep, z))); o++;
    } }                         /* store items in string array */
  else {                        /* if items are integer numbers */
    cons = PROTECT(allocVector(INTSXP, 1));
    INTEGER(cons)[0] = (int)(ptrdiff_t)isr_itemobj(rep, item);
    ante = PROTECT(allocVector(INTSXP, (R_xlen_t)(m-1)));
    iset = INTEGER(ante);         /* create object for antecedent */
    for (i = o = 0; i < m; i++) { /* map identifiers to items */
      z = isr_itemx(rep, i);      /* get the next item and skip it */
      if (z == item) continue;    /* if it is the head of the rule */
      iset[o++] = (int)(ptrdiff_t)isr_itemobj(rep, z);
    }                             /* store items in integer array */
  }
  info = PROTECT(allocVector(REALSXP, (R_xlen_t)rd->len));
  r    = REAL(info);            /* create an information array */
  supp = isr_supp(rep);         /* get the item set support */
  base = isr_suppx(rep, 0);     /* and the total transaction weight */
  for (v = 0; v < rd->len; v++){/* traverse the values to store */
    switch (rd->rep[v]) {       /* evaluate the value indicator */
      case 'a': r[v] = (double)supp;                      break;
      case 'b': r[v] = (double)body;                      break;
      case 'h': r[v] = (double)head;                      break;
      case 's': r[v] = (double)supp /(double)base;        break;
      case 'S': r[v] = (double)supp /(double)base *100;   break;
      case 'x': r[v] = (double)body /(double)base;        break;
      case 'X': r[v] = (double)body /(double)base *100;   break;
      case 'y': r[v] = (double)head /(double)base;        break;
      case 'Y': r[v] = (double)head /(double)base *100;   break;
      case 'c': r[v] = (double)supp /(double)body;        break;
      case 'C': r[v] = (double)supp /(double)body *100;   break;
      case 'l': r[v] = lift(supp, body, head, base);      break;
      case 'L': r[v] = lift(supp, body, head, base) *100; break;
      case 'e': r[v] = isr_eval(rep);                     break;
      case 'E': r[v] = isr_eval(rep) *100;                break;
      case 'Q': r[v] = (double)base;                      break;
      default : r[v] = 0.0;                               break;
    }                           /* get the requested value and */
  }                             /* store it in the information array */
  relt = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(relt,0,cons);  /* build element for rule */
  SET_VECTOR_ELT(relt,1,ante);  /* and store head and body */
  SET_VECTOR_ELT(relt,2,info);  /* and rule information */
  SET_VECTOR_ELT(rd->res, (R_xlen_t)rd->cnt, relt);
  rd->cnt += 1;                 /* store and count created ass. rule */
  UNPROTECT(4);                 /* release the sub-objects */
}  /* isr_rule2RObj() */

/*--------------------------------------------------------------------*/

static SEXP psp_toRObj (PATSPEC *psp, double scale, int format)
{                               /* --- report pattern spectrum */
  ITEM   z;                     /* loop variable for sizes */
  SUPP   s;                     /* loop variable for support values */
  size_t i, n;                  /* list size, result list index */
  size_t k;                     /* number of occurrences */
  double *r;                    /* to store the values */
  SEXP   res;                   /* created R list object */
  SEXP   row;                   /* R object  for row-wise    rep. */
  SEXP   col[3];                /* R objects for column-wise rep. */
  int    *sizes;                /* pattern sizes */
  int    *supps;                /* support values */
  double *freqs;                /* occurrence frequencies */

  assert(psp);                  /* check the function arguments */
  n = psp_sigcnt(psp);          /* get the number of signatures */
  if ((format == '=')           /* if list/vector of triplets */
  ||  (format == '-')) {        /* (size, support, frequency) */
    res = PROTECT(allocVector(VECSXP, (R_xlen_t)n));
    for (i = 0, z = psp_min(psp); z <= psp_max(psp); z++) {
      for (s = psp_min4sz(psp, z); s <= psp_max4sz(psp, z); s++) {
        k = psp_getfrq(psp,z,s);/* traverse the pattern signatures */
        if (k <= 0) continue;   /* skip empty counters */
        row  = PROTECT(allocVector(REALSXP, 3));
        r    = REAL(row);       /* create a row with three elements */
        r[0] = (double)z;       /* store the patter size, */
        r[1] = (double)s;       /* the pattern support, */
        r[2] = (double)k *scale;/* and the frequency */
        SET_VECTOR_ELT(res, (R_xlen_t)i, row); i += 1;
        UNPROTECT(1);           /* set the row in the output vector */
      }                         /* and count and unprotect the row */
    } }                         /* (protected as part of 'res') */
  else {                        /* if three arrays (size,supp,freq) */
    res    = PROTECT(allocVector(VECSXP,  (R_xlen_t)3));
    col[0] = PROTECT(allocVector(INTSXP,  (R_xlen_t)n));
    SET_VECTOR_ELT(res, 0, col[0]); UNPROTECT(1);
    col[1] = PROTECT(allocVector(INTSXP,  (R_xlen_t)n));
    SET_VECTOR_ELT(res, 1, col[1]); UNPROTECT(1);
    col[2] = PROTECT(allocVector(REALSXP, (R_xlen_t)n));
    SET_VECTOR_ELT(res, 2, col[2]); UNPROTECT(1);
    sizes = INTEGER(col[0]);    /* create two integer and */
    supps = INTEGER(col[1]);    /* one real-valued array for */
    freqs = REAL   (col[2]);    /* size, support and frequency */
    for (i = 0, z = psp_min(psp); z <= psp_max(psp); z++) {
      for (s = psp_min4sz(psp, z); s <= psp_max4sz(psp, z); s++) {
        k = psp_getfrq(psp,z,s);/* traverse the pattern signatures */
        if (k <= 0) continue;   /* skip empty counters */
        sizes[i]   = z; supps[i] = s;
        freqs[i++] = (double)k *scale;
      }                         /* store the signature elements */
    }                           /* in the corresponding arrays */
  }
  assert(i == n);               /* check signature counter */
  return res;                   /* return created pattern spectrum */
}  /* psp_toRObj() */

/*--------------------------------------------------------------------*/

static int repinit (REPDATA *data, ISREPORT *isrep, CCHAR *report,
                    int target)
{                               /* --- initialize reporting */
  assert(data && isrep && report); /* check the function arguments */
  data->err = 0;                /* initialize the error indicator */
  if ((report[0] == '#')        /* if to get a pattern spectrum */
  ||  (report[0] == '|')        /* "#", "|": column-wise */
  ||  (report[0] == '=')        /* "=", "-": row-wise */
  ||  (report[0] == '-'))
    return isr_addpsp(isrep, NULL);
  data->res  = NULL;            /* initialize the report structure */
  data->size = data->cnt = 0;   /* and the array parameters */
  data->istr = ((ib_mode(isr_base(isrep)) & IB_OBJNAMES) == 0);
  data->len  = (int)strlen(data->rep = report);
  if (target & ISR_RULES) isr_setrule(isrep, isr_rule2RObj, data);
  else                    isr_setrepo(isrep, isr_iset2RObj, data);
  return 0;                     /* set the report function */
}  /* repinit() */              /* and return 'ok' */

/*--------------------------------------------------------------------*/

static int repterm (REPDATA *data, ISREPORT *isrep, CCHAR *report)
{                               /* --- terminate reporting */
  size_t k;                     /* loop variable */
  SEXP   p;                     /* resized R vector/list */

  assert(data && isrep && report); /* check the function arguments */
  if (data->err) return data->err -1;
  if ((report[0] == '#')        /* if to get a pattern spectrum */
  ||  (report[0] == '|')        /* "#", "|": column-wise */
  ||  (report[0] == '=')        /* "=", "-": row-wise */
  ||  (report[0] == '-')) {
    data->res = psp_toRObj(isr_getpsp(isrep), 1.0, report[0]);
    return data->err = (data->res) ? 0 : -1;
  }                             /* make R pattern spectrum */
  if (data->cnt != data->size) {/* if result list has wrong size */
    p = PROTECT(allocVector(VECSXP, (R_xlen_t)data->cnt));
    for (k = 0; k < data->cnt; k++) /* create a new result list/array */
      SET_VECTOR_ELT(p, (R_xlen_t)k, VECTOR_ELT(data->res,(R_xlen_t)k));
    UNPROTECT(2); PROTECT(p);   /* copy the existing item sets */
    data->res = p; data->size = data->cnt;
  }                             /* set array/list and its size */
  return data->err;             /* return the error status */
}  /* repterm() */

/*--------------------------------------------------------------------*/

static void repfn (long int cnt, void *data)
{                               /* --- progress reporting function */
  if ((cnt > *(long int*)data) && ((cnt % 20) == 0))
    fprintf(stderr, "%10ld\b\b\b\b\b\b\b\b\b\b",
                    *(long int*)data = cnt);
}  /* repfn() */

/*--------------------------------------------------------------------*/
/* fim (tracts, wgts=NULL, target="s", supp=10.0, zmin=0, zmax=-1,    */
/*      report="a", eval="x", agg="x", thresh=10.0, border=NULL)      */
/*--------------------------------------------------------------------*/

SEXP f4r_fim (SEXP ptracts, SEXP pwgts, SEXP ptarget, SEXP psupp,
              SEXP pzmin, SEXP pzmax, SEXP preport, SEXP peval,
              SEXP pagg, SEXP pthresh, SEXP pborder)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = FPG_SIMPLE;   /* algorithm variant */
  int      mode    = FPG_DEFAULT;  /* operation mode/flags */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  FPGROWTH *fpgrowth;           /* fpgrowth miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "ascmg");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_evalx(peval);    /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  agg  = get_agg(pagg);         /* get aggregation mode */
  if (agg    < 0)    error("invalid 'agg' argument");
  thresh = get_dbl(pthresh, thresh);
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  fpgrowth = fpg_create(target, supp, 100.0, 100.0,
                        (ITEM)zmin, (ITEM)zmax,
                        eval, agg, thresh, algo, mode);
  if (!fpgrowth) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = fpg_data(fpgrowth, tabag, 0, +2);
  if (r) fpg_delete(fpgrowth,1);/* prepare data for fpgrowth */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (fpg_report(fpgrowth, isrep)           != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    fpg_delete(fpgrowth, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = fpg_mine(fpgrowth, ITEM_MIN, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  fpg_delete(fpgrowth, 1);      /* delete the fpgrowth miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_fim() */              /* return the created R object */

/*--------------------------------------------------------------------*/
/* arules (tracts, wgts=NULL, supp=10.0, conf=80.0,                   */
/*         zmin=0, zmax=-1, report="aC", eval="x", thresh=10.0,       */
/*         appear=NULL)                                               */
/*--------------------------------------------------------------------*/

SEXP f4r_arules (SEXP ptracts, SEXP pwgts, SEXP psupp,   SEXP pconf,
                 SEXP pzmin,   SEXP pzmax, SEXP preport, SEXP peval,
                 SEXP pthresh, SEXP pmode, SEXP pappear)
{                               /* --- association rule induction */
  double   supp    = 10;        /* minimum support of a rule */
  double   conf    = 80;        /* minimum confidence of a rule */
  int      zmin    =  0;        /* minimum size of a rule */
  int      zmax    = -1;        /* maximum size of a rule */
  CCHAR    *report = "aC";      /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = FPG_SIMPLE;   /* algorithm variant */
  int      mode    = FPG_DEFAULT;  /* operation mode/flags */
  CCHAR    *smode  = "";        /* operation mode as a string */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  FPGROWTH *fpgrowth;           /* fpgrowth miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, pappear);
  if (r < -3) error("invalid 'appear' argument");
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  supp   = get_dbl(psupp, supp);
  conf   = get_dbl(pconf, conf);
  if (conf   < 0)    error("invalid 'conf' argument (must be >= 0)");
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_evalx(peval);    /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  thresh = get_dbl(pthresh, thresh);
  smode  = get_str(pmode, smode);

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, pappear);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  fpgrowth = fpg_create(FPG_RULES, supp, 100.0, conf,
                        (ITEM)zmin, (ITEM)zmax,
                        eval, FPG_NONE, thresh, algo, mode);
  if (!fpgrowth) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = fpg_data(fpgrowth, tabag, 0, +2);
  if (r) fpg_delete(fpgrowth,1);/* prepare data for fpgrowth */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (fpg_report(fpgrowth, isrep)              != 0)
  ||  (repinit(&data, isrep, report, ISR_RULES) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    fpg_delete(fpgrowth, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = fpg_mine(fpgrowth, ITEM_MIN, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  fpg_delete(fpgrowth, 1);      /* delete the fpgrowth miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_arules() */           /* return the created R object */

/*--------------------------------------------------------------------*/
/* apriori (tracts, wgts=NULL, target="s", supp=10, zmin=0, zmax=-1,  */
/*          report="a", eval="x", agg="x", thresh=10.0, prune=NA,     */
/*          algo="a", mode="", border=NULL, appear=NULL)              */
/*--------------------------------------------------------------------*/

SEXP f4r_apriori (SEXP ptracts, SEXP pwgts, SEXP ptarget,
                  SEXP psupp, SEXP pconf, SEXP pzmin, SEXP pzmax,
                  SEXP preport, SEXP peval, SEXP pagg, SEXP pthresh,
                  SEXP pprune,  SEXP palgo, SEXP pmode,
                  SEXP pborder, SEXP pappear)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  double   conf    = 80;        /* minimum confidence of a rule */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = APR_AUTO;     /* algorithm variant */
  int      mode    = APR_DEFAULT;  /* operation mode/flags */
  CCHAR    *s      = "";        /* to access the operation mode/flags */
  long int prune   = LONG_MIN;  /* min. size for evaluation filtering */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  APRIORI  *apriori;            /* apriori miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, pappear);
  if (r < -3) error("invalid 'appear' argument "
                    "(must be pair of item and string array)");
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "ascmgr");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  conf   = get_dbl(pconf, conf);
  if (conf   < 0)    error("invalid 'conf' argument (must be >= 0)");
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_evalx(peval);    /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  agg    = get_agg(pagg);       /* get aggregation mode */
  if (agg    < 0)    error("invalid 'agg' argument");
  thresh = get_dbl(pthresh, thresh);
  prune  = get_lng(pprune,  prune);
  if (eval <= RE_NONE) prune = LONG_MIN;
  if (TYPEOF(palgo) != STRSXP)
    error("invalid Apriori algorithm variant");
  s = (length(palgo) > 0) ? CHAR(STRING_ELT(palgo, 0)) : "auto";
  if (s[0] && s[1]) {           /* if textual identifier */
    if      (strcmp(s, "auto")   == 0) s    = "a";
    else if (strcmp(s, "basic")  == 0) s    = "b";
    else                               algo = -1;
  }                             /* translate the algorithm string */
  if (s[0] && !s[1]) {          /* if single character, */
    switch (s[0]) {             /* evaluate the algorithm code */
      case 'a': algo = APR_AUTO;  break;
      case 'b': algo = APR_BASIC; break;
      default : algo = -1;        break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) error("invalid Apriori algorithm");
  if (TYPEOF(pmode) != STRSXP) error("invalid 'mode' argument");
  s = (length(pmode) > 0) ? CHAR(STRING_ELT(pmode, 0)) : "";
  for ( ; *s; s++) {            /* traverse the mode characters */
    if      (*s == 'o') mode |=  APR_ORIGSUPP;
    else if (*s == 'x') mode &= ~APR_PERFECT;
    else if (*s == 't') mode &= ~APR_TATREE;
    else if (*s == 'T') mode &= ~APR_TATREE;
    else if (*s == 'y') mode |=  APR_POST;
    else if (*s == 'z') eval |=  IST_INVBXS;
  }                             /* adapt the operation mode */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- get and prepare transactions --- */
  sig_install();                /* install the signal handler */
  if (!(target & ISR_RULES)) pappear = R_NilValue;
  tabag = tbg_fromRObj(ptracts, pwgts, pappear);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  apriori = apriori_create(target, supp, 100.0, conf,
                           (ITEM)zmin, (ITEM)zmax,
                           eval, agg, thresh, algo, mode);
  if (!apriori) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = apriori_data(apriori, tabag, 0, +2);
  if (r) apriori_delete(apriori, 1);   /* prepare data for apriori */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (apriori_report(apriori, isrep)        != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    apriori_delete(apriori, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = apriori_mine(apriori, (ITEM)prune, 0.01, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  apriori_delete(apriori, 1);   /* delete the apriori miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  if (r != 0) ERR_MEM();        /* check for an error */
  sig_remove();                 /* remove the signal handler */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_apriori() */          /* return the created R object */

/*--------------------------------------------------------------------*/
/* eclat (tracts, wgts=NULL, target="s", supp=10, zmin=0, zmax=-1,    */
/*        report="a", eval="x", agg="x", thresh=10.0, prune=NA,       */
/*        algo="a", mode="", border=NULL, appear=NULL)                */
/*--------------------------------------------------------------------*/

SEXP f4r_eclat (SEXP ptracts, SEXP pwgts, SEXP ptarget,
                SEXP psupp, SEXP pconf, SEXP pzmin, SEXP pzmax,
                SEXP preport, SEXP peval, SEXP pagg, SEXP pthresh,
                SEXP pprune,  SEXP palgo, SEXP pmode,
                SEXP pborder, SEXP pappear)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  double   conf    = 80;        /* minimum confidence of a rule */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = ECL_OCCDLV;   /* algorithm variant */
  int      mode    = ECL_DEFAULT;  /* operation mode/flags */
  CCHAR    *s      = "";        /* to access the operation mode/flags */
  long int prune   = LONG_MIN;  /* min. size for evaluation filtering */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  ECLAT    *eclat;              /* eclat miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, pappear);
  if (r < -3) error("invalid 'appear' argument "
                    "(must be pair of item and string array)");
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "ascmgr");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  conf   = get_dbl(pconf, conf);
  if (conf   < 0)    error("invalid 'conf' argument (must be >= 0)");
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_evalx(peval);    /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  agg    = get_agg(pagg);       /* get aggregation mode */
  if (agg    < 0)    error("invalid 'agg' argument");
  thresh = get_dbl(pthresh, thresh);
  prune  = get_lng(pprune,  prune);
  if (eval <= RE_NONE) prune = LONG_MIN;
  if (TYPEOF(palgo) != STRSXP) error("invalid Eclat algorithm");
  s = (length(palgo) > 0) ? CHAR(STRING_ELT(palgo, 0)) : "auto";
  if (s[0] && s[1]) {           /* if textual identifier */
    if      (strcmp(s, "auto")   == 0) s    = "a";
    else if (strcmp(s, "basic")  == 0) s    = "e";
    else if (strcmp(s, "lists")  == 0) s    = "i";
    else if (strcmp(s, "tids")   == 0) s    = "i";
    else if (strcmp(s, "bits")   == 0) s    = "b";
    else if (strcmp(s, "table")  == 0) s    = "t";
    else if (strcmp(s, "simple") == 0) s    = "s";
    else if (strcmp(s, "ranges") == 0) s    = "r";
    else if (strcmp(s, "occdlv") == 0) s    = "o";
    else if (strcmp(s, "diff")   == 0) s    = "d";
    else                               algo = -1;
  }                             /* translate the algorithm string */
  if (s[0] && !s[1]) {          /* if single character, */
    switch (s[0]) {             /* evaluate the algorithm code */
      case 'a': algo = ECL_AUTO;   break;
      case 'e': algo = ECL_BASIC;  break;
      case 'i': algo = ECL_LISTS;  break;
      case 'b': algo = ECL_BITS;   break;
      case 't': algo = ECL_TABLE;  break;
      case 's': algo = ECL_SIMPLE; break;
      case 'r': algo = ECL_RANGES; break;
      case 'o': algo = ECL_OCCDLV; break;
      case 'd': algo = ECL_DIFFS;  break;
      default : algo = -1;         break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) error("invalid Eclat algorithm");
  if (TYPEOF(pmode) != STRSXP) error("invalid 'mode' argument");
  s = (length(pmode) > 0) ? CHAR(STRING_ELT(pmode, 0)) : "";
  for ( ; *s; s++) {            /* traverse the mode characters */
    if      (*s == 'o') mode |=  ECL_ORIGSUPP;
    else if (*s == 'l') mode &= ~ECL_FIM16;
    else if (*s == 'x') mode &= ~ECL_PERFECT;
    else if (*s == 'i') mode &= ~ECL_REORDER;
    else if (*s == 'u') mode &= ~ECL_TAIL;
    else if (*s == 'y') mode |=  ECL_HORZ;
    else if (*s == 'Y') mode |=  ECL_VERT;
    else if (*s == 'z') eval |=  IST_INVBXS;
  }                             /* adapt the operation mode */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- get and prepare transactions --- */
  sig_install();                /* install the signal handler */
  if (!(target & ISR_RULES)) pappear = R_NilValue;
  tabag = tbg_fromRObj(ptracts, pwgts, pappear);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  eclat = eclat_create(target, supp, 100.0, conf,
                       (ITEM)zmin, (ITEM)zmax,
                       eval, agg, thresh, algo, mode);
  if (!eclat) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = eclat_data(eclat, tabag, 0, +2);
  if (r) eclat_delete(eclat, 1);/* prepare data for eclat */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (eclat_report(eclat, isrep)            != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    eclat_delete(eclat, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = eclat_mine(eclat, (ITEM)prune, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  eclat_delete(eclat, 1);       /* delete the eclat miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_eclat() */            /* return the created R object */

/*--------------------------------------------------------------------*/
/* fpgrowth (tracts, wgts=NULL, target="s", supp=10, zmin=0, zmax=-1, */
/*           report="a", eval="x", agg="x", thresh=10.0, prune=NA,    */
/*           algo="a", mode="", border=NULL, appear=NULL)             */
/*--------------------------------------------------------------------*/

SEXP f4r_fpgrowth (SEXP ptracts, SEXP pwgts, SEXP ptarget,
                   SEXP psupp, SEXP pconf, SEXP pzmin, SEXP pzmax,
                   SEXP preport, SEXP peval, SEXP pagg, SEXP pthresh,
                   SEXP pprune,  SEXP palgo, SEXP pmode,
                   SEXP pborder, SEXP pappear)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  double   conf    = 80;        /* minimum confidence of a rule */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  int      agg     =  0;        /* aggregation mode */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = FPG_SIMPLE;   /* algorithm variant */
  int      mode    = FPG_DEFAULT;  /* operation mode/flags */
  CCHAR    *s      = "";        /* to access the operation mode/flags */
  long int prune   = LONG_MIN;  /* min. size for evaluation filtering */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  FPGROWTH *fpgrowth;           /* fpgrowth miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, pappear);
  if (r < -3) error("invalid 'appear' argument "
                    "(must be pair of item and string array)");
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "ascmgr");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  conf   = get_dbl(pconf, conf);
  if (conf   < 0)    error("invalid 'conf' argument (must be >= 0)");
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_evalx(peval);    /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  agg    = get_agg(pagg);       /* get aggregation mode */
  if (agg    < 0)    error("invalid 'agg' argument");
  thresh = get_dbl(pthresh, thresh);
  prune  = get_lng(pprune,  prune);
  if (eval <= RE_NONE) prune = LONG_MIN;
  if (TYPEOF(palgo) != STRSXP) error("invalid FP-growth algorithm");
  s = (length(palgo) > 0) ? CHAR(STRING_ELT(palgo, 0)) : "auto";
  if (s[0] && s[1]) {           /* if textual identifier */
    if      (strcmp(s, "auto")    == 0) s    = "a";
    else if (strcmp(s, "simple")  == 0) s    = "s";
    else if (strcmp(s, "complex") == 0) s    = "c";
    else if (strcmp(s, "single")  == 0) s    = "d";
    else if (strcmp(s, "topdown") == 0) s    = "t";
    else                                algo = -1;
  }                             /* translate the algorithm string */
  if (s[0] && !s[1]) {          /* if single character, */
    switch (s[0]) {             /* evaluate the algorithm code */
      case 'a': algo = FPG_SIMPLE;  break;
      case 's': algo = FPG_SIMPLE;  break;
      case 'c': algo = FPG_COMPLEX; break;
      case 'd': algo = FPG_SINGLE;  break;
      case 't': algo = FPG_TOPDOWN; break;
      default : algo = -1;          break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) error("invalid FP-growth algorithm");
  if (TYPEOF(pmode) != STRSXP) error("invalid 'mode' argument");
  s = (length(pmode) > 0) ? CHAR(STRING_ELT(pmode, 0)) : "";
  for ( ; *s; s++) {            /* traverse the mode characters */
    if      (*s == 'o') mode |=  FPG_ORIGSUPP;
    else if (*s == 'l') mode &= ~FPG_FIM16;
    else if (*s == 'x') mode &= ~FPG_PERFECT;
    else if (*s == 'i') mode &= ~FPG_REORDER;
    else if (*s == 'u') mode &= ~FPG_TAIL;
    else if (*s == 'z') eval |=  FPG_INVBXS;
  }                             /* adapt the operation mode */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  if (!(target & ISR_RULES)) pappear = R_NilValue;
  tabag = tbg_fromRObj(ptracts, pwgts, pappear);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  fpgrowth = fpg_create(target, supp, 100.0, conf,
                        (ITEM)zmin, (ITEM)zmax,
                        eval, agg, thresh, algo, mode);
  if (!fpgrowth) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = fpg_data(fpgrowth, tabag, 0, +2);
  if (r) fpg_delete(fpgrowth,1);/* prepare data for fpgrowth */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (fpg_report(fpgrowth, isrep)           != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    fpg_delete(fpgrowth, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = fpg_mine(fpgrowth, (ITEM)prune, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  fpg_delete(fpgrowth, 1);      /* delete the fpgrowth miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_fpgrowth() */         /* return the created R object */

/*--------------------------------------------------------------------*/
/* sam (tracts, wgts=NULL, target="s", supp=10, zmin=0, zmax=-1,      */
/*      report="a", eval="x", thresh=10.0, algo="a", mode="",         */
/*      border=NULL)                                                  */
/*--------------------------------------------------------------------*/

SEXP f4r_sam (SEXP ptracts, SEXP pwgts, SEXP ptarget, SEXP psupp,
              SEXP pzmin,   SEXP pzmax, SEXP preport, SEXP peval,
              SEXP pthresh, SEXP palgo, SEXP pmode, SEXP pborder)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = SAM_BSEARCH;  /* algorithm variant */
  int      mode    = SAM_DEFAULT;  /* operation mode/flags */
  CCHAR    *s      = "";        /* to access the operation mode/flags */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  SAM      *sam;                /* split and merge miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "ascmgr");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_eval(peval);     /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  thresh = get_dbl(pthresh, thresh);
  if (TYPEOF(palgo) != STRSXP) error("invalid SaM algorithm");
  s = (length(palgo) > 0) ? CHAR(STRING_ELT(palgo, 0)) : "auto";
  if (s[0] && s[1]) {           /* if textual identifier */
    if      (strcmp(s, "auto")    == 0) s    = "a";
    else if (strcmp(s, "basic")   == 0) s    = "s";
    else if (strcmp(s, "simple")  == 0) s    = "s";
    else if (strcmp(s, "bsearch") == 0) s    = "b";
    else if (strcmp(s, "double")  == 0) s    = "d";
    else if (strcmp(s, "tree")    == 0) s    = "t";
    else                                algo = -1;
  }                             /* translate the algorithm string */
  if (s[0] && !s[1]) {          /* if single character, */
    switch (s[0]) {             /* evaluate the algorithm code */
      case 'a': algo = SAM_BSEARCH; break;
      case 's': algo = SAM_BASIC;   break;
      case 'b': algo = SAM_BSEARCH; break;
      case 'd': algo = SAM_DOUBLE;  break;
      case 't': algo = SAM_TREE;    break;
      default : algo = -1;          break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) error("invalid SaM algorithm");
  if (TYPEOF(pmode) != STRSXP) error("invalid 'mode' argument");
  s = (length(pmode) > 0) ? CHAR(STRING_ELT(pmode, 0)) : "";
  for ( ; *s; s++) {            /* traverse the mode characters */
    if      (*s == 'l') mode &= ~SAM_FIM16;
    else if (*s == 'x') mode &= ~SAM_PERFECT;
  }                             /* adapt the operation mode */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  sam = sam_create(target, supp, 0.0, (ITEM)zmin, (ITEM)zmax,
                   0, -1.0, eval, thresh, algo, mode);
  if (!sam) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = sam_data(sam, tabag, +2); /* create a split and merge miner */
  if (r) sam_delete(sam, 1);    /* prepare data for split and merge */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (sam_report(sam, isrep)                != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    sam_delete(sam, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = sam_mine(sam, 8192);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  sam_delete(sam, 1);           /* delete the split and merge miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_sam() */              /* return the created R object */

/*--------------------------------------------------------------------*/
/* relim (tracts, wgts=NULL, target="s", supp=10, zmin=0, zmax=-1,    */
/*        report="a", eval="x", thresh=10.0, algo="a", mode="",       */
/*        border=NULL)                                                */
/*--------------------------------------------------------------------*/

SEXP f4r_relim (SEXP ptracts, SEXP pwgts, SEXP ptarget, SEXP psupp,
                SEXP pzmin,   SEXP pzmax, SEXP preport, SEXP peval,
                SEXP pthresh, SEXP palgo, SEXP pmode, SEXP pborder)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = REL_BASIC;              /* algorithm variant */
  int      mode    = REL_DEFAULT|REL_FIM16;  /* operation mode/flags */
  CCHAR    *s      = "";        /* to access the operation mode/flags */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  RELIM    *relim;              /* relim miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "ascmgr");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_eval(peval);     /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  thresh = get_dbl(pthresh, thresh);
  if (TYPEOF(palgo) != STRSXP) error("invalid RElim algorithm");
  s = (length(palgo) > 0) ? CHAR(STRING_ELT(palgo, 0)) : "auto";
  if (s[0] && s[1]) {           /* if textual identifier */
    if      (strcmp(s, "auto")    == 0) s    = "a";
    else if (strcmp(s, "basic")   == 0) s    = "s";
    else if (strcmp(s, "simple")  == 0) s    = "s";
    else                                algo = -1;
  }                             /* translate the algorithm string */
  if (s[0] && !s[1]) {          /* if single character, */
    switch (s[0]) {             /* evaluate the algorithm code */
      case 'a': algo = REL_BASIC; break;
      case 's': algo = REL_BASIC; break;
      default : algo = -1;        break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) error("invalid RElim algorithm");
  if (TYPEOF(pmode) != STRSXP) error("invalid 'mode' argument");
  s = (length(pmode) > 0) ? CHAR(STRING_ELT(pmode, 0)) : "";
  for ( ; *s; s++) {            /* traverse the mode characters */
    if      (*s == 'l') mode &= ~REL_FIM16;
    else if (*s == 'x') mode &= ~REL_PERFECT;
  }                             /* adapt the operation mode */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  relim = relim_create(target, supp, 0.0, (ITEM)zmin, (ITEM)zmax,
                       0, -1.0, eval, thresh, algo, mode);
  if (!relim) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = relim_data(relim, tabag, +2);
  if (r) relim_delete(relim,1); /* prepare data for relim */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (relim_report(relim, isrep)            != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    relim_delete(relim, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = relim_mine(relim, 32);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  relim_delete(relim, 1);       /* delete the relim miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_relim() */            /* return the created R object */

/*--------------------------------------------------------------------*/
/* carpenter (tracts, wgts=NULL, target="c", supp=10, zmin=0, zmax=-1,*/
/*            report="a", eval="x", thresh=10.0, algo="a", mode="",   */
/*            border=NULL)                                            */
/*--------------------------------------------------------------------*/

SEXP f4r_carpenter (SEXP ptracts, SEXP pwgts, SEXP ptarget, SEXP psupp,
                    SEXP pzmin,   SEXP pzmax, SEXP preport, SEXP peval,
                    SEXP pthresh, SEXP palgo, SEXP pmode, SEXP pborder)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = CARP_AUTO;    /* algorithm variant */
  int      mode    = CARP_DEFAULT; /* operation mode/flags */
  CCHAR    *s      = "";        /* to access the operation mode/flags */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  CARP     *carp;               /* carpenter miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "cm");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_eval(peval);     /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  thresh = get_dbl(pthresh, thresh);
  if (TYPEOF(palgo) != STRSXP) error("invalid Carpenter algorithm");
  s = (length(palgo) > 0) ? CHAR(STRING_ELT(palgo, 0)) : "auto";
  if (s[0] && s[1]) {           /* if textual identifier */
    if      (strcmp(s, "auto")    == 0) s    = "a";
    else if (strcmp(s, "table")   == 0) s    = "t";
    else if (strcmp(s, "table")   == 0) s    = "t";
    else if (strcmp(s, "tids")    == 0) s    = "l";
    else if (strcmp(s, "tidlist") == 0) s    = "l";
    else if (strcmp(s, "list")    == 0) s    = "l";
    else                                algo = -1;
  }                             /* translate the algorithm string */
  if (s[0] && !s[1]) {          /* if single character, */
    switch (s[0]) {             /* evaluate the algorithm code */
      case 'a': algo = CARP_AUTO;    break;
      case 't': algo = CARP_TABLE;   break;
      case 'l': algo = CARP_TIDLIST; break;
      default : algo = -1;           break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) error("invalid Carpenter algorithm");
  if (TYPEOF(pmode) != STRSXP) error("invalid 'mode' argument");
  s = (length(pmode) > 0) ? CHAR(STRING_ELT(pmode, 0)) : "";
  for ( ; *s; s++) {            /* traverse the mode characters */
    if      (*s == 'x') mode &= ~CARP_PERFECT;
    else if (*s == 'z') mode |=  CARP_FILTER;
    else if (*s == 'y') mode &= ~CARP_MAXONLY;
    else if (*s == 'p') mode &= ~CARP_COLLATE;
  }                             /* adapt the operation mode */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  carp = carp_create(target, supp, 100.0, (ITEM)zmin, (ITEM)zmax,
                     eval, thresh, algo, mode);
  if (!carp) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = carp_data(carp, tabag,-2);/* create a carpenter miner */
  if (r) carp_delete(carp, 1);  /* prepare data for carpenter */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (carp_report(carp, isrep)              != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    carp_delete(carp, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = carp_mine(carp);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  carp_delete(carp, 1);         /* delete the carpenter miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_carpenter() */        /* return the created R object */

/*--------------------------------------------------------------------*/
/* ista (tracts, wgts=NULL, target="c", supp=10, zmin=0, zmax=-1,     */
/*       report="a", eval="x", thresh=10.0, algo="a", mode="",        */
/*       border=NULL)                                                 */
/*--------------------------------------------------------------------*/

SEXP f4r_ista (SEXP ptracts, SEXP pwgts, SEXP ptarget, SEXP psupp,
               SEXP pzmin,   SEXP pzmax, SEXP preport, SEXP peval,
               SEXP pthresh, SEXP palgo, SEXP pmode, SEXP pborder)
{                               /* --- frequent item set mining */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  int      eval    = 'x';       /* evaluation measure */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = ISTA_AUTO;    /* algorithm variant */
  int      mode    = ISTA_DEFAULT; /* operation mode/flags */
  CCHAR    *s      = "";        /* to access the operation mode/flags */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  ISTA     *ista;               /* ista miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -2) error("invalid 'wgts' argument "
                    "(must be numeric array same length as 'tracts')");
  if (r < -1) error("invalid 'tracts' argument "
                    "(must be list of integer or string arrays)");
  if (r <  0) error("invalid 'tracts' argument "
                    "(must be non-empty list)");
  target = get_target(ptarget, "cm");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  eval   = get_eval(peval);     /* get evaluation measure */
  if (eval   < 0)    error("invalid 'eval' argument");
  thresh = get_dbl(pthresh, thresh);
  if (TYPEOF(palgo) != STRSXP) error("invalid IsTa algorithm");
  s = (length(palgo) > 0) ? CHAR(STRING_ELT(palgo, 0)) : "auto";
  if (s[0] && s[1]) {           /* if textual identifier */
    if      (strcmp(s, "auto")     == 0) s    = "a";
    else if (strcmp(s, "pfx")      == 0) s    = "x";
    else if (strcmp(s, "prefix")   == 0) s    = "x";
    else if (strcmp(s, "pat")      == 0) s    = "p";
    else if (strcmp(s, "patricia") == 0) s    = "p";
    else                                 algo = -1;
  }                             /* translate the algorithm string */
  if (s[0] && !s[1]) {          /* if single character, */
    switch (s[0]) {             /* evaluate the algorithm code */
      case 'a': algo = ISTA_AUTO;     break;
      case 'x': algo = ISTA_PREFIX;   break;
      case 'p': algo = ISTA_PATRICIA; break;
      default : algo = -1;            break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) error("invalid IsTa algorithm");
  if (TYPEOF(pmode) != STRSXP) error("invalid 'mode' argument");
  s = (length(pmode) > 0) ? CHAR(STRING_ELT(pmode, 0)) : "";
  for ( ; *s; s++) {            /* traverse the mode characters */
    if      (*s == 'p') mode &= ~ISTA_PRUNE;
    else if (*s == 'z') mode |=  ISTA_FILTER;
  }                             /* adapt the operation mode */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  ista = ista_create(target, supp, 100.0, (ITEM)zmin, (ITEM)zmax,
                     eval, thresh, algo, mode);
  if (!ista) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = ista_data(ista, tabag,-2);/* create an ista miner */
  if (r) ista_delete(ista, 1);  /* prepare data for ista */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (ista_report(ista, isrep)              != 0)
  ||  (isr_Rborder(isrep, pborder)           != 0)
  ||  (repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    ista_delete(ista, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = ista_mine(ista);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  ista_delete(ista, 1);         /* delete the ista miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_ista() */             /* return the created R object */

/*--------------------------------------------------------------------*/
/* apriacc (tracts, wgts=NULL, supp=-2, zmin=2, zmax=-1, report="aP", */
/*          stat="c", siglvl=1.0, prune=NA, mode="", border=NULL)     */
/*--------------------------------------------------------------------*/

SEXP f4r_apriacc (SEXP ptracts, SEXP pwgts, SEXP psupp,
                  SEXP pzmin, SEXP pzmax,   SEXP preport,
                  SEXP pstat, SEXP psiglvl, SEXP pprune,
                  SEXP pmode, SEXP pborder)
{                               /* --- frequent item set mining */
  double   supp    = -2;        /* minimum support of an item set */
  int      zmin    =  2;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aP";      /* indicators of values to report */
  int      stat    = 'c';       /* test statistic */
  double   siglvl  =  1;        /* minimum evaluation measure value */
  int      mode    = APR_DEFAULT;  /* operation mode/flags */
  long int prune   = 0;         /* min. size for evaluation filtering */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  APRIORI  *apriori;            /* apriori miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -1)        error("invalid 'wgts' argument");
  if (r <  0)        error("invalid 'tracts' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  stat   = get_stat(pstat);     /* get statistic */
  if (stat   < 0)    error("invalid 'stat' argument");
  siglvl = get_dbl(psiglvl, siglvl);
  prune  = get_lng(pprune,  prune);
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- get and prepare transactions --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  apriori = apriori_create(ISR_MAXIMAL, supp, 100.0, 100.0,
                           (ITEM)zmin, (ITEM)zmax,
                           stat, APR_MAX, siglvl, APR_AUTO, mode);
  if (!apriori) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = apriori_data(apriori, tabag, 0, +2);
  if (r) apriori_delete(apriori, 1);   /* prepare data for apriori */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (apriori_report(apriori, isrep)          != 0)
  ||  (isr_Rborder(isrep, pborder)             != 0)
  ||  (repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    apriori_delete(apriori, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = apriori_mine(apriori, (ITEM)prune, 0.01, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  apriori_delete(apriori, 1);   /* delete the apriori miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  if (r != 0) ERR_MEM();        /* check for an error */
  sig_remove();                 /* remove the signal handler */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_apriacc() */          /* return the created R object */

/*--------------------------------------------------------------------*/
/* apriacc (tracts, wgts=NULL, supp=-2, zmin=2, zmax=-1, report="aP", */
/*          stat="c", siglvl=1.0, prune=NA, mode="", border=NULL)     */
/*--------------------------------------------------------------------*/

SEXP f4r_accretion (SEXP ptracts, SEXP pwgts, SEXP psupp,
                    SEXP pzmin, SEXP pzmax,   SEXP preport,
                    SEXP pstat, SEXP psiglvl, SEXP pmaxext,
                    SEXP pmode, SEXP pborder)
{                               /* --- frequent item set mining */
  double   supp    =  1;        /* minimum support of an item set */
  int      zmin    =  2;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aP";      /* indicators of values to report */
  int      stat    = 'c';       /* test statistic */
  double   siglvl  =  1;        /* minimum evaluation measure value */
  int      mode    = ACC_DEFAULT;  /* operation mode/flags */
  long int maxext  =  2;        /* maximum number of extension items */
  TABAG    *tabag;              /* created transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  ACCRET   *accret;             /* accretion miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -1)        error("invalid 'wgts' argument");
  if (r <  0)        error("invalid 'tracts' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin, zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax, zmax);/* check the size range */
  if (zmax   < 0)    zmax = ITEM_MAX;
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  stat   = get_stat(pstat);     /* get statistic */
  if (stat   < 0)    error("invalid 'stat' argument");
  siglvl = get_dbl(psiglvl, siglvl); /* get significance level and */
  maxext = get_lng(pmaxext, maxext); /* maximum number of extensions */
  if (maxext < 0)               /* a negative value means that */
    maxext = LONG_MAX;          /* there is no limit on extensions */
  if ((pborder != R_NilValue)   /* check the filtering border */
  &&  (TYPEOF(pborder) != INTSXP) && (TYPEOF(pborder) != REALSXP))
    error("invalid 'border' argument (must be numeric)");

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create & init. transaction bag */
  accret = accret_create(ISR_MAXIMAL, supp, 100.0,
                         (ITEM)zmin, (ITEM)zmax, stat, siglvl, mode);
  if (!accret) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = accret_data(accret, tabag, +2);
  if (r) accret_delete(accret, 1);  /* prepare data for accretion */
  if (r == -1) ERR_MEM();       /* check for an error and no items */
  if (r <   0) { sig_remove(); return allocVector(VECSXP, 0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (accret_report(accret, isrep)            != 0)
  ||  (isr_Rborder(isrep, pborder)             != 0)
  ||  (repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    accret_delete(accret, 1); ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (maxext > ITEM_MAX) maxext = ITEM_MAX;
  r = accret_mine(accret, (ITEM)maxext);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  accret_delete(accret, 1);     /* delete the accretion miner */
  if (data.res) UNPROTECT(1);   /* unprotect the result object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  if (r != 0) ERR_MEM();        /* check for an error */
  return (data.res) ? data.res : allocVector(VECSXP, 0);
}  /* f4r_accretion() */        /* return the created R object */

/*--------------------------------------------------------------------*/
/* genpsp (tracts, wgts=NULL, target="s",                             */
/*         supp=10.0, zmin=0, zmax=-1, report="|",                    */
/*         cnt=1000, surr="p", seed=0, cpus=0)                        */
/*--------------------------------------------------------------------*/

SEXP f4r_genpsp (SEXP ptracts, SEXP pwgts, SEXP ptarget,
                 SEXP psupp,   SEXP pzmin, SEXP pzmax, SEXP preport,
                 SEXP pcnt,    SEXP psurr, SEXP pseed, SEXP pcpus)
{                               /* --- generate a pattern spectrum */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support of an item set */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "|";       /* indicators of reporting format */
  int      cnt     = 1000;      /* number of data sets to generate */
  int      surr    =  5;        /* surrogate method identifier */
  long int seed    =  0;        /* seed for random number generator */
  int      cpus    =  0;        /* number of cpus */
  PATSPEC  *psp    = NULL;      /* created pattern spectrum */
  SEXP     rpsp    = NULL;      /* created R pattern spectrum */
  long int done    = 0;         /* number of completed data sets */
  TABAG    *tabag;              /* transaction bag (C) */
  double   wgt;                 /* total transaction weight */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -1) error("invalid 'wgts' argument");
  if (r <  0) error("invalid 'tracts' argument");
  target = get_target(ptarget, "ascmg");
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin,  zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax,  zmax);
  if (zmax   < 0)    zmax = ITEM_MAX;  /* check the size range */
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  cnt    = get_int(pcnt, cnt);
  if (cnt   <= 0)    cnt = 1;   /* check the number of data sets */
  surr   = get_surr(psurr);     /* translate the surrogate string */
  if (surr  <  0)    error("invalid 'surr' argument");;
  if (surr  == 0)    cnt = 1;   /* only one surrogate for identity */
  seed   = get_lng(pseed, seed);
  cpus   = get_int(pcpus, cpus);

  /* --- generate pattern spectrum --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create a transaction bag */
  if ((surr == FPG_SHUFFLE) && !tbg_istab(tabag)) {
    tbg_delete(tabag, 1);       /* if shuffle surrogates requested */
    error("for shuffle surrogates transactions must form a table");
  }                             /* check for table-derived data */
  wgt  = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? supp/100.0 *(double)wgt *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  psp  = fpg_genpsp(tabag, target, (SUPP)smin, (ITEM)zmin, (ITEM)zmax,
                    FPG_SIMPLE, FPG_DEFAULT, (size_t)cnt, surr, seed,
                    cpus, repfn, &done);
  if (psp) { rpsp = psp_toRObj(psp, 1.0/(double)cnt, report[0]);
             psp_delete(psp); } /* generate a pattern spectrum */
  tbg_delete(tabag, 1);         /* delete the transaction bag */
  if (!rpsp) ERR_MEM();         /* check for an error, otherwise */
  else       UNPROTECT(1);      /* unprotect the created object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  return rpsp;                  /* return the created R object */
}  /* f4r_genpsp() */

/*--------------------------------------------------------------------*/
/* estpsp (tracts, wgts=NULL, target="a",                             */
/*         supp=10.0, zmin=0, zmax=-1, report="|",                    */
/*         equiv=10000, alpha=0.5, smpls=1000, seed=0)                */
/*--------------------------------------------------------------------*/

SEXP f4r_estpsp (SEXP ptracts, SEXP pwgts,   SEXP ptarget,
                 SEXP psupp,   SEXP pzmin,   SEXP pzmax,  SEXP preport,
                 SEXP pequiv,  SEXP palpha,  SEXP psmpls, SEXP pseed)
{                               /* --- pattern spectrum estimation */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    = 10.0;      /* minimum support of an item set */
  SUPP     smin    =  1;        /* minimum support of an item set */
  int      zmin    =  0;        /* minimum size of an item set */
  int      zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "|";       /* indicators of reporting format */
  int      equiv   = 10000;     /* equivalent number of surrogates */
  double   alpha   =  0.5;      /* probability dispersion factor */
  int      smpls   = 1000;      /* number of samples per set size */
  long int seed    =  0;        /* seed for random number generator */
  PATSPEC  *psp    = NULL;      /* created pattern spectrum */
  SEXP     rpsp    = NULL;      /* created R pattern spectrum */
  TABAG    *tabag;              /* created transaction bag */
  double   wgt;                 /* total transaction weight */
  int      r;                   /* result of function call */

  /* --- evaluate function arguments --- */
  r = chktracts(ptracts, pwgts, R_NilValue);
  if (r < -1) error("invalid 'wgts' argument");
  if (r <  0) error("invalid 'tracts' argument");
  target = get_target(ptarget, "as");  /* translate the target string */
  if (target < 0)    error("invalid 'target' argument");
  supp   = get_dbl(psupp, supp);
  zmin   = get_int(pzmin,  zmin);
  if (zmin   < 0)    error("invalid 'zmin' argument (must be >= 0)");
  zmax   = get_int(pzmax,  zmax);
  if (zmax   < 0)    zmax = ITEM_MAX;  /* check the size range */
  if (zmax   < zmin) error("invalid 'zmax' argument (must be >= zmin)");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  report = get_str(preport, report);
  equiv  = get_int(pequiv,  equiv);
  if (equiv <= 0)    equiv = 1; /* check the number of data sets */
  alpha  = get_dbl(palpha,  alpha);
  if (alpha <= 0)    error("invalid 'alpha' argument");
  smpls  = get_int(psmpls,  smpls);
  if (smpls <= 0)    error("invalid 'smpls' argument");
  seed   = get_lng(pseed,   seed);

  /* --- estimate pattern spectrum --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromRObj(ptracts, pwgts, R_NilValue);
  if (!tabag) ERR_MEM();        /* create a transaction bag and */
  wgt  = tbg_wgt(tabag);        /* get the total transaction weight */
  supp = (supp >= 0) ? supp/100.0 *(double)wgt *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support value */
  if (tbg_recode(tabag, smin, -1, -1, -2) < 0) {
    tbg_delete(tabag, 1); ERR_MEM(); }
  tbg_filter(tabag, (ITEM)zmin, NULL, 0);
  psp = fpg_estpsp(tabag, target, (SUPP)smin, (ITEM)zmin, (ITEM)zmax,
                   (size_t)equiv, alpha, (size_t)smpls, seed);
  if (psp) { rpsp = psp_toRObj(psp, 1.0/(double)equiv, report[0]);
             psp_delete(psp); } /* generate a pattern spectrum */
  tbg_delete(tabag, 1);         /* delete the transaction bag */
  if (!rpsp) ERR_MEM();         /* check for an error, otherwise */
  else       UNPROTECT(1);      /* unprotect the created object */
  if (sig_aborted()) { sig_abort(0); ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  return rpsp;                  /* return the created R object */
}  /* f4r_estpsp() */

/*--------------------------------------------------------------------*/
/* patred (pats, method="S", border=NULL, addis=TRUE)                 */
/*--------------------------------------------------------------------*/

SEXP f4r_patred (SEXP ppats, SEXP pmethod, SEXP pborder, SEXP paddis)
{                               /* --- pattern set reduction */
  int    method = PSR_COVER1;   /* pattern set reduction method */
  int    addis  = -1;           /* whether to add pairwise isects. */
  IDMAP  *map;                  /* item to identifier map */
  PATSET *patset;               /* internal pattern set */
  size_t i, n;                  /* loop variable, number of patterns */
  size_t k, z;                  /* (maximum) size of a pattern */
  size_t x;                     /* pattern extent (item instances) */
  SEXP   rset;                  /* to traverse the R item sets */
  SEXP   items;                 /* to traverse the item sets */
  SEXP   supp;                  /* (vector of) support value(s) */
  int    ti, ts;                /* type of items and support */
  RSUPP  s;                     /* support of item set or in border */
  CCHAR  *p;                    /* to traverse string items */
  int    *b;                    /* to access an integer border */
  double *d;                    /* to access a real-valued border */

  /* --- evaluate function arguments --- */
  if (ppats == R_NilValue) return ppats;
  if (TYPEOF(ppats) != VECSXP)  /* pattern set must be a vector/array */
    error("invalid 'pats' argument (must be a list)");
  n = (size_t)length(ppats);    /* get the number of patterns */
  if (n <= 0) return ppats;     /* and return if there are none */
  method = get_red(pmethod);    /* get and check the reduction method */
  if (method < 0) error("invalid 'method' argument");
  addis = get_lgl(paddis,addis);/* get flag for adding intersections */
  rset = VECTOR_ELT(ppats, 0);  /* get and check first pattern */
  if ((TYPEOF(rset) != VECSXP) || (length(rset) < 2))
    error("invalid 'pats' argument (wrong pattern type)");
  items = VECTOR_ELT(rset, 0);  /* get items of first pattern */
  ti    = TYPEOF(items);        /* get and check the item type */
  if ((ti != INTSXP) && (ti != STRSXP))
    error("invalid 'pats' argument (wrong item type)");
  supp  = VECTOR_ELT(rset, 1);  /* get support of first pattern */
  ts    = TYPEOF(supp);         /* get and check the support type */
  if (((ts != INTSXP) && (ts != REALSXP)) || (length(supp) <= 0))
    error("invalid 'pats' argument (wrong support type)");
  k = 1; x = 0;                 /* init. pattern size and extent */
  for (i = 0; i < n; i++) {     /* traverse the patterns */
    rset = VECTOR_ELT(ppats, (R_xlen_t)i);
    if ((TYPEOF(rset) != VECSXP) || (length(rset) < 2))
      error("invalid 'pats' argument (wrong pattern type)");
    items = VECTOR_ELT(rset,0); /* get items of next pattern */
    if (TYPEOF(items) != ti)    /* and check the item type */
      error("invalid 'pats' argument (wrong item type)");
    supp  = VECTOR_ELT(rset,1); /* get support of next pattern */
    if ((TYPEOF(supp) != ts) || (length(supp) <= 0))
      error("invalid 'pats' argument (wrong support type)");
    x += z = (size_t)length(items);
    if (z > k) k = z;           /* get the pattern size and */
  }                             /* update extent and maximal size */
  /* n: number of patterns (item set/support pairs) */
  /* k: maximal size of a pattern (number of items) */
  /* x: total number of item istances (extent) */

  /* --- create pattern set --- */
  z   = (n > 255) ? n : 255;    /* compute identifier map size */
  map = (ti == INTSXP)          /* according to the item type */
      ? idm_create(z, 0, ST_INTFN, (OBJFN*)0)
      : idm_create(z, 0, ST_STRFN, (OBJFN*)0);
  if (!map) ERR_MEM();          /* create item map and pattern set */
  patset = psr_create(n, (ITEM)k, x, map);
  if (!patset) { idm_delete(map); ERR_MEM(); }
  if ((pborder != R_NilValue) && (length(pborder) > 0)) {
    z = (size_t)length(pborder);/* get the border length */
    if (z > k+1) z = k+1;       /* limit border to pattern size */
    if      (TYPEOF(pborder) == INTSXP) {
      b = INTEGER(pborder);     /* get the support values */
      for (i = 2; i < z; i++) { /* traverse the pattern sizes */
        s = (b[i] < 0)    ? RSUPP_MAX : (RSUPP)b[i];
        psr_setbdr(patset, (ITEM)i, s);
      } }                       /* store support values in border */
    else if (TYPEOF(pborder) == REALSXP) {
      d = REAL(pborder);        /* get the support values */
      for (i = 2; i < z; i++) { /* traverse the pattern sizes */
        s = (IS_NA(d[i])) ? RSUPP_MAX : (RSUPP)d[i];
        psr_setbdr(patset, (ITEM)i, s);
      } }                       /* store support values in border */
    else error("invalid 'border' argument (must be numeric or NULL)");
  }                             /* check the border type */
  if (sig_aborted()) { sig_abort(0); psr_delete(patset,1); ERR_ABORT();}

  /* --- collect patterns --- */
  for (i = 0; i < n; i++) {     /* collect the item sets */
    rset  = VECTOR_ELT(ppats, (R_xlen_t)i);
    psr_addorig(patset, rset);  /* note the input pattern */
    items = VECTOR_ELT(rset,0); /* get the items in the set */
    supp  = VECTOR_ELT(rset,1); /* and the support value */
    z = (size_t)length(items);  /* get the number of items */
    if (ti == INTSXP) {         /* if items are integers */
      for (k = 0; k < z; k++) { /* traverse the items */
        if (psr_additem(patset, INTEGER(items) +k) != 0) {
          psr_delete(patset, 1); ERR_MEM(); }
      } }                       /* add items to pattern set */
    else {                      /* if items are strings */
      for (k = 0; k < z; k++) { /* traverse the items */
        p = CHAR(STRING_ELT(items, (R_xlen_t)k));
        if (psr_additem(patset, p) != 0) {
          psr_delete(patset, 1); ERR_MEM(); }
      }                         /* add items to pattern set */
    }                           /* (convert to integers internally) */
    s = (ts == INTSXP) ? (RSUPP)INTEGER(supp)[0] : (RSUPP)REAL(supp)[0];
    psr_addsupp(patset, s);     /* get the pattern support and */
  }                             /* store pattern in pattern set */
  if (sig_aborted()) { sig_abort(0); psr_delete(patset,1); ERR_ABORT();}

  /* --- pattern set reduction --- */
  k     = psr_reduce(patset, method, addis);
  ppats = PROTECT(allocVector(VECSXP, (R_xlen_t)k));
  for (i = k = 0; i < n; i++) { /* traverse the reduced patterns */
    rset = (SEXP)psr_getorig(patset, i);
    if (rset) { SET_VECTOR_ELT(ppats, (R_xlen_t)k, rset); k++; }
  }                             /* collect the kept patterns */
  psr_delete(patset, 0);        /* clean up allocated memory */
  UNPROTECT(1);                 /* release the new pattern list */
  return ppats;                 /* return modified pattern set */
}  /* f4r_patred() */
