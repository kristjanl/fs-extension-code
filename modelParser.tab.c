/* A Bison parser, made by GNU Bison 3.0.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "modelParser.y" /* yacc.c:339  */

	/*---
	Flow*: A Verification Tool for Cyber-Physical Systems.
	Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
	Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

	The code is released as is under the GNU General Public License (GPL).
	---*/


	#include "modelParser.h"
	#include "MyLogger.h"

	extern int yyerror(const char *);
	extern int yyerror(string);
	extern int yylex();
	extern int yyparse();
	bool err;

#line 86 "modelParser.tab.c" /* yacc.c:339  */

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "modelParser.tab.h".  */
#ifndef YY_YY_MODELPARSER_TAB_H_INCLUDED
# define YY_YY_MODELPARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NUM = 258,
    IDENT = 259,
    STATEVAR = 260,
    TMVAR = 261,
    TM = 262,
    EQ = 263,
    GEQ = 264,
    LEQ = 265,
    ASSIGN = 266,
    END = 267,
    MODE = 268,
    INIT = 269,
    BELONGSTO = 270,
    POLYODE1 = 271,
    POLYODE2 = 272,
    POLYODE3 = 273,
    VISUALIZE = 274,
    PARAAGGREG = 275,
    INTAGGREG = 276,
    TMAGGREG = 277,
    OUTPUT = 278,
    CONTINUOUS = 279,
    HYBRID = 280,
    SETTING = 281,
    FIXEDST = 282,
    FIXEDORD = 283,
    ADAPTIVEST = 284,
    ADAPTIVEORD = 285,
    MIN = 286,
    MAX = 287,
    REMEST = 288,
    INTERVAL = 289,
    OCTAGON = 290,
    GRID = 291,
    QRPRECOND = 292,
    IDPRECOND = 293,
    SHRINRWRAPPING = 294,
    TIME = 295,
    MODES = 296,
    JUMPS = 297,
    INV = 298,
    GUARD = 299,
    RESET = 300,
    START = 301,
    MAXJMPS = 302,
    PRINTON = 303,
    PRINTOFF = 304,
    UNSAFESET = 305,
    CONTINUOUSFLOW = 306,
    HYBRIDFLOW = 307,
    TAYLOR_PICARD = 308,
    TAYLOR_REMAINDER = 309,
    TAYLOR_POLYNOMIAL = 310,
    NONPOLY_CENTER = 311,
    EXP = 312,
    SIN = 313,
    COS = 314,
    LOG = 315,
    SQRT = 316,
    NPODE_TAYLOR = 317,
    CUTOFF = 318,
    PRECISION = 319,
    GNUPLOT = 320,
    MATLAB = 321,
    COMPUTATIONPATHS = 322,
    LINEARODE = 323,
    PAR = 324,
    METHOD = 325,
    ALGORITHM = 326,
    ALG_FLOW = 327,
    ALG_SIMPLE_IMPL = 328,
    ALG_SIMPLE_COMP = 329,
    ALG_SMALL_COMP = 330,
    DECOMPOSITION = 331,
    NODECOMPOSITION = 332,
    uminus = 333
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 22 "modelParser.y" /* yacc.c:355  */

	double dblVal;
	string *identifier;
	vector<Interval> *intVec;
	vector<int> *iVec;
	vector<double> *dVec;
	vector<Monomial> *monoVec;
	vector<Polynomial> *polyVec;
	Monomial *mono;
	Polynomial *poly;
	TaylorModelVec *tmVec;
	Matrix *mat;
	vector<vector<double> > *dVecVec;
	vector<PolynomialConstraint> *vecConstraints;
	ResetMap *resetMap;
	Flowpipe *pFlowpipe;
	TaylorModel *ptm;
	Interval *pint;
	vector<string> *strVec;
	TreeNode *pNode;

#line 227 "modelParser.tab.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_MODELPARSER_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 242 "modelParser.tab.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef __attribute__
/* This feature is available in gcc versions 2.5 and later.  */
# if (! defined __GNUC__ || __GNUC__ < 2 \
      || (__GNUC__ == 2 && __GNUC_MINOR__ < 5))
#  define __attribute__(Spec) /* empty */
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  18
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1511

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  95
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  60
/* YYNRULES -- Number of rules.  */
#define YYNRULES  277
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1046

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   333

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    93,
      87,    91,    80,    78,    88,    79,     2,    81,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    94,    86,
       2,     2,    92,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    89,     2,    90,    83,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    84,     2,    85,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    82
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   117,   117,   166,   236,   283,   352,   392,   408,   441,
     450,   457,   463,   469,   476,   485,   495,   514,   534,   555,
     577,   581,   586,   621,   655,   661,   668,   675,   680,   696,
     706,   722,   741,   763,   782,   813,   847,   884,   935,   940,
     956,   975,   997,  1016,  1047,  1081,  1118,  1169,  1174,  1190,
    1200,  1210,  1220,  1230,  1240,  1250,  1260,  1270,  1281,  1292,
    1304,  1339,  1358,  1368,  1378,  1388,  1398,  1408,  1418,  1428,
    1438,  1448,  1461,  1465,  1469,  1481,  1497,  1502,  1508,  1513,
    1529,  1541,  1561,  1586,  1613,  1654,  1658,  1689,  1713,  1718,
    1745,  1766,  1782,  1788,  1793,  1834,  1851,  1871,  1876,  1895,
    1919,  1924,  1949,  1977,  2031,  2091,  2154,  2232,  2291,  2363,
    2382,  2400,  2424,  2451,  2468,  2492,  2500,  2507,  2517,  2528,
    2539,  2545,  2551,  2558,  2562,  2572,  2580,  2582,  2584,  2610,
    2635,  2655,  2691,  2727,  2765,  2801,  2837,  2876,  2885,  2896,
    2901,  2935,  2974,  3002,  3013,  3038,  3068,  3096,  3119,  3141,
    3156,  3176,  3193,  3215,  3235,  3243,  3251,  3256,  3264,  3282,
    3289,  3317,  3335,  3344,  3352,  3357,  3365,  3383,  3390,  3434,
    3442,  3456,  3464,  3472,  3477,  3485,  3503,  3510,  3552,  3563,
    3571,  3579,  3584,  3592,  3610,  3617,  3659,  3671,  3698,  3743,
    3751,  3759,  3764,  3772,  3790,  3797,  3825,  3832,  3859,  3866,
    3873,  3889,  3894,  3913,  3922,  3931,  3940,  3949,  3958,  4005,
    4012,  4035,  4058,  4065,  4072,  4088,  4093,  4112,  4121,  4130,
    4139,  4148,  4157,  4213,  4219,  4240,  4262,  4270,  4278,  4288,
    4293,  4306,  4315,  4324,  4333,  4342,  4351,  4361,  4367,  4385,
    4405,  4416,  4427,  4438,  4449,  4460,  4471,  4482,  4493,  4504,
    4515,  4527,  4537,  4565,  4584,  4618,  4629,  4640,  4651,  4662,
    4673,  4684,  4695,  4706,  4717,  4728,  4740,  4750,  4765,  4782,
    4813,  4821,  4829,  4846,  4851,  4858,  4900,  4907
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "IDENT", "STATEVAR", "TMVAR",
  "TM", "EQ", "GEQ", "LEQ", "ASSIGN", "END", "MODE", "INIT", "BELONGSTO",
  "POLYODE1", "POLYODE2", "POLYODE3", "VISUALIZE", "PARAAGGREG",
  "INTAGGREG", "TMAGGREG", "OUTPUT", "CONTINUOUS", "HYBRID", "SETTING",
  "FIXEDST", "FIXEDORD", "ADAPTIVEST", "ADAPTIVEORD", "MIN", "MAX",
  "REMEST", "INTERVAL", "OCTAGON", "GRID", "QRPRECOND", "IDPRECOND",
  "SHRINRWRAPPING", "TIME", "MODES", "JUMPS", "INV", "GUARD", "RESET",
  "START", "MAXJMPS", "PRINTON", "PRINTOFF", "UNSAFESET", "CONTINUOUSFLOW",
  "HYBRIDFLOW", "TAYLOR_PICARD", "TAYLOR_REMAINDER", "TAYLOR_POLYNOMIAL",
  "NONPOLY_CENTER", "EXP", "SIN", "COS", "LOG", "SQRT", "NPODE_TAYLOR",
  "CUTOFF", "PRECISION", "GNUPLOT", "MATLAB", "COMPUTATIONPATHS",
  "LINEARODE", "PAR", "METHOD", "ALGORITHM", "ALG_FLOW", "ALG_SIMPLE_IMPL",
  "ALG_SIMPLE_COMP", "ALG_SMALL_COMP", "DECOMPOSITION", "NODECOMPOSITION",
  "'+'", "'-'", "'*'", "'/'", "uminus", "'^'", "'{'", "'}'", "';'", "'('",
  "','", "'['", "']'", "')'", "'>'", "'\\''", "':'", "$accept", "model",
  "continuous_flowpipes", "modeDecls", "hybrid_flowpipes",
  "computation_paths", "computation_path", "print", "unsafe_continuous",
  "unsafe_hybrid", "hybrid_constraints", "polynomial_constraints",
  "linear_constraints", "continuous", "hybrid", "hybrid_init", "modes",
  "local_setting", "parameters", "jumps", "reset", "real_valued_vectors",
  "real_valued_vector", "vector_components", "stateVarDecls",
  "stateIdDeclList", "parDecls", "parDeclList", "settings",
  "remainder_estimation", "remainders", "orders", "precondition",
  "algorithm", "decomposition", "components", "component", "compVarIds",
  "plotting", "init", "tmVarDecls", "tmIdDeclList", "taylor_model",
  "taylor_model_domain", "intervals", "ode", "npode", "linear_ode",
  "polynomial", "ODEpolynomial", "constraint_polynomial",
  "reset_polynomial", "interval_taylor_model", "interval_polynomial",
  "non_polynomial_rhs_picard", "non_polynomial_rhs_remainder",
  "non_polynomial_rhs_no_remainder", "non_polynomial_rhs_string",
  "non_polynomial_rhs_center", "linear_polynomial", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,    43,    45,
      42,    47,   333,    94,   123,   125,    59,    40,    44,    91,
      93,    41,    62,    39,    58
};
# endif

#define YYPACT_NINF -595

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-595)))

#define YYTABLE_NINF -189

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     475,    10,   -61,   -33,   -23,   -19,     0,     4,   118,    24,
    -595,    71,   117,   117,   112,   121,   128,    12,  -595,    80,
     798,   871,    14,   134,   210,   176,    38,   223,   126,  -595,
     225,   240,   255,   270,   293,   112,   112,   328,   456,  -595,
     313,   319,   339,   356,   379,   121,   121,   333,   760,  -595,
     431,   440,   446,   467,   473,   128,   128,   490,   768,  -595,
    -595,   484,   486,   506,   508,   510,    12,    12,   559,   776,
     262,   573,   604,   607,   612,   618,   621,   542,   553,   636,
    -595,   591,   568,   572,   628,   617,   587,   643,   112,   112,
     112,   112,   112,   599,   472,   596,   112,   112,   112,   112,
     683,  -595,   121,   121,   121,   121,   121,   614,   505,   611,
     121,   121,   121,   121,   698,  -595,   128,   128,   128,   128,
     128,   629,   511,   626,   128,   128,   128,   128,   713,  -595,
      12,    12,    12,    12,    12,   644,   520,   641,    12,    12,
      12,    12,   728,  -595,   679,   656,   658,   753,   671,   673,
     770,   688,   772,   764,   707,  -595,   256,   785,   718,   722,
    -595,   256,   726,   534,   540,   549,   555,   564,  -595,   801,
     321,   321,   599,   599,  -595,   570,   579,   585,   594,   600,
    -595,   818,   384,   384,   614,   614,  -595,   609,   615,   624,
     630,   639,  -595,   839,   443,   443,   629,   629,  -595,   645,
     654,   660,   669,   675,  -595,   841,   487,   487,   644,   644,
    -595,   847,   848,   854,   778,   856,   898,   820,   870,  -595,
       9,   180,   930,  -595,   932,   852,     7,   929,    13,   256,
    -595,   891,   256,  -595,  -595,  -595,  -595,  -595,   849,  -595,
    -595,  -595,  -595,  -595,   850,  -595,  -595,  -595,  -595,  -595,
     851,  -595,  -595,  -595,  -595,  -595,   853,   857,  -595,  -595,
     940,  -595,  -595,   941,   943,   323,   348,  -595,   944,   398,
      41,   908,   921,  -595,  -595,   864,   947,   945,  -595,     7,
      15,   948,   907,  -595,  -595,  -595,  -595,   873,  -595,  -595,
     874,   933,  -595,   867,   876,   910,  -595,  -595,   220,  -595,
     220,   375,   955,   959,   457,  -595,   960,   879,   881,  -595,
       7,   963,  -595,   883,   964,   880,   965,   886,   889,   -31,
     453,   507,   541,   884,   220,   220,   220,   971,   942,   888,
     893,   894,   895,   896,   897,  -595,   470,  -595,   899,     7,
     111,  -595,   254,   979,   981,   902,   965,  -595,  -595,  -595,
    -595,  -595,  -595,  -595,   868,   -44,   -44,   889,  -595,    19,
     872,   951,  -595,  -595,  -595,  -595,  -595,   903,   904,   905,
     906,   911,   141,   950,   912,   913,   146,   915,   949,   914,
    -595,   916,   988,   793,   902,   917,   923,  -595,   989,  -595,
    -595,   991,   323,   993,    20,    26,    39,    42,    46,  -595,
    -595,  -595,  -595,  -595,  -595,   925,   962,  -595,   927,   996,
     926,  1003,  1004,   999,  1009,   988,  -595,   799,  1011,  1012,
     928,   151,  -595,   291,   934,   931,  1001,  1006,  1007,   935,
    1013,   936,  1016,    55,    62,    82,   156,   170,  1014,   939,
    -595,   952,   173,   996,   953,  -595,    30,    77,   179,  1009,
    -595,   954,   956,   937,  -595,  1021,   181,   957,   977,  1023,
     958,   961,   966,  1024,   967,  1025,   968,  1020,  1026,  1033,
    1034,  1035,  1027,   190,  1014,   965,   970,  -595,   226,   946,
    -595,  -595,    30,    30,  1032,   526,    30,   969,  1040,  -595,
     230,  -595,  -595,  1036,   972,   975,  1052,  1028,  1054,    34,
     965,   965,   965,    17,   965,    50,   965,   976,   978,   980,
     983,   984,   985,   986,  1027,  1019,   244,   902,   965,  -595,
     992,   982,   478,   987,    54,    30,    30,  1060,   790,  1069,
     990,  -595,   994,   995,  1070,   997,   228,   249,   942,  -595,
    -595,    34,    34,  1071,   796,   998,  -595,  1072,  1000,  1002,
    -595,  -595,  1005,  1008,  1010,  1015,  1017,    17,    17,  1074,
     784,  1018,  -595,  -595,    50,    50,  1075,   834,  1022,   965,
     965,   965,   965,   965,  -595,   497,  1029,  1030,  1038,   831,
     902,  1031,    30,  -595,  1078,  1083,   252,   252,   982,  -595,
      59,  1037,  1085,  1086,  1087,  1041,  1090,  1043,  1092,  1039,
    1052,   872,  1045,   729,  1042,    34,    34,    34,  1091,  -595,
    1094,  1084,  -595,  -595,    17,    17,    17,    17,    17,  1046,
     684,  1044,    17,    17,    17,    17,  1097,  -595,  -595,   496,
    1047,    50,    50,    50,  -595,  1048,  1049,  1051,  1053,  1055,
     479,  1057,  1058,  1059,  1061,  1062,   504,  -595,  1063,  -595,
     833,  1104,  1064,  1056,  1098,  1106,  1065,  1066,  1067,  1107,
    -595,  1108,  1068,  1080,   280,   323,  -595,  1112,   460,   460,
    1045,  -595,   198,  1113,  1050,   690,   699,   705,   714,   720,
    -595,  1115,   736,   736,  1046,  1046,  -595,  -595,  1116,  1077,
    1077,  -595,  -595,  -595,  -595,  -595,  -595,  1117,   253,  1076,
    1079,   257,  -595,  -595,  1118,  -595,  -595,  -595,  -595,  -595,
    -595,  1081,  1082,  1088,  1089,  1093,   258,  -595,  -595,  -595,
    -595,  1119,  1073,  1095,  1121,  -595,  1123,  1114,  1100,  1128,
    1145,  1096,  1122,  1099,   241,  1146,  -595,  -595,  -595,  -595,
    -595,  1101,  1102,  -595,  -595,  1052,  1120,  1127,  -595,   989,
    -595,   259,   265,   266,   267,   268,  -595,  -595,  -595,  -595,
    -595,  1103,  1138,   269,  1105,  1156,  -595,  1109,  1110,  1163,
    1165,  -595,  1111,  1137,   274,  -595,  -595,  -595,   241,   241,
     802,  1124,  -595,  -595,   743,  1167,   276,   804,  1131,  1132,
    1133,  1135,  1136,   277,   278,   284,   287,   288,  1125,  1126,
    1157,  1176,  1129,  -595,  -595,   343,  1158,  1130,  1134,  1139,
    1052,  1140,   735,    76,   241,   241,  1180,  1183,  -595,  1141,
    1142,  1052,  -595,  1143,  1144,  1147,  1148,  1149,  1151,  1154,
    1155,  1160,  1161,  1184,  1186,  1150,  1197,  -595,  -595,  -595,
    -595,   843,  1201,  1203,  1052,  1204,   805,   241,  -595,  1205,
     814,   814,  1140,  -595,  1152,  1177,  1179,   810,  -595,  -595,
    -595,  -595,  -595,  1153,  1159,  1162,  1164,  1166,  1169,  1168,
    1170,  1186,  1171,  -595,  -595,  -595,  1172,   811,  1174,  1178,
    1173,  -595,  1211,  1212,  1175,   150,   153,   172,   189,     8,
    -595,  -595,  -595,  -595,  -595,  1181,   965,  -595,  1182,  1185,
    1213,  1187,  1217,  1218,  1219,  1188,  1190,  1192,  1191,  1193,
    1194,  1195,  1196,   330,   199,   208,   215,   218,    23,  -595,
    -595,   292,  -595,  1221,   235,  -595,  1216,  1199,  1222,  1189,
    1198,  -595,  -595,  1200,  -595,  -595,  -595,  -595,  -595,   918,
     920,   922,  1202,  1207,  1208,  1209,  1210,  1214,   239,  1094,
    -595,  -595,   271,  1185,  -595,  1231,  1223,  1236,  1241,  -595,
    1052,  -595,  -595,  -595,  -595,  -595,  -595,   924,  -595,  -595,
    -595,  -595,  -595,  1206,   294,  1243,  -595,  -595,  -595,  1225,
    -595,  1226,   815,  1215,  1220,  1227,  -595,  -595,  1251,  1252,
    1224,  1254,  1255,  -595,  1235,  -595,  -595,  1228,  1229,   297,
    1260,  -595,  -595,  1230,   909,  -595,  1257,  1232,  -595,   243,
    -595,  -595,  -595,   243,   243,   808,   -54,  1234,   744,    79,
     243,   243,  1263,  -595,  1265,  -595,   243,  -595,  1267,   821,
     821,  1234,  -595,  1233,   305,  1237,  1268,  1273,  -595,  1279,
    -595,  1238,  1239,  1280,  -595,  -595
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      99,    97,     0,     0,     0,     0,     0,     0,     1,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   210,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   224,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   238,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   269,
     267,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      98,     2,     0,     0,     0,     4,     0,     0,     0,     0,
       0,     0,     0,   209,     0,     0,     0,     0,     0,     0,
       0,    10,     0,     0,     0,     0,     0,   223,     0,     0,
       0,     0,     0,     0,     0,    11,     0,     0,     0,     0,
       0,   237,     0,     0,     0,     0,     0,     0,     0,    12,
       0,     0,     0,     0,     0,   266,     0,     0,     0,     0,
       0,     0,     0,    13,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     3,     0,     0,     0,     0,
       5,     0,     0,     0,     0,     0,     0,     0,   201,     0,
     198,   199,   200,   202,   208,     0,     0,     0,     0,     0,
     215,     0,   212,   213,   214,   216,   222,     0,     0,     0,
       0,     0,   229,     0,   226,   227,   228,   230,   236,     0,
       0,     0,     0,     0,   258,     0,   255,   256,   257,   259,
     265,     0,     0,     0,     0,     0,     0,     0,     0,    23,
       0,     0,     0,    38,     0,     0,     0,     0,     0,     0,
      29,     0,     0,   203,   204,   205,   206,   207,     0,   217,
     218,   219,   220,   221,     0,   231,   232,   233,   234,   235,
       0,   260,   261,   262,   263,   264,     0,     0,   131,   132,
       0,   134,   135,     0,     0,     0,     0,    21,     0,     0,
       0,     0,     0,    24,    25,     0,     0,     0,   100,     0,
       0,     0,     0,   211,   225,   239,   268,     0,   133,   136,
       0,     0,    20,     0,     0,     0,   178,   177,     0,    26,
       0,     0,     0,     0,     0,   102,     0,     0,     0,    27,
       0,     0,    38,     0,     0,     0,     0,     0,   176,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   101,     0,    38,     0,     0,
       0,    38,     0,     0,     0,     0,     0,   173,    32,    36,
      31,    35,    30,    34,     0,   171,   172,   174,   175,     0,
       0,     0,   149,   149,   149,   151,   153,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     141,   139,     0,     0,     0,     0,     0,   109,     0,   115,
     116,     0,     0,     0,     0,     0,     0,     0,     0,   149,
     149,   149,   151,   153,    28,     0,     0,    17,     0,     0,
       0,     0,     0,     0,     0,     0,     7,     0,     0,     0,
       0,     0,   117,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      16,     0,     0,     0,     0,   140,     0,     0,     0,     0,
       6,     0,     0,     0,   110,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    73,     0,     0,     0,     0,     9,     0,     0,
     197,   195,     0,     0,     0,     0,     0,     0,     0,    15,
       0,    33,    37,     0,     0,     0,     0,     0,     0,     0,
     147,   147,   147,     0,   147,     0,   147,     0,     0,     0,
       0,     0,     0,     0,    73,     0,     0,     0,     0,     8,
       0,   194,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    14,     0,     0,     0,     0,     0,     0,     0,   170,
     168,     0,     0,     0,   148,     0,   143,   138,     0,     0,
     254,   252,     0,     0,     0,     0,     0,     0,     0,     0,
     150,     0,   276,   275,     0,     0,     0,   152,     0,   147,
     147,   147,   147,   147,    85,     0,     0,     0,     0,     0,
       0,     0,     0,   191,     0,     0,   189,   190,   192,   193,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   167,     0,     0,     0,     0,     0,     0,    48,
       0,     0,    49,    50,     0,     0,     0,     0,     0,   251,
       0,     0,     0,     0,     0,     0,     0,    51,   274,     0,
       0,     0,     0,     0,    52,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    88,     0,    19,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     114,     0,     0,     0,     0,     0,   164,     0,   162,   163,
     165,   166,     0,   137,     0,     0,     0,     0,     0,     0,
     243,     0,   240,   241,   242,   244,   250,   273,     0,   270,
     271,   272,    53,    54,    55,    56,    57,     0,     0,     0,
       0,     0,    77,    78,     0,    72,   149,   149,   149,   151,
     153,     0,     0,     0,     0,     0,     0,    88,    18,    22,
     196,     0,     0,     0,     0,   112,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   245,   246,   247,   248,
     249,     0,     0,    74,    79,     0,     0,     0,    75,     0,
      80,     0,     0,     0,     0,     0,   149,   149,   149,   151,
     153,     0,     0,     0,     0,     0,   145,     0,     0,     0,
       0,   113,     0,     0,     0,   169,   161,   160,     0,     0,
       0,     0,   253,   277,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   196,     0,   144,   111,   122,     0,     0,     0,     0,
       0,   159,     0,     0,     0,     0,     0,     0,    83,     0,
       0,     0,    76,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   196,   121,   118,   119,
     120,   125,     0,     0,     0,     0,     0,     0,   156,     0,
     154,   155,   157,   158,     0,     0,     0,     0,    38,    38,
      38,    38,    47,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   124,   103,   105,     0,     0,     0,     0,
       0,   146,     0,     0,     0,     0,     0,     0,     0,     0,
      38,    38,    38,    38,    47,     0,   147,    58,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    38,
     143,     0,    59,     0,     0,   127,     0,     0,     0,     0,
       0,    81,    82,     0,    63,    65,    67,    69,    71,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      60,   130,     0,     0,   123,     0,     0,     0,     0,   142,
       0,    41,    45,    40,    44,    39,    43,     0,    62,    64,
      66,    68,    70,     0,     0,     0,   128,   126,   104,     0,
     107,     0,     0,     0,     0,     0,    61,   129,     0,     0,
       0,     0,     0,    91,     0,   108,    84,     0,     0,     0,
       0,    42,    46,     0,     0,   106,     0,     0,    87,     0,
      93,   186,   185,     0,     0,    90,     0,   184,     0,     0,
       0,     0,     0,    86,     0,    92,     0,   181,     0,   179,
     180,   182,   183,     0,     0,     0,     0,     0,    94,     0,
      96,     0,     0,     0,    89,    95
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -595,  -595,  -375,  -595,   842,  -595,  1240,   -30,  1242,   774,
    -595,  -312,   392,  -595,  -595,   316,   813,   782,  -595,   580,
    -595,  -595,  -595,  -595,   919,  -595,  1262,  -595,   149,   762,
     552,  -594,   697,  -595,  -595,  -595,   349,  -595,  -262,  -110,
    -314,  -595,   385,  -448,   408,  -359,  -395,  -393,  -525,  -155,
     107,  -184,   892,  -112,   406,   366,   320,  -133,   289,  -497
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     8,   383,    22,   442,   220,   221,   275,   155,   160,
     280,   270,   889,    25,    27,   870,   473,   513,   640,   716,
     999,  1016,  1025,  1034,     9,    11,    84,   228,   226,   360,
     421,   536,   392,   841,   874,   924,   925,   952,    23,   545,
     546,   381,   610,   448,   547,   394,   397,   398,   850,   544,
     301,  1029,   414,   586,    38,    48,    58,   560,    69,   913
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     340,   490,   345,   291,   395,   396,   664,   436,   567,   417,
     437,   562,   563,   219,    10,    59,    60,   277,    77,   308,
     550,   551,   387,    12,   425,   372,   562,   563,    19,   376,
     425,  1023,   384,   480,   481,  1024,   326,   539,   540,   327,
     433,   434,   435,   425,   296,   297,   429,   324,   325,   326,
     431,    13,   327,   562,   563,   273,   274,   480,   481,   425,
     347,    14,   480,   481,    82,    15,   425,   628,   629,    61,
      62,    63,    64,    65,   552,   553,   554,   555,   556,   776,
     777,    78,  1011,  1012,    16,   486,   425,   564,    17,    20,
      21,    66,   487,   912,   265,   565,   557,   566,   278,    67,
     309,    68,   564,   388,   558,   426,   559,    83,   947,   482,
     565,   427,   566,   541,   296,   297,    29,   483,    18,   484,
     298,   542,     1,   543,   428,    39,   299,   430,   300,   564,
     423,   432,    49,   482,   689,   690,   691,   565,   482,   566,
     467,   483,   579,   585,   296,   297,   483,   468,   654,   296,
     297,   784,    86,   296,   297,   778,   296,   297,  1013,    24,
     429,   517,   673,   779,    70,   849,  1014,   469,  1028,    30,
      31,    32,    33,    34,   431,   296,   297,   476,    40,    41,
      42,    43,    44,   488,   495,    50,    51,    52,    53,    54,
     298,    35,   296,   297,   514,    83,   375,    79,   300,    36,
      45,    37,   296,   297,   580,   650,   734,    55,    46,   780,
      47,   296,   297,   487,    80,    56,   846,    57,   296,   297,
     298,   296,   297,   296,   297,   298,   404,   857,   300,   298,
     476,   408,   298,   300,   488,   908,   454,   300,   909,   455,
     300,   470,   296,   297,   776,   777,  1011,  1012,   514,   307,
     877,   298,   599,   811,   812,   471,   744,   910,   477,   300,
     748,    81,   761,   425,   489,   496,   267,   268,   298,   425,
     425,   429,   431,   761,   911,   515,   300,   809,   298,   820,
     338,   425,   425,   224,   943,   225,   300,   298,   425,   851,
     852,   429,   431,   944,   298,   300,   611,   298,   488,   298,
     945,  1003,   300,   946,   159,   300,   377,   300,    85,   374,
     231,   519,    88,   597,   754,   531,   598,   755,   298,   456,
     778,   457,  1013,   953,   973,   954,   300,    89,   779,   578,
    1014,    95,   526,   600,   485,   527,   109,   745,   939,   940,
     941,   749,    90,   762,   788,   942,   144,   751,   752,   753,
     789,   790,   791,   792,   800,   135,   136,    91,   810,   975,
     821,   976,   828,   829,   796,   731,   982,   797,   598,   830,
     521,   522,   831,   832,   528,   121,   122,   950,   279,   986,
      92,   282,  1004,   320,   321,   322,   602,   603,    20,    21,
     323,   548,   549,  1037,   561,  1038,   568,   793,   794,   795,
     102,    98,    99,   732,   100,   318,   103,   319,   631,   632,
     633,   107,   108,   587,   588,   837,   838,   839,   840,   199,
     200,   201,   202,   203,   619,   620,   104,   206,   207,   208,
     209,   355,   356,   357,   292,   268,   187,   188,   189,   190,
     191,    93,    94,   105,   194,   195,   196,   197,   154,   294,
     668,   669,   670,   324,   325,   326,   348,   349,   327,   635,
     636,   637,   638,   639,   112,   113,   106,   114,   175,   176,
     177,   178,   179,   330,   331,   332,   182,   183,   184,   185,
       1,   675,   676,   677,   678,   679,   367,   368,   369,   682,
     683,   684,   685,   123,   163,   164,   165,   166,   167,     2,
       3,   974,   170,   171,   172,   173,   697,   698,   699,   700,
     350,   351,   701,   641,   642,   643,   702,   703,   116,   333,
     711,   712,   713,   126,   127,   334,   128,   117,     4,     5,
       6,     7,   370,   118,    96,    97,    98,    99,   371,   100,
     607,   101,   704,   608,   352,   353,   885,   886,   887,   888,
      96,    97,    98,    99,   119,   100,   582,   525,   526,   644,
     120,   527,   137,   168,   705,   645,   714,   140,   141,   583,
     142,   130,   715,   131,   631,   632,   633,   145,   914,   915,
     916,   917,   920,   110,   111,   112,   113,   687,   114,   124,
     125,   126,   127,   132,   128,   133,   180,   134,   138,   139,
     140,   141,   192,   142,   524,   525,   526,   948,   146,   527,
     147,   204,    96,    97,    98,    99,   148,   100,    96,    97,
      98,    99,   149,   100,   150,   233,   151,    96,    97,    98,
      99,   234,   100,    96,    97,    98,    99,   152,   100,   153,
     235,   154,    96,    97,    98,    99,   236,   100,   110,   111,
     112,   113,   156,   114,   158,   237,   157,   110,   111,   112,
     113,   239,   114,   110,   111,   112,   113,   159,   114,   162,
     240,   161,   110,   111,   112,   113,   241,   114,   110,   111,
     112,   113,   100,   114,   169,   242,   174,   124,   125,   126,
     127,   243,   128,   124,   125,   126,   127,   114,   128,   181,
     245,   186,   124,   125,   126,   127,   246,   128,   124,   125,
     126,   127,   128,   128,   193,   247,   198,   124,   125,   126,
     127,   248,   128,   138,   139,   140,   141,   142,   142,   205,
     249,   210,   138,   139,   140,   141,   251,   142,   138,   139,
     140,   141,   211,   142,   212,   252,   213,   138,   139,   140,
     141,   253,   142,   138,   139,   140,   141,   214,   142,   215,
     254,   216,   622,   623,   624,   625,   255,   626,   622,   623,
     624,   625,   218,   626,   217,   680,   219,   622,   623,   624,
     625,   736,   626,   622,   623,   624,   625,   222,   626,   227,
     737,   223,   622,   623,   624,   625,   738,   626,   622,   623,
     624,   625,   229,   626,   238,   739,   230,   605,   606,   607,
     232,   740,   608,   847,   814,   815,   624,   625,   816,   626,
     666,   244,  1026,  1020,  1021,  1015,   848,  1022,   818,  1017,
    1018,   598,    71,    72,    73,  1027,  1030,  1031,   110,   111,
     112,   113,   250,   114,   256,   115,   124,   125,   126,   127,
     257,   128,   258,   129,   138,   139,   140,   141,   259,   142,
     261,   143,   622,   623,   624,   625,   260,   626,   590,   525,
     526,   385,   386,   527,   605,   606,   607,   415,   416,   608,
     813,   814,   815,   415,   450,   816,  1019,  1020,  1021,   822,
     879,  1022,   455,   598,   815,   884,   901,   816,   598,   598,
     990,  1021,   262,   598,  1022,    74,    75,    76,   263,   389,
     390,   391,   631,   632,   633,   415,   649,   415,   718,   872,
     873,   961,   962,   963,   964,   965,   966,   983,   984,  1007,
    1008,    26,    28,   264,   269,   271,   272,   276,   281,   283,
     284,   285,   287,   286,   288,   289,   290,   293,   302,   304,
     305,   310,   303,   306,   311,   315,   314,   312,   328,   313,
     316,   317,   329,   335,   336,   337,   339,   341,   342,   343,
     346,   344,   327,   354,   358,   359,   361,   362,   363,   364,
     365,   366,   379,   393,   373,   380,   382,   399,   400,   401,
     402,   405,   413,   420,   422,   403,   424,   406,   407,   409,
     441,   410,   411,   439,   412,   418,   444,   446,   445,   438,
     443,   419,   440,   447,   451,   460,   452,   498,   472,   458,
     461,   462,   453,   474,   459,   494,   493,   464,   463,   465,
     466,   499,   503,   505,   507,   523,   475,   520,   534,   532,
     508,   497,   500,   479,   491,   501,   492,   509,   510,   511,
     502,   504,   506,   512,   518,   530,   535,   538,   529,   537,
     569,   577,   570,   589,   571,   527,   533,   572,   573,   574,
     575,   581,   591,   595,   604,   584,   611,   621,   630,   592,
     648,   652,   593,   609,   594,   612,   653,   613,   656,   657,
     658,   596,   614,   660,   671,   615,   662,   616,   672,   674,
     686,   722,   617,   627,   618,   659,   661,   634,   719,   723,
     727,   728,   730,   646,   647,   733,   378,   488,   741,   742,
     743,   750,   764,   651,   767,   655,   768,   663,   608,   626,
     667,   771,   681,   692,   693,   688,   694,   769,   695,   735,
     696,   706,   707,   708,   721,   709,   710,   717,   772,   781,
     774,   785,   799,   724,   720,   726,   725,   633,   786,   802,
     746,   765,   729,   747,   770,   756,   757,   805,   806,   808,
     819,   835,   758,   759,   823,   824,   825,   760,   826,   827,
    -188,   842,   798,   853,   773,   766,   854,   898,   868,   775,
     869,   782,   783,   843,   863,   801,   807,   864,   865,   803,
     804,  -187,   845,   866,   867,   875,   876,   878,   880,   882,
     834,   883,   817,   895,   905,   906,   926,   833,   844,   836,
     928,   929,   930,   816,   933,   951,   979,   858,   859,   855,
     856,   860,   861,   862,   871,   978,   900,   890,   902,   955,
     980,   903,   881,   891,   981,   957,   892,   987,   893,   989,
     894,   985,   896,   958,   994,   897,   995,   997,  1000,   998,
     899,   904,   956,   907,  1005,   919,  1032,   922,  1009,  1033,
    1035,  1040,   927,   931,   923,   932,   934,  1041,   935,   936,
     937,   938,  1042,  1045,   960,   478,   918,   516,   959,   988,
      87,   967,   968,   969,   970,   971,   576,   763,   665,   972,
     601,   787,   977,   991,   921,   949,     0,   449,   992,   996,
       0,   993,     0,     0,     0,     0,  1010,  1022,  1001,  1002,
       0,     0,     0,  1006,     0,  1039,     0,  1036,     0,  1044,
       0,     0,  1043,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     266,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   295
};

static const yytype_int16 yycheck[] =
{
     312,   449,   316,   265,   363,   364,   600,   402,   505,   384,
     403,     3,     4,     4,     4,     3,     4,     4,     4,     4,
       3,     4,     3,    84,     4,   337,     3,     4,     4,   341,
       4,    85,   346,     3,     4,    89,    80,     3,     4,    83,
     399,   400,   401,     4,     3,     4,     4,    78,    79,    80,
       4,    84,    83,     3,     4,    48,    49,     3,     4,     4,
      91,    84,     3,     4,    26,    84,     4,   564,   565,    57,
      58,    59,    60,    61,    57,    58,    59,    60,    61,     3,
       4,    67,     3,     4,    84,     8,     4,    79,    84,    65,
      66,    79,    15,    85,    85,    87,    79,    89,    85,    87,
      85,    89,    79,    84,    87,    85,    89,    69,    85,    79,
      87,    85,    89,    79,     3,     4,     4,    87,     0,    89,
      79,    87,     5,    89,    85,     4,    85,    85,    87,    79,
     392,    85,     4,    79,   631,   632,   633,    87,    79,    89,
      85,    87,   517,    89,     3,     4,    87,    85,    89,     3,
       4,   745,    26,     3,     4,    79,     3,     4,    79,    88,
       4,   475,   610,    87,    84,    89,    87,    85,    89,    57,
      58,    59,    60,    61,     4,     3,     4,     4,    57,    58,
      59,    60,    61,     4,     3,    57,    58,    59,    60,    61,
      79,    79,     3,     4,     4,    69,    85,    63,    87,    87,
      79,    89,     3,     4,   518,   580,     8,    79,    87,   734,
      89,     3,     4,    15,     4,    87,   810,    89,     3,     4,
      79,     3,     4,     3,     4,    79,    85,   821,    87,    79,
       4,    85,    79,    87,     4,    85,    85,    87,    85,    88,
      87,    85,     3,     4,     3,     4,     3,     4,     4,   279,
     844,    79,     3,   778,   779,    85,     3,    85,    85,    87,
       3,    85,     4,     4,    85,    84,    86,    87,    79,     4,
       4,     4,     4,     4,    85,    85,    87,     3,    79,     3,
     310,     4,     4,    27,    85,    29,    87,    79,     4,   814,
     815,     4,     4,    85,    79,    87,     4,    79,     4,    79,
      85,     4,    87,    85,    50,    87,    52,    87,    85,   339,
     161,    85,    87,    85,   709,    85,    88,   710,    79,    28,
      79,    30,    79,    88,    85,    90,    87,    87,    87,    85,
      87,     3,    80,    84,   446,    83,     3,    84,     8,     9,
      10,    84,    87,    85,    85,    15,    84,   706,   707,   708,
      85,    85,    85,    85,    85,    66,    67,    87,    84,    88,
      84,    90,    85,    85,   759,    85,   960,   760,    88,    85,
     482,   483,    85,    85,   486,    55,    56,    85,   229,    85,
      87,   232,    85,     8,     9,    10,   541,   542,    65,    66,
      15,   501,   502,    88,   504,    90,   506,   756,   757,   758,
      87,    80,    81,   665,    83,   298,    87,   300,    78,    79,
      80,    45,    46,   525,   526,    72,    73,    74,    75,   130,
     131,   132,   133,   134,   557,   558,    87,   138,   139,   140,
     141,   324,   325,   326,    86,    87,   116,   117,   118,   119,
     120,    35,    36,    87,   124,   125,   126,   127,    50,    51,
     605,   606,   607,    78,    79,    80,     3,     4,    83,   569,
     570,   571,   572,   573,    80,    81,    87,    83,   102,   103,
     104,   105,   106,    16,    17,    18,   110,   111,   112,   113,
       5,   614,   615,   616,   617,   618,    16,    17,    18,   622,
     623,   624,   625,     3,    88,    89,    90,    91,    92,    24,
      25,   949,    96,    97,    98,    99,    27,    28,    29,    30,
       3,     4,    33,    16,    17,    18,    37,    38,    87,    62,
      16,    17,    18,    80,    81,    68,    83,    87,    53,    54,
      55,    56,    62,    87,    78,    79,    80,    81,    68,    83,
      80,    85,    63,    83,     3,     4,   858,   859,   860,   861,
      78,    79,    80,    81,    87,    83,    78,    79,    80,    62,
      87,    83,     3,    91,    85,    68,    62,    80,    81,    91,
      83,    87,    68,    87,    78,    79,    80,     4,   890,   891,
     892,   893,   896,    78,    79,    80,    81,    91,    83,    78,
      79,    80,    81,    87,    83,    87,    91,    87,    78,    79,
      80,    81,    91,    83,    78,    79,    80,   919,     4,    83,
       3,    91,    78,    79,    80,    81,     4,    83,    78,    79,
      80,    81,     4,    83,     3,    91,    84,    78,    79,    80,
      81,    91,    83,    78,    79,    80,    81,    84,    83,     3,
      91,    50,    78,    79,    80,    81,    91,    83,    78,    79,
      80,    81,    84,    83,    26,    91,    84,    78,    79,    80,
      81,    91,    83,    78,    79,    80,    81,    50,    83,    26,
      91,    84,    78,    79,    80,    81,    91,    83,    78,    79,
      80,    81,    83,    83,    88,    91,     3,    78,    79,    80,
      81,    91,    83,    78,    79,    80,    81,    83,    83,    88,
      91,     3,    78,    79,    80,    81,    91,    83,    78,    79,
      80,    81,    83,    83,    88,    91,     3,    78,    79,    80,
      81,    91,    83,    78,    79,    80,    81,    83,    83,    88,
      91,     3,    78,    79,    80,    81,    91,    83,    78,    79,
      80,    81,    63,    83,    88,    91,    88,    78,    79,    80,
      81,    91,    83,    78,    79,    80,    81,     4,    83,    88,
      91,    88,    78,    79,    80,    81,    91,    83,    78,    79,
      80,    81,    84,    83,     4,    91,     4,    78,    79,    80,
      81,    91,    83,    78,    79,    80,    81,    23,    83,     4,
      91,    84,    78,    79,    80,    81,    91,    83,    78,    79,
      80,    81,    84,    83,     3,    91,    84,    78,    79,    80,
      84,    91,    83,    78,    79,    80,    80,    81,    83,    83,
      91,     3,    78,    79,    80,  1009,    91,    83,    85,  1013,
    1014,    88,    34,    35,    36,    91,  1020,  1021,    78,    79,
      80,    81,     3,    83,     3,    85,    78,    79,    80,    81,
       3,    83,     4,    85,    78,    79,    80,    81,     4,    83,
       4,    85,    78,    79,    80,    81,    88,    83,    78,    79,
      80,     3,     4,    83,    78,    79,    80,    84,    85,    83,
      78,    79,    80,    84,    85,    83,    78,    79,    80,    85,
      85,    83,    88,    88,    80,    85,    85,    83,    88,    88,
      85,    80,     4,    88,    83,    34,    35,    36,    88,    37,
      38,    39,    78,    79,    80,    84,    85,    84,    85,    76,
      77,     3,     4,     3,     4,     3,     4,     3,     4,    20,
      21,    12,    13,    63,     4,     3,    84,     8,    47,    90,
      90,    90,    85,    90,     4,     4,     3,     3,    40,    85,
       3,     3,    31,     8,    47,    88,    23,    84,     3,    85,
      84,    51,     3,     3,    85,    84,     3,    84,     4,    89,
      84,     6,    83,    89,     3,    33,    88,    84,    84,    84,
      84,    84,     3,    32,    85,     4,    84,    84,    84,    84,
      84,    41,     4,     4,     3,    84,     3,    85,    85,    84,
       4,    52,    88,    41,    88,    88,     3,     8,     4,    84,
      84,    88,    85,     4,     3,    14,     4,    40,     4,    85,
      14,    14,    94,    84,    93,     4,    89,    14,    93,    93,
      14,     8,     8,     8,    14,     3,    84,    91,    63,     3,
      14,    84,    84,    90,    90,    84,    90,    14,    14,    14,
      84,    84,    84,    26,    84,    15,     4,     3,    89,    31,
      84,    42,    84,     3,    84,    83,    94,    84,    84,    84,
      84,    79,     3,     3,     3,    88,     4,     3,     3,    89,
      42,     3,    88,    85,    89,    85,     3,    85,     3,     3,
       3,    94,    87,     3,     3,    87,     4,    87,     4,    15,
       3,     3,    87,    85,    87,    64,    63,    85,     4,     3,
       3,     3,    32,    84,    84,     3,   342,     4,     3,     3,
       3,     3,     3,    92,     3,    88,     3,    88,    83,    83,
      88,     3,    88,    85,    85,    88,    85,    23,    85,    89,
      85,    84,    84,    84,    88,    84,    84,    84,     3,     3,
      28,    31,    14,    88,    90,    88,    90,    80,    31,     3,
      84,    88,    94,    84,    64,    84,    84,     4,     3,    32,
       3,    14,    84,    84,    43,    43,    43,    84,    43,    43,
       4,    23,    79,     3,    88,    90,     3,   871,     4,    90,
       4,    90,    90,    63,    43,    90,    85,    43,    43,    90,
      90,     4,    63,    43,    43,     4,     3,     3,     3,    32,
      84,    32,    88,    44,     3,     3,     3,    92,    84,    90,
       3,     3,     3,    83,    32,     4,     3,    84,    84,    88,
      88,    84,    84,    84,    84,     4,    64,    84,    64,    23,
       4,    63,    90,    84,     3,    23,    84,     4,    84,    23,
      84,    45,    84,    64,     3,    85,     4,     3,    23,     4,
      89,    88,    63,    88,     4,    84,     3,    85,    11,     4,
       3,     3,    85,    85,    89,    85,    85,     4,    85,    85,
      85,    85,     3,     3,    84,   443,   894,   474,    90,    64,
      28,    89,    85,    85,    85,    85,   514,   717,   601,    85,
     538,   749,   953,    88,   896,   920,    -1,   415,    88,    85,
      -1,    84,    -1,    -1,    -1,    -1,    84,    83,    90,    90,
      -1,    -1,    -1,    93,    -1,    88,    -1,    94,    -1,    90,
      -1,    -1,    94,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     220,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   269
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     5,    24,    25,    53,    54,    55,    56,    96,   119,
       4,   120,    84,    84,    84,    84,    84,    84,     0,     4,
      65,    66,    98,   133,    88,   108,   119,   109,   119,     4,
      57,    58,    59,    60,    61,    79,    87,    89,   149,     4,
      57,    58,    59,    60,    61,    79,    87,    89,   150,     4,
      57,    58,    59,    60,    61,    79,    87,    89,   151,     3,
       4,    57,    58,    59,    60,    61,    79,    87,    89,   153,
      84,    34,    35,    36,    34,    35,    36,     4,    67,    63,
       4,    85,    26,    69,   121,    85,    26,   121,    87,    87,
      87,    87,    87,   149,   149,     3,    78,    79,    80,    81,
      83,    85,    87,    87,    87,    87,    87,   150,   150,     3,
      78,    79,    80,    81,    83,    85,    87,    87,    87,    87,
      87,   151,   151,     3,    78,    79,    80,    81,    83,    85,
      87,    87,    87,    87,    87,   153,   153,     3,    78,    79,
      80,    81,    83,    85,    84,     4,     4,     3,     4,     4,
       3,    84,    84,     3,    50,   103,    84,    84,    26,    50,
     104,    84,    26,   149,   149,   149,   149,   149,    91,    88,
     149,   149,   149,   149,     3,   150,   150,   150,   150,   150,
      91,    88,   150,   150,   150,   150,     3,   151,   151,   151,
     151,   151,    91,    88,   151,   151,   151,   151,     3,   153,
     153,   153,   153,   153,    91,    88,   153,   153,   153,   153,
       3,    63,    88,    88,     4,    88,    88,     4,    84,     4,
     100,   101,    23,    84,    27,    29,   123,     4,   122,    84,
      84,   123,    84,    91,    91,    91,    91,    91,     3,    91,
      91,    91,    91,    91,     3,    91,    91,    91,    91,    91,
       3,    91,    91,    91,    91,    91,     3,     3,     4,     4,
      88,     4,     4,    88,    63,    85,   101,    86,    87,     4,
     106,     3,    84,    48,    49,   102,     8,     4,    85,   123,
     105,    47,   123,    90,    90,    90,    90,    85,     4,     4,
       3,   133,    86,     3,    51,   103,     3,     4,    79,    85,
      87,   145,    40,    31,    85,     3,     8,   102,     4,    85,
       3,    47,    84,    85,    23,    88,    84,    51,   145,   145,
       8,     9,    10,    15,    78,    79,    80,    83,     3,     3,
      16,    17,    18,    62,    68,     3,    85,    84,   102,     3,
     106,    84,     4,    89,     6,   135,    84,    91,     3,     4,
       3,     4,     3,     4,    89,   145,   145,   145,     3,    33,
     124,    88,    84,    84,    84,    84,    84,    16,    17,    18,
      62,    68,   106,    85,   102,    85,   106,    52,   104,     3,
       4,   136,    84,    97,   135,     3,     4,     3,    84,    37,
      38,    39,   127,    32,   140,   140,   140,   141,   142,    84,
      84,    84,    84,    84,    85,    41,    85,    85,    85,    84,
      52,    88,    88,     4,   147,    84,    85,    97,    88,    88,
       4,   125,     3,   133,     3,     4,    85,    85,    85,     4,
      85,     4,    85,   140,   140,   140,   141,   142,    84,    41,
      85,     4,    99,    84,     3,     4,     8,     4,   138,   147,
      85,     3,     4,    94,    85,    88,    28,    30,    85,    93,
      14,    14,    14,    93,    14,    93,    14,    85,    85,    85,
      85,    85,     4,   111,    84,    84,     4,    85,    99,    90,
       3,     4,    79,    87,    89,   148,     8,    15,     4,    85,
     138,    90,    90,    89,     4,     3,    84,    84,    40,     8,
      84,    84,    84,     8,    84,     8,    84,    14,    14,    14,
      14,    14,    26,   112,     4,    85,   111,   135,    84,    85,
      91,   148,   148,     3,    78,    79,    80,    83,   148,    89,
      15,    85,     3,    94,    63,     4,   126,    31,     3,     3,
       4,    79,    87,    89,   144,   134,   135,   139,   134,   134,
       3,     4,    57,    58,    59,    60,    61,    79,    87,    89,
     152,   134,     3,     4,    79,    87,    89,   154,   134,    84,
      84,    84,    84,    84,    84,    84,   112,    42,    85,    97,
     135,    79,    78,    91,    88,    89,   148,   148,   148,     3,
      78,     3,    89,    88,    89,     3,    94,    85,    88,     3,
      84,   124,   144,   144,     3,    78,    79,    80,    83,    85,
     137,     4,    85,    85,    87,    87,    87,    87,    87,   152,
     152,     3,    78,    79,    80,    81,    83,    85,   154,   154,
       3,    78,    79,    80,    85,   134,   134,   134,   134,   134,
     113,    16,    17,    18,    62,    68,    84,    84,    42,    85,
      97,    92,     3,     3,    89,    88,     3,     3,     3,    64,
       3,    63,     4,    88,   126,   127,    91,    88,   144,   144,
     144,     3,     4,   138,    15,   152,   152,   152,   152,   152,
      91,    88,   152,   152,   152,   152,     3,    91,    88,   154,
     154,   154,    85,    85,    85,    85,    85,    27,    28,    29,
      30,    33,    37,    38,    63,    85,    84,    84,    84,    84,
      84,    16,    17,    18,    62,    68,   114,    84,    85,     4,
      90,    88,     3,     3,    88,    90,    88,     3,     3,    94,
      32,    85,   133,     3,     8,    89,    91,    91,    91,    91,
      91,     3,     3,     3,     3,    84,    84,    84,     3,    84,
       3,   140,   140,   140,   141,   142,    84,    84,    84,    84,
      84,     4,    85,   114,     3,    88,    90,     3,     3,    23,
      64,     3,     3,    88,    28,    90,     3,     4,    79,    87,
     143,     3,    90,    90,   126,    31,    31,   125,    85,    85,
      85,    85,    85,   140,   140,   140,   141,   142,    79,    14,
      85,    90,     3,    90,    90,     4,     3,    85,    32,     3,
      84,   143,   143,    78,    79,    80,    83,    88,    85,     3,
       3,    84,    85,    43,    43,    43,    43,    43,    85,    85,
      85,    85,    85,    92,    84,    14,    90,    72,    73,    74,
      75,   128,    23,    63,    84,    63,   126,    78,    91,    89,
     143,   143,   143,     3,     3,    88,    88,   126,    84,    84,
      84,    84,    84,    43,    43,    43,    43,    43,     4,     4,
     110,    84,    76,    77,   129,     4,     3,   126,     3,    85,
       3,    90,    32,    32,    85,   106,   106,   106,   106,   107,
      84,    84,    84,    84,    84,    44,    84,    85,   110,    89,
      64,    85,    64,    63,    88,     3,     3,    88,    85,    85,
      85,    85,    85,   154,   106,   106,   106,   106,   107,    84,
     135,   139,    85,    89,   130,   131,     3,    85,     3,     3,
       3,    85,    85,    32,    85,    85,    85,    85,    85,     8,
       9,    10,    15,    85,    85,    85,    85,    85,   106,   137,
      85,     4,   132,    88,    90,    23,    63,    23,    64,    90,
      84,     3,     4,     3,     4,     3,     4,    89,    85,    85,
      85,    85,    85,    85,   138,    88,    90,   131,     4,     3,
       4,     3,   126,     3,     4,    45,    85,     4,    64,    23,
      85,    88,    88,    84,     3,     4,    85,     3,     4,   115,
      23,    90,    90,     4,    85,     4,    93,    20,    21,    11,
      84,     3,     4,    79,    87,   146,   116,   146,   146,    78,
      79,    80,    83,    85,    89,   117,    78,    91,    89,   146,
     146,   146,     3,     4,   118,     3,    94,    88,    90,    88,
       3,     4,     3,    94,    90,     3
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    95,    96,    96,    96,    96,    96,    96,    96,    96,
      96,    96,    96,    96,    97,    97,    98,    98,    99,    99,
     100,   100,   101,   101,   102,   102,   103,   104,   105,   105,
     106,   106,   106,   106,   106,   106,   106,   106,   106,   107,
     107,   107,   107,   107,   107,   107,   107,   107,   108,   108,
     108,   108,   108,   108,   108,   108,   108,   108,   109,   109,
     110,   110,   111,   111,   111,   111,   111,   111,   111,   111,
     111,   111,   112,   112,   113,   113,   113,   113,   113,   113,
     113,   113,   113,   113,   113,   113,   114,   114,   114,   115,
     115,   115,   116,   116,   117,   118,   118,   119,   120,   120,
     121,   122,   122,   123,   123,   123,   123,   123,   123,   124,
     124,   125,   125,   126,   126,   127,   127,   127,   128,   128,
     128,   128,   128,   129,   129,   129,   130,   130,   131,   132,
     132,   133,   133,   133,   133,   133,   133,   134,   134,   135,
     136,   136,   137,   137,   138,   138,   139,   139,   140,   140,
     141,   141,   142,   142,   143,   143,   143,   143,   143,   143,
     143,   143,   144,   144,   144,   144,   144,   144,   144,   144,
     144,   145,   145,   145,   145,   145,   145,   145,   145,   146,
     146,   146,   146,   146,   146,   146,   146,   147,   147,   148,
     148,   148,   148,   148,   148,   148,   148,   148,   149,   149,
     149,   149,   149,   149,   149,   149,   149,   149,   149,   149,
     149,   149,   150,   150,   150,   150,   150,   150,   150,   150,
     150,   150,   150,   150,   150,   150,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   151,   151,   151,   151,   151,
     152,   152,   152,   152,   152,   152,   152,   152,   152,   152,
     152,   152,   152,   152,   152,   153,   153,   153,   153,   153,
     153,   153,   153,   153,   153,   153,   153,   153,   153,   153,
     154,   154,   154,   154,   154,   154,   154,   154
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     4,     5,     4,     5,    12,    11,    14,    13,
       4,     4,     4,     4,     5,     4,    11,    10,     6,     5,
       3,     2,    13,     1,     1,     1,     4,     4,     5,     0,
       4,     4,     4,     8,     4,     4,     4,     8,     0,     4,
       4,     4,     8,     4,     4,     4,     8,     0,    14,    14,
      14,    14,    14,    15,    15,    15,    15,    15,    20,    21,
       4,     6,    13,    12,    13,    12,    13,    12,    13,    12,
      13,    12,     4,     0,     3,     3,     5,     2,     2,     3,
       3,     9,     9,     5,    13,     0,    17,    14,     0,    11,
       5,     0,     2,     0,     3,     5,     3,     2,     3,     1,
       4,     4,     3,    17,    21,    17,    25,    21,    23,     2,
       4,     9,     7,     5,     3,     1,     1,     2,     1,     1,
       1,     1,     0,     4,     1,     0,     3,     1,     3,     3,
       1,     5,     5,     6,     5,     5,     6,     3,     1,     2,
       3,     1,    10,     0,     8,     7,     8,     0,     5,     0,
       5,     0,     5,     0,     3,     3,     3,     3,     3,     2,
       1,     1,     3,     3,     3,     3,     3,     2,     1,     5,
       1,     3,     3,     3,     3,     3,     2,     1,     1,     3,
       3,     3,     3,     3,     2,     1,     1,    10,     9,     3,
       3,     3,     3,     3,     2,     1,     5,     1,     3,     3,
       3,     3,     3,     4,     4,     4,     4,     4,     3,     2,
       1,     5,     3,     3,     3,     3,     3,     4,     4,     4,
       4,     4,     3,     2,     1,     5,     3,     3,     3,     3,
       3,     4,     4,     4,     4,     4,     3,     2,     1,     5,
       3,     3,     3,     3,     3,     4,     4,     4,     4,     4,
       3,     2,     1,     5,     1,     3,     3,     3,     3,     3,
       4,     4,     4,     4,     4,     3,     2,     1,     5,     1,
       3,     3,     3,     3,     2,     1,     1,     5
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 118 "modelParser.y" /* yacc.c:1646  */
    {
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	logger.log(sbuilder() << "print: " << continuousProblem.bPrint);
	continuousProblem.run();
	end = clock();
	printf("%ld flowpipes computed4.\n", continuousProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	
	continuousProblem.bSafetyChecking = false;

	printf("Preparing for plotting and dumping...\n");
	continuousProblem.composition();
	printf("Done.\n");

	continuousProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, continuousProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	continuousProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);
}
#line 2085 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 3:
#line 167 "modelParser.y" /* yacc.c:1646  */
    {
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	continuousProblem.run();
	end = clock();
	printf("%ld flowpipes computed1.\n", continuousProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	continuousProblem.bSafetyChecking = true;

	printf("Preparing for plotting and dumping...\n");
	continuousProblem.composition();
	printf("Done.\n");

	continuousProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, continuousProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	continuousProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = continuousProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
#line 2158 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 4:
#line 237 "modelParser.y" /* yacc.c:1646  */
    {
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	hybridProblem.run();
	end = clock();
	printf("%ld flowpipes computed2.\n", hybridProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	hybridProblem.bSafetyChecking = false;

	printf("Preparing for plotting and dumping...\n");
	printf("Done.\n");

	hybridProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, hybridProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}

	printf("Dumping the Taylor model flowpipes...\n");
	hybridProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);
}
#line 2208 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 5:
#line 284 "modelParser.y" /* yacc.c:1646  */
    {
	int mkres = mkdir(outputDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	clock_t begin, end;
	begin = clock();
	hybridProblem.run();
	end = clock();
	printf("%ld flowpipes computed3.\n", hybridProblem.numOfFlowpipes());
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	hybridProblem.bSafetyChecking = true;

	printf("Preparing for plotting and dumping...\n");
	printf("Done.\n");

	hybridProblem.plot_2D();

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s%s.flow", outputDir, hybridProblem.outputFileName);
	FILE *fpDumping = fopen(filename, "w");

	if(fpDumping == NULL)
	{
		printf("Can not create the dumping file.\n");
		exit(1);
	}
	
	printf("Dumping the Taylor model flowpipes...\n");
	hybridProblem.dump(fpDumping);
	printf("Done.\n");

	fclose(fpDumping);

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = hybridProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
#line 2280 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 6:
#line 353 "modelParser.y" /* yacc.c:1646  */
    {
	clock_t begin, end;

	if((yyvsp[-8].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-(yyvsp[-8].dblVal),(yyvsp[-8].dblVal));
	continuousProblem.cutoff_threshold = I;

	strcpy(continuousProblem.outputFileName, (yyvsp[-6].identifier)->c_str());

	continuousProblem.plot_2D();

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = continuousProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}
}
#line 2323 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 7:
#line 393 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-7].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-(yyvsp[-7].dblVal),(yyvsp[-7].dblVal));
	continuousProblem.cutoff_threshold = I;
	
	strcpy(continuousProblem.outputFileName, (yyvsp[-5].identifier)->c_str());

	continuousProblem.plot_2D();
}
#line 2342 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 409 "modelParser.y" /* yacc.c:1646  */
    {
	clock_t begin, end;
	strcpy(hybridProblem.outputFileName, (yyvsp[-5].identifier)->c_str());
	generateNodeSeq(hybridProblem.traceNodes, hybridProblem.traceTree);

	hybridProblem.plot_2D();

	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();
	checkingResult = hybridProblem.safetyChecking();
	end = clock();
	printf("Done.\n");
	printf("time cost for safety checking: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("Result: ");

	switch(checkingResult)
	{
	case UNSAFE:
		printf("UNSAFE\n");
		break;
	case SAFE:
		printf("SAFE\n");
		break;
	case UNKNOWN:
		printf("UNKNOWN\n");
		break;
	}

	delete (yyvsp[-5].identifier);
}
#line 2378 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 442 "modelParser.y" /* yacc.c:1646  */
    {
	strcpy(hybridProblem.outputFileName, (yyvsp[-4].identifier)->c_str());
	generateNodeSeq(hybridProblem.traceNodes, hybridProblem.traceTree);
	hybridProblem.plot_2D();

	delete (yyvsp[-4].identifier);
}
#line 2390 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 451 "modelParser.y" /* yacc.c:1646  */
    {
	(yyvsp[-1].ptm)->getExpansion(parseResult.expansion);
	parseResult.remainder = (yyvsp[-1].ptm)->getRemainder();
	delete (yyvsp[-1].ptm);
}
#line 2400 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 458 "modelParser.y" /* yacc.c:1646  */
    {
	parseResult.remainder = (*(yyvsp[-1].pint));
	delete (yyvsp[-1].pint);
}
#line 2409 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 12:
#line 464 "modelParser.y" /* yacc.c:1646  */
    {
	parseResult.expansion = (*(yyvsp[-1].poly));
	delete (yyvsp[-1].poly);
}
#line 2418 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 13:
#line 470 "modelParser.y" /* yacc.c:1646  */
    {
	parseResult.strExpansion = (*(yyvsp[-1].identifier));
	delete (yyvsp[-1].identifier);
}
#line 2427 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 14:
#line 477 "modelParser.y" /* yacc.c:1646  */
    {
	continuousProblem.flowpipesCompo.push_back(*(yyvsp[-2].tmVec));
	continuousProblem.domains.push_back(*(yyvsp[-1].intVec));

	delete (yyvsp[-2].tmVec);
	delete (yyvsp[-1].intVec);
}
#line 2439 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 15:
#line 486 "modelParser.y" /* yacc.c:1646  */
    {
	continuousProblem.flowpipesCompo.push_back(*(yyvsp[-2].tmVec));
	continuousProblem.domains.push_back(*(yyvsp[-1].intVec));

	delete (yyvsp[-2].tmVec);
	delete (yyvsp[-1].intVec);
}
#line 2451 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 16:
#line 496 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModelVec tmvDummy;

	if((yyvsp[-5].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-(yyvsp[-5].dblVal),(yyvsp[-5].dblVal));
	mode_local_setting.cutoff_threshold = I;
	
	hybridProblem.declareMode(*(yyvsp[-9].identifier), tmvDummy, *(yyvsp[-2].vecConstraints), 0, mode_local_setting);

	delete (yyvsp[-9].identifier);
	delete (yyvsp[-2].vecConstraints);
}
#line 2473 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 17:
#line 515 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModelVec tmvDummy;

	if((yyvsp[-5].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-(yyvsp[-5].dblVal),(yyvsp[-5].dblVal));
	mode_local_setting.cutoff_threshold = I;

	hybridProblem.declareMode(*(yyvsp[-9].identifier), tmvDummy, *(yyvsp[-2].vecConstraints), 0, mode_local_setting);

	delete (yyvsp[-9].identifier);
	delete (yyvsp[-2].vecConstraints);
}
#line 2495 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 18:
#line 535 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForMode(*(yyvsp[-4].identifier));
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-4].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.modeIDs.push_back(id);
	hybridProblem.flowpipesCompo.push_back(continuousProblem.flowpipesCompo);
	hybridProblem.domains.push_back(continuousProblem.domains);

	continuousProblem.flowpipesCompo.clear();
	continuousProblem.domains.clear();
	continuousProblem.tmVarTab.clear();
	continuousProblem.tmVarNames.clear();
}
#line 2519 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 19:
#line 556 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForMode(*(yyvsp[-4].identifier));
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-4].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.modeIDs.push_back(id);
	hybridProblem.flowpipesCompo.push_back(continuousProblem.flowpipesCompo);
	hybridProblem.domains.push_back(continuousProblem.domains);

	continuousProblem.flowpipesCompo.clear();
	continuousProblem.domains.clear();
	continuousProblem.tmVarTab.clear();
	continuousProblem.tmVarNames.clear();
}
#line 2543 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 20:
#line 578 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 2550 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 582 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 2557 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 22:
#line 587 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForMode(*(yyvsp[0].identifier));
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	list<TreeNode *>::iterator iter = (yyval.pNode)->children.begin();
	bool found = false;
	for(; iter!=(yyval.pNode)->children.end(); ++iter)
	{
		if((*iter)->jumpID == (yyvsp[-10].dblVal) && (*iter)->modeID == id)
		{
			(yyval.pNode) = *iter;
			found = true;
			break;
		}
	}

	if(!found)
	{
		Interval I((yyvsp[-7].dblVal), (yyvsp[-5].dblVal));
		TreeNode *tmp = new TreeNode((int)(yyvsp[-10].dblVal), id, I);
		tmp->parent = (yyval.pNode);
		(yyval.pNode)->children.push_back(tmp);
		(yyval.pNode) = tmp;
	}

	delete (yyvsp[0].identifier);
}
#line 2595 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 622 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForMode(*(yyvsp[0].identifier));
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(hybridProblem.traceTree == NULL)
	{
		Interval intZero;
		hybridProblem.traceTree = new TreeNode(0, id, intZero);
		(yyval.pNode) = hybridProblem.traceTree;
	}
	else
	{
		if(hybridProblem.traceTree->modeID == id)
		{
			(yyval.pNode) = hybridProblem.traceTree;
		}
		else
		{
			parseError("Invalid computation path.", lineNum);
			exit(1);
		}
	}

	delete (yyvsp[0].identifier);
}
#line 2631 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 656 "modelParser.y" /* yacc.c:1646  */
    {
	continuousProblem.bPrint = true;
	hybridProblem.bPrint = true;
}
#line 2640 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 25:
#line 662 "modelParser.y" /* yacc.c:1646  */
    {
	continuousProblem.bPrint = false;
	hybridProblem.bPrint = false;
}
#line 2649 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 669 "modelParser.y" /* yacc.c:1646  */
    {
	continuousProblem.unsafeSet = *(yyvsp[-1].vecConstraints);
	delete (yyvsp[-1].vecConstraints);
}
#line 2658 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 27:
#line 676 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 2665 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 681 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForMode(*(yyvsp[-3].identifier));
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.unsafeSet[id] = *(yyvsp[-1].vecConstraints);
	hybridProblem.bVecUnderCheck[id] = true;
	delete (yyvsp[-1].vecConstraints);
}
#line 2684 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 29:
#line 696 "modelParser.y" /* yacc.c:1646  */
    {
	vector<PolynomialConstraint> vecEmpty;
	for(int i=0; i<hybridProblem.modeNames.size(); ++i)
	{
		hybridProblem.unsafeSet.push_back(vecEmpty);
		hybridProblem.bVecUnderCheck.push_back(false);
	}
}
#line 2697 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 30:
#line 707 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval B((yyvsp[0].dblVal));
	PolynomialConstraint pc(*(yyvsp[-2].poly), B);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
}
#line 2716 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 723 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);

	Interval B(-(yyvsp[0].dblVal));
	PolynomialConstraint pc(*(yyvsp[-2].poly), B);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
}
#line 2738 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 32:
#line 742 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval B((yyvsp[0].dblVal));
	PolynomialConstraint pc1(*(yyvsp[-2].poly), B);
	(yyval.vecConstraints)->push_back(pc1);

	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);
	Interval mB(-(yyvsp[0].dblVal));
	PolynomialConstraint pc2(*(yyvsp[-2].poly), mB);
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-2].poly);
}
#line 2763 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 33:
#line 764 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-6].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*(yyvsp[-6].poly), (yyvsp[-1].dblVal));
	(yyval.vecConstraints)->push_back(pc1);

	Interval I(-1);
	(yyvsp[-6].poly)->mul_assign(I);
	PolynomialConstraint pc2(*(yyvsp[-6].poly), -(yyvsp[-3].dblVal));
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-6].poly);
}
#line 2785 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 34:
#line 783 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	PolynomialConstraint pc(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
	delete (yyvsp[0].identifier);
}
#line 2819 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 814 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);
	range *= I;
	PolynomialConstraint pc(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
	delete (yyvsp[0].identifier);
}
#line 2856 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 36:
#line 848 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	PolynomialConstraint pc1(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc1);

	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-2].poly);
	delete (yyvsp[0].identifier);
}
#line 2896 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 37:
#line 885 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-6].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*(yyvsp[-1].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[-1].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[-1].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*(yyvsp[-6].poly), range);
	(yyval.vecConstraints)->push_back(pc1);

	id = continuousProblem.getIDForPar(*(yyvsp[-3].identifier));

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[-3].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval I(-1);
	(yyvsp[-6].poly)->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*(yyvsp[-6].poly), range);
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-6].poly);
	delete (yyvsp[-3].identifier);
	delete (yyvsp[-1].identifier);
}
#line 2950 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 38:
#line 935 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.vecConstraints) = new vector<PolynomialConstraint>(0);
}
#line 2958 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 941 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval B((yyvsp[0].dblVal));
	PolynomialConstraint pc(*(yyvsp[-2].poly), B);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
}
#line 2977 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 40:
#line 957 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);

	Interval B(-(yyvsp[0].dblVal));
	PolynomialConstraint pc(*(yyvsp[-2].poly), B);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
}
#line 2999 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 41:
#line 976 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval B((yyvsp[0].dblVal));
	PolynomialConstraint pc1(*(yyvsp[-2].poly), B);
	(yyval.vecConstraints)->push_back(pc1);

	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);
	Interval mB(-(yyvsp[0].dblVal));
	PolynomialConstraint pc2(*(yyvsp[-2].poly), mB);
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-2].poly);
}
#line 3024 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 998 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-6].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*(yyvsp[-6].poly), (yyvsp[-1].dblVal));
	(yyval.vecConstraints)->push_back(pc1);

	Interval I(-1);
	(yyvsp[-6].poly)->mul_assign(I);
	PolynomialConstraint pc2(*(yyvsp[-6].poly), -(yyvsp[-3].dblVal));
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-6].poly);
}
#line 3046 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 43:
#line 1017 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}
	
	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	PolynomialConstraint pc(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
	delete (yyvsp[0].identifier);
}
#line 3080 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 44:
#line 1048 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);
	range *= I;
	PolynomialConstraint pc(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc);

	delete (yyvsp[-2].poly);
	delete (yyvsp[0].identifier);
}
#line 3117 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 1082 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}

	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.vecConstraints) = (yyvsp[-3].vecConstraints);
	PolynomialConstraint pc1(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc1);

	Interval I(-1);
	(yyvsp[-2].poly)->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*(yyvsp[-2].poly), range);
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-2].poly);
	delete (yyvsp[0].identifier);
}
#line 3157 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 46:
#line 1119 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-6].poly)->degree() < 1)
	{
		parseError("Invalid constraint.", lineNum);
		exit(1);
	}
	
	int id = continuousProblem.getIDForPar(*(yyvsp[-1].identifier));
	Interval range;

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[-1].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[-1].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	PolynomialConstraint pc1(*(yyvsp[-6].poly), range);
	(yyval.vecConstraints)->push_back(pc1);

	id = continuousProblem.getIDForPar(*(yyvsp[-3].identifier));

	if(id >= 0)
	{
		continuousProblem.getRangeForPar(range, *(yyvsp[-3].identifier));
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Parameter %s is not declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval I(-1);
	(yyvsp[-6].poly)->mul_assign(I);
	range *= I;
	PolynomialConstraint pc2(*(yyvsp[-6].poly), range);
	(yyval.vecConstraints)->push_back(pc2);

	delete (yyvsp[-6].poly);
	delete (yyvsp[-3].identifier);
	delete (yyvsp[-1].identifier);
}
#line 3211 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 1169 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.vecConstraints) = new vector<PolynomialConstraint>(0);
}
#line 3219 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 1175 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("CONT1");
	logger.log(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	//logger.log(sbuilder() << *$9); // const TaylorModelVec’
	//logger.log(sbuilder() << *$13); // const Flowpipe
	logger.log(sbuilder() << "ONLY_PICARD: " << ONLY_PICARD);
	
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = ONLY_PICARD;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3238 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 1191 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LOW_DEGREE;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3251 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 1201 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = HIGH_DEGREE;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3264 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 1211 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].strVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = NONPOLY_TAYLOR;

	delete (yyvsp[-5].strVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3277 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 52:
#line 1221 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LINEAR;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3290 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 1231 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = ONLY_PICARD;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3303 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 1241 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LOW_DEGREE;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3316 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 1251 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = HIGH_DEGREE;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3329 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 1261 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].strVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = NONPOLY_TAYLOR;

	delete (yyvsp[-5].strVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3342 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 57:
#line 1271 "modelParser.y" /* yacc.c:1646  */
    {
	ContinuousSystem system(*(yyvsp[-5].tmVec), *(yyvsp[-1].pFlowpipe));
	continuousProblem.system = system;
	continuousProblem.integrationScheme = LINEAR;

	delete (yyvsp[-5].tmVec);
	delete (yyvsp[-1].pFlowpipe);
}
#line 3355 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 1282 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-14].dblVal) < 0)
	{
		parseError("The maximum jump depth should be a nonnegative integer.", lineNum);
		exit(1);
	}

	hybridProblem.maxJumps = (int)(yyvsp[-14].dblVal);
}
#line 3369 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 1293 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-14].dblVal) < 0)
	{
		parseError("The maximum jump depth should be a nonnegative integer.", lineNum);
		exit(1);
	}

	hybridProblem.maxJumps = (int)(yyvsp[-14].dblVal);
}
#line 3383 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 1305 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForMode(*(yyvsp[-3].identifier));
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	Flowpipe initialSet(*(yyvsp[-1].intVec), intZero);
	hybridProblem.initialConfig(id, initialSet);

	int numVars = hybridProblem.stateVarNames.size();

	string tVar("local_t");
	hybridProblem.declareTMVar(tVar);
	continuousProblem.declareTMVar(tVar);

	char name[NAME_SIZE];

	for(int i=0; i<numVars; ++i)
	{
		sprintf(name, "%s%d", local_var_name, i+1);
		string tmVarName(name);
		hybridProblem.declareTMVar(tmVarName);
		continuousProblem.declareTMVar(tmVarName);
	}

	delete (yyvsp[-3].identifier);
	delete (yyvsp[-1].intVec);
}
#line 3421 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 1340 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForMode(*(yyvsp[-5].identifier));
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-5].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Flowpipe initialSet(*(yyvsp[-2].tmVec), *(yyvsp[-1].intVec));
	hybridProblem.initialConfig(id, initialSet);

	delete (yyvsp[-2].tmVec);
	delete (yyvsp[-1].intVec);
}
#line 3442 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 62:
#line 1359 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), ONLY_PICARD, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3455 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 63:
#line 1369 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), ONLY_PICARD, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3468 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 64:
#line 1379 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), LOW_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3481 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 65:
#line 1389 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), LOW_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3494 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 66:
#line 1399 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), HIGH_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3507 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 67:
#line 1409 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), HIGH_DEGREE, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3520 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 68:
#line 1419 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].strVec), *(yyvsp[-2].vecConstraints), NONPOLY_TAYLOR, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].strVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3533 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 69:
#line 1429 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].strVec), *(yyvsp[-2].vecConstraints), NONPOLY_TAYLOR, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].strVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3546 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 70:
#line 1439 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), LINEAR, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3559 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 71:
#line 1449 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareMode(*(yyvsp[-11].identifier), *(yyvsp[-6].tmVec), *(yyvsp[-2].vecConstraints), LINEAR, mode_local_setting);
	mode_local_setting.clear();

	delete (yyvsp[-11].identifier);
	delete (yyvsp[-6].tmVec);
	delete (yyvsp[-2].vecConstraints);
}
#line 3572 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 72:
#line 1462 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 3579 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 73:
#line 1465 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 3586 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 74:
#line 1470 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[0].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveSteps = false;
	mode_local_setting.step = (yyvsp[0].dblVal);
}
#line 3601 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 75:
#line 1482 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[0].dblVal) <= 0)
	{
		parseError("Remainder estimation should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-(yyvsp[0].dblVal), (yyvsp[0].dblVal));

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		mode_local_setting.estimation.push_back(I);
	}
}
#line 3620 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 76:
#line 1498 "modelParser.y" /* yacc.c:1646  */
    {
	mode_local_setting.estimation = *(yyvsp[-1].intVec);
}
#line 3628 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 77:
#line 1503 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("QR1.\n");
	mode_local_setting.precondition = QR_PRE;
}
#line 3637 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 78:
#line 1509 "modelParser.y" /* yacc.c:1646  */
    {
	mode_local_setting.precondition = ID_PRE;
}
#line 3645 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 79:
#line 1514 "modelParser.y" /* yacc.c:1646  */
    {
	int order = (int)(yyvsp[0].dblVal);

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveOrders = false;
	mode_local_setting.orderType = UNIFORM;
	mode_local_setting.orders.push_back(order);
	mode_local_setting.globalMaxOrder = order;
}
#line 3664 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 80:
#line 1530 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[0].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-(yyvsp[0].dblVal), (yyvsp[0].dblVal));
	mode_local_setting.cutoff_threshold = I;
}
#line 3679 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 81:
#line 1542 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-4].dblVal) <= 0 || (yyvsp[-1].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	if((yyvsp[-4].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveSteps = true;
	mode_local_setting.step = (yyvsp[-1].dblVal);
	mode_local_setting.miniStep = (yyvsp[-4].dblVal);
	mode_local_setting.bAdaptiveOrders = false;
}
#line 3702 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 82:
#line 1562 "modelParser.y" /* yacc.c:1646  */
    {
	int minOrder = (int)(yyvsp[-4].dblVal);
	int maxOrder = (int)(yyvsp[-1].dblVal);

	if(minOrder <= 0 || maxOrder <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	if(minOrder > maxOrder)
	{
		parseError("MAX order should be no smaller than MIN order.", lineNum);
		exit(1);
	}

	mode_local_setting.bAdaptiveSteps = false;
	mode_local_setting.bAdaptiveOrders = true;
	mode_local_setting.orderType = UNIFORM;
	mode_local_setting.orders.push_back(minOrder);
	mode_local_setting.maxOrders.push_back(maxOrder);
	mode_local_setting.globalMaxOrder = maxOrder;
}
#line 3730 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 83:
#line 1587 "modelParser.y" /* yacc.c:1646  */
    {
	mode_local_setting.bAdaptiveOrders = false;
	mode_local_setting.orderType = MULTI;
	mode_local_setting.orders = *(yyvsp[-1].iVec);

	for(int i=0; i<(yyvsp[-1].iVec)->size(); ++i)
	{
		if((*(yyvsp[-1].iVec))[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*(yyvsp[-1].iVec))[0];
	for(int i=1; i<(yyvsp[-1].iVec)->size(); ++i)
	{
		if(maxOrder < (*(yyvsp[-1].iVec))[i])
		{
			maxOrder = (*(yyvsp[-1].iVec))[i];
		}
	}

	mode_local_setting.globalMaxOrder = maxOrder;
}
#line 3760 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 84:
#line 1614 "modelParser.y" /* yacc.c:1646  */
    {
	mode_local_setting.bAdaptiveSteps = false;
	mode_local_setting.bAdaptiveOrders = true;
	mode_local_setting.orderType = MULTI;
	mode_local_setting.orders = *(yyvsp[-7].iVec);
	mode_local_setting.maxOrders = *(yyvsp[-2].iVec);

	if((yyvsp[-7].iVec)->size() != (yyvsp[-2].iVec)->size())
	{
		parseError("Orders are not properly specified.", lineNum);
		exit(1);
	}

	for(int i=0; i<(yyvsp[-2].iVec)->size(); ++i)
	{
		if((*(yyvsp[-7].iVec))[i] <= 0 || (*(yyvsp[-2].iVec))[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}

		if((*(yyvsp[-7].iVec))[i] > (*(yyvsp[-2].iVec))[i])
		{
			parseError("MAX order should be no smaller than MIN order.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*(yyvsp[-2].iVec))[0];
	for(int i=1; i<(yyvsp[-2].iVec)->size(); ++i)
	{
		if(maxOrder < (*(yyvsp[-2].iVec))[i])
		{
			maxOrder = (*(yyvsp[-2].iVec))[i];
		}
	}

	mode_local_setting.globalMaxOrder = maxOrder;
}
#line 3804 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 85:
#line 1654 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 3811 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 86:
#line 1659 "modelParser.y" /* yacc.c:1646  */
    {
	int startID = hybridProblem.getIDForMode(*(yyvsp[-15].identifier));
	if(startID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-15].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int endID = hybridProblem.getIDForMode(*(yyvsp[-12].identifier));
	if(endID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-12].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-1].dVecVec)->size() > 0)
	{
		hybridProblem.declareTrans(startID, endID, *(yyvsp[-9].vecConstraints), *(yyvsp[-5].resetMap), PARA_AGGREG, *(yyvsp[-1].dVecVec));
	}
	else
	{
		vector<vector<double> > emptyVec;
		hybridProblem.declareTrans(startID, endID, *(yyvsp[-9].vecConstraints), *(yyvsp[-5].resetMap), PARA_AGGREG, emptyVec);
	}
}
#line 3845 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 87:
#line 1690 "modelParser.y" /* yacc.c:1646  */
    {
	int startID = hybridProblem.getIDForMode(*(yyvsp[-12].identifier));
	if(startID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-12].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int endID = hybridProblem.getIDForMode(*(yyvsp[-9].identifier));
	if(endID < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Mode %s has not been declared.", (*(yyvsp[-9].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	vector<vector<double> > empty;
	hybridProblem.declareTrans(startID, endID, *(yyvsp[-6].vecConstraints), *(yyvsp[-2].resetMap), INTERVAL_AGGREG, empty);
}
#line 3872 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 88:
#line 1713 "modelParser.y" /* yacc.c:1646  */
    {
	hybridProblem.declareTrans();
}
#line 3880 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 89:
#line 1719 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.resetMap) = (yyvsp[-10].resetMap);

	int id = hybridProblem.getIDForStateVar(*(yyvsp[-9].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-9].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("Invalid remainder interval.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal), (yyvsp[-1].dblVal));
	TaylorModel tmTemp(*(yyvsp[-6].poly), I);
	(yyval.resetMap)->tmvReset.tms[id] = tmTemp;

	delete (yyvsp[-6].poly);
}
#line 3910 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 90:
#line 1746 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.resetMap) = (yyvsp[-4].resetMap);

	int id = hybridProblem.getIDForStateVar(*(yyvsp[-3].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*(yyvsp[0].poly), intZero);
	(yyval.resetMap)->tmvReset.tms[id] = tmTemp;

	delete (yyvsp[0].poly);
}
#line 3934 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 91:
#line 1766 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = hybridProblem.stateVarNames.size();

	Matrix coefficients_identity_reset(numVars, numVars+1);

	for(int i=0; i<numVars; ++i)
	{
		coefficients_identity_reset.set(1, i, i+1);
	}

	TaylorModelVec tmvReset(coefficients_identity_reset);

	(yyval.resetMap) = new ResetMap(tmvReset);
}
#line 3953 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 92:
#line 1783 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.dVecVec)->push_back(*(yyvsp[0].dVec));
	delete (yyvsp[0].dVec);
}
#line 3962 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 93:
#line 1788 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.dVecVec) = new vector<vector<double> >(0);
}
#line 3970 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 94:
#line 1794 "modelParser.y" /* yacc.c:1646  */
    {
	int rangeDim = (yyvsp[-1].dVec)->size();

	if(rangeDim != hybridProblem.stateVarNames.size())
	{
		parseError("The vector dimension should be equivalent to the system dimension.", lineNum);
		exit(1);
	}

	(yyval.dVec) = new vector<double>(0);

	for(int i=0; i<rangeDim; ++i)
	{
		(yyval.dVec)->push_back(0);
	}

	bool bZero = true;
	for(int i=0; i<rangeDim; ++i)
	{
		if((*(yyvsp[-1].dVec))[i] < -THRESHOLD_LOW || (*(yyvsp[-1].dVec))[i] > THRESHOLD_LOW)
		{
			if(bZero)
			{
				bZero = false;
			}
		}

		(*(yyval.dVec))[i] = (*(yyvsp[-1].dVec))[i];
	}

	if(bZero)
	{
		parseError("A template vector should not be zero.", lineNum);
		exit(1);
	}

	delete (yyvsp[-1].dVec);
}
#line 4013 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 95:
#line 1835 "modelParser.y" /* yacc.c:1646  */
    {
	int id = hybridProblem.getIDForStateVar(*(yyvsp[-2].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.dVec) = (yyvsp[-4].dVec);
	(*(yyval.dVec))[id] = (yyvsp[0].dblVal);
	delete (yyvsp[-2].identifier);
}
#line 4033 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 96:
#line 1852 "modelParser.y" /* yacc.c:1646  */
    {
	int num = hybridProblem.stateVarNames.size();
	(yyval.dVec) = new vector<double>(num);

	int id = hybridProblem.getIDForStateVar(*(yyvsp[-2].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*(yyval.dVec))[id] = (yyvsp[0].dblVal);
	delete (yyvsp[-2].identifier);
}
#line 4055 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 97:
#line 1872 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 4062 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 98:
#line 1877 "modelParser.y" /* yacc.c:1646  */
    {
	//logger.log("stateIdDeclList1");
  logger.inc();
	//logger.log(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	if(!continuousProblem.declareStateVar(*(yyvsp[0].identifier)))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	//logger.log(sbuilder() << "*$3: " << *$3);
  logger.dec();
	hybridProblem.declareStateVar(*(yyvsp[0].identifier));
	delete (yyvsp[0].identifier);
}
#line 4084 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 99:
#line 1896 "modelParser.y" /* yacc.c:1646  */
    {
	//logger.log("stateIdDeclList2");
  logger.inc();
	//logger.log(sbuilder() << "*$1: " << *$1);
	//logger.log(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	if(!continuousProblem.declareStateVar(*(yyvsp[0].identifier)))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
	
	//logger.log(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	hybridProblem.declareStateVar(*(yyvsp[0].identifier));
	//logger.log(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
  logger.dec();
	delete (yyvsp[0].identifier);
}
#line 4108 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 100:
#line 1920 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 4115 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 101:
#line 1925 "modelParser.y" /* yacc.c:1646  */
    {
	if(continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier)) != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
	else
	{
		Interval range((yyvsp[0].dblVal));

		if(!continuousProblem.declarePar(*(yyvsp[-2].identifier), range))
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "Parameter %s has already been declared.", (*(yyvsp[-2].identifier)).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		hybridProblem.declarePar(*(yyvsp[-2].identifier), range);
	}
}
#line 4143 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 102:
#line 1950 "modelParser.y" /* yacc.c:1646  */
    {
	if(continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier)) != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
	else
	{
		Interval range((yyvsp[0].dblVal));

		if(!continuousProblem.declarePar(*(yyvsp[-2].identifier), range))
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "Parameter %s has already been declared.", (*(yyvsp[-2].identifier)).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		hybridProblem.declarePar(*(yyvsp[-2].identifier), range);
	}
}
#line 4171 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 103:
#line 1978 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("settings1");
  logger.log(continuousProblem.algorithm);
	if((yyvsp[-15].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	int order = (int)(yyvsp[-8].dblVal);

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}
	logger.log(sbuilder() << "step: " << (yyvsp[-15].dblVal) << ", time: " << time << ", order: " << order);
	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = (yyvsp[-15].dblVal);
	continuousProblem.time = (yyvsp[-13].dblVal);
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(order);
	continuousProblem.globalMaxOrder = order;

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = (yyvsp[-15].dblVal);
	hybridProblem.time = (yyvsp[-13].dblVal);
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = UNIFORM;
	hybridProblem.global_setting.orders.push_back(order);
	hybridProblem.global_setting.globalMaxOrder = order;

	if((yyvsp[-6].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-(yyvsp[-6].dblVal),(yyvsp[-6].dblVal));
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)(yyvsp[-4].dblVal);

	strcpy(continuousProblem.outputFileName, (*(yyvsp[-2].identifier)).c_str());
	strcpy(hybridProblem.outputFileName, (*(yyvsp[-2].identifier)).c_str());
	logger.log(sbuilder() << "cutoff_thrs: " << (yyvsp[-6].dblVal) << ", intnumprec: " << (yyvsp[-4].dblVal) << ", outFname: " << *(yyvsp[-2].identifier));
	

	delete (yyvsp[-2].identifier);
}
#line 4228 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 104:
#line 2032 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("settings2");
	if((yyvsp[-19].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	int minOrder = (int)(yyvsp[-10].dblVal);
	int maxOrder = (int)(yyvsp[-7].dblVal);

	if(minOrder <= 0 || maxOrder <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	if(minOrder > maxOrder)
	{
		parseError("MAX order should be no smaller than MIN order.", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = (yyvsp[-19].dblVal);
	continuousProblem.time = (yyvsp[-17].dblVal);
	continuousProblem.bAdaptiveOrders = true;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(minOrder);
	continuousProblem.maxOrders.push_back(maxOrder);
	continuousProblem.globalMaxOrder = maxOrder;

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = (yyvsp[-19].dblVal);
	hybridProblem.time = (yyvsp[-17].dblVal);
	hybridProblem.global_setting.bAdaptiveOrders = true;
	hybridProblem.global_setting.orderType = UNIFORM;
	hybridProblem.global_setting.orders.push_back(minOrder);
	hybridProblem.global_setting.maxOrders.push_back(maxOrder);
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if((yyvsp[-4].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-(yyvsp[-4].dblVal),(yyvsp[-4].dblVal));
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)(yyvsp[-2].dblVal);

	strcpy(continuousProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());
	strcpy(hybridProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());

	delete (yyvsp[0].identifier);
}
#line 4291 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 105:
#line 2092 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("settings3");
	if((yyvsp[-15].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = (yyvsp[-15].dblVal);
	continuousProblem.time = (yyvsp[-13].dblVal);
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *(yyvsp[-7].iVec);

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = (yyvsp[-15].dblVal);
	hybridProblem.time = (yyvsp[-13].dblVal);
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = MULTI;
	hybridProblem.global_setting.orders = *(yyvsp[-7].iVec);

	for(int i=0; i<(yyvsp[-7].iVec)->size(); ++i)
	{
		if((*(yyvsp[-7].iVec))[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*(yyvsp[-7].iVec))[0];
	for(int i=1; i<(yyvsp[-7].iVec)->size(); ++i)
	{
		if(maxOrder < (*(yyvsp[-7].iVec))[i])
		{
			maxOrder = (*(yyvsp[-7].iVec))[i];
		}
	}

	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if((yyvsp[-4].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-(yyvsp[-4].dblVal),(yyvsp[-4].dblVal));
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)(yyvsp[-2].dblVal);

	strcpy(continuousProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());
	strcpy(hybridProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());

	delete (yyvsp[-7].iVec);
	delete (yyvsp[0].identifier);
}
#line 4357 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 106:
#line 2155 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("settings4");
	if((yyvsp[-23].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = (yyvsp[-23].dblVal);
	continuousProblem.time = (yyvsp[-21].dblVal);
	continuousProblem.bAdaptiveOrders = true;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *(yyvsp[-13].iVec);
	continuousProblem.maxOrders = *(yyvsp[-8].iVec);

	hybridProblem.global_setting.bAdaptiveSteps = false;
	hybridProblem.global_setting.step = (yyvsp[-23].dblVal);
	hybridProblem.time = (yyvsp[-21].dblVal);
	hybridProblem.global_setting.bAdaptiveOrders = true;
	hybridProblem.global_setting.orderType = MULTI;
	hybridProblem.global_setting.orders = *(yyvsp[-13].iVec);
	hybridProblem.global_setting.maxOrders = *(yyvsp[-8].iVec);

	if((yyvsp[-13].iVec)->size() != (yyvsp[-8].iVec)->size())
	{
		parseError("Orders are not properly specified.", lineNum);
		exit(1);
	}

	for(int i=0; i<(yyvsp[-8].iVec)->size(); ++i)
	{
		if((*(yyvsp[-13].iVec))[i] <= 0 || (*(yyvsp[-8].iVec))[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}

		if((*(yyvsp[-13].iVec))[i] > (*(yyvsp[-8].iVec))[i])
		{
			parseError("MAX order should be no smaller than MIN order.", lineNum);
			exit(1);
		}
	}

	int maxOrder = (*(yyvsp[-8].iVec))[0];
	for(int i=1; i<(yyvsp[-8].iVec)->size(); ++i)
	{
		if(maxOrder < (*(yyvsp[-8].iVec))[i])
		{
			maxOrder = (*(yyvsp[-8].iVec))[i];
		}
	}

	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if((yyvsp[-4].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-(yyvsp[-4].dblVal),(yyvsp[-4].dblVal));
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)(yyvsp[-2].dblVal);

	strcpy(continuousProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());
	strcpy(hybridProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());

	delete (yyvsp[-13].iVec);
	delete (yyvsp[-8].iVec);
	delete (yyvsp[0].identifier);
}
#line 4438 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 107:
#line 2233 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("settings5");
	if((yyvsp[-17].dblVal) <= 0 || (yyvsp[-14].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	if((yyvsp[-17].dblVal) > (yyvsp[-14].dblVal))
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	int order = (int)(yyvsp[-6].dblVal);

	if(order <= 0)
	{
		parseError("Orders should be larger than zero.", lineNum);
		exit(1);
	}

	continuousProblem.bAdaptiveSteps = true;
	continuousProblem.step = (yyvsp[-14].dblVal);
	continuousProblem.miniStep = (yyvsp[-17].dblVal);
	continuousProblem.time = (yyvsp[-11].dblVal);
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(order);
	continuousProblem.globalMaxOrder = order;

	hybridProblem.global_setting.bAdaptiveSteps = true;
	hybridProblem.global_setting.step = (yyvsp[-14].dblVal);
	hybridProblem.global_setting.miniStep = (yyvsp[-17].dblVal);
	hybridProblem.time = (yyvsp[-11].dblVal);
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = UNIFORM;
	hybridProblem.global_setting.orders.push_back(order);
	hybridProblem.global_setting.globalMaxOrder = order;

	if((yyvsp[-4].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-(yyvsp[-4].dblVal),(yyvsp[-4].dblVal));
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)(yyvsp[-2].dblVal);

	strcpy(continuousProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());
	strcpy(hybridProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());

	delete (yyvsp[0].identifier);
}
#line 4500 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 108:
#line 2292 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("settings6");
	if((yyvsp[-19].dblVal) <= 0 || (yyvsp[-16].dblVal) <= 0)
	{
		parseError("A time step-size should be larger than 0", lineNum);
		exit(1);
	}

	if((yyvsp[-19].dblVal) > (yyvsp[-16].dblVal))
	{
		parseError("MIN step should be no larger than MAX step.", lineNum);
		exit(1);
	}

	for(int i=0; i<(yyvsp[-7].iVec)->size(); ++i)
	{
		if((*(yyvsp[-7].iVec))[i] <= 0)
		{
			parseError("Orders should be larger than zero.", lineNum);
			exit(1);
		}
	}

	continuousProblem.bAdaptiveSteps = true;
	continuousProblem.step = (yyvsp[-16].dblVal);
	continuousProblem.miniStep = (yyvsp[-19].dblVal);
	continuousProblem.time = (yyvsp[-13].dblVal);
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = MULTI;
	continuousProblem.orders = *(yyvsp[-7].iVec);

	hybridProblem.global_setting.bAdaptiveSteps = true;
	hybridProblem.global_setting.step = (yyvsp[-16].dblVal);
	hybridProblem.global_setting.miniStep = (yyvsp[-19].dblVal);
	hybridProblem.time = (yyvsp[-13].dblVal);
	hybridProblem.global_setting.bAdaptiveOrders = false;
	hybridProblem.global_setting.orderType = MULTI;
	hybridProblem.global_setting.orders = *(yyvsp[-7].iVec);

	int maxOrder = (*(yyvsp[-7].iVec))[0];
	for(int i=1; i<(yyvsp[-7].iVec)->size(); ++i)
	{
		if(maxOrder < (*(yyvsp[-7].iVec))[i])
		{
			maxOrder = (*(yyvsp[-7].iVec))[i];
		}
	}
	
	continuousProblem.globalMaxOrder = maxOrder;
	hybridProblem.global_setting.globalMaxOrder = maxOrder;

	if((yyvsp[-4].dblVal) <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-(yyvsp[-4].dblVal),(yyvsp[-4].dblVal));
	continuousProblem.cutoff_threshold = cutoff_threshold;
	hybridProblem.global_setting.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = (int)(yyvsp[-2].dblVal);

	strcpy(continuousProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());
	strcpy(hybridProblem.outputFileName, (*(yyvsp[0].identifier)).c_str());

	delete (yyvsp[-7].iVec);
	delete (yyvsp[0].identifier);
}
#line 4574 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 109:
#line 2364 "modelParser.y" /* yacc.c:1646  */
    {
	//logger.log("remainder_estimation1");
	if((yyvsp[0].dblVal) <= 0)
	{
		parseError("Remainder estimation should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-(yyvsp[0].dblVal), (yyvsp[0].dblVal));

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		//logger.log(sbuilder() << "i: " << i << " - " << $2);
		continuousProblem.estimation.push_back(I);
		hybridProblem.global_setting.estimation.push_back(I);
	}
}
#line 4596 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 110:
#line 2383 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("remainder_estimation2");
	for(int i=0; i<(yyvsp[-1].intVec)->size(); ++i)
	{
		if((*(yyvsp[-1].intVec))[i].inf() >= (*(yyvsp[-1].intVec))[i].sup() - THRESHOLD_LOW)
		{
			parseError("Invalid remainder estimation.", lineNum);
			exit(1);
		}
	}

	continuousProblem.estimation = *(yyvsp[-1].intVec);
	hybridProblem.global_setting.estimation = *(yyvsp[-1].intVec);
	delete (yyvsp[-1].intVec);
}
#line 4616 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 111:
#line 2401 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.intVec) = (yyvsp[-8].intVec);
	int id = continuousProblem.getIDForStateVar(*(yyvsp[-6].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-6].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) >= (yyvsp[-1].dblVal))
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	(*(yyval.intVec))[id] = I;
	delete (yyvsp[-6].identifier);
}
#line 4643 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 112:
#line 2425 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size();
	(yyval.intVec) = new vector<Interval>(numVars);

	int id = continuousProblem.getIDForStateVar(*(yyvsp[-6].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-6].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) >= (yyvsp[-1].dblVal))
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	(*(yyval.intVec))[id] = I;
	delete (yyvsp[-6].identifier);
}
#line 4672 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 113:
#line 2452 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.iVec) = (yyvsp[-4].iVec);
	int id = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*(yyval.iVec))[id] = (int)(yyvsp[0].dblVal);
	delete (yyvsp[-2].identifier);
}
#line 4692 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 114:
#line 2469 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size();
	(yyval.iVec) = new vector<int>(numVars);
	for(int i=0; i<numVars; ++i)
	{
		(*(yyval.iVec))[i] = 0;
	}

	int id = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*(yyval.iVec))[id] = (int)(yyvsp[0].dblVal);
	delete (yyvsp[-2].identifier);
}
#line 4718 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 115:
#line 2493 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("QR1");
	logger.log(sbuilder() << "QR_PRE: " << QR_PRE);
	continuousProblem.precondition = QR_PRE;
	hybridProblem.global_setting.precondition = QR_PRE;
}
#line 4729 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 116:
#line 2501 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("QR2");
	continuousProblem.precondition = ID_PRE;
	hybridProblem.global_setting.precondition = ID_PRE;
}
#line 4739 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 117:
#line 2508 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("shrink");
	int everyNthStep = (yyvsp[0].dblVal);
	//logger.log(everyNthStep);
	continuousProblem.precondition = SHRINK_WRAPPING;
	continuousProblem.sw_step = everyNthStep;
}
#line 4751 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 118:
#line 2518 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ALG_SIMPLE_IMPL");
	if(continuousProblem.precondition == SHRINK_WRAPPING) {
    parseError(
        "Only small component algorith supports shrink wrapping", lineNum);
    exit(1);
	}
	continuousProblem.algorithm = ALGORITHM_SIMPLE_IMPL;
}
#line 4765 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 119:
#line 2529 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ALG_SIMPLE_COMP");
	if(continuousProblem.precondition == SHRINK_WRAPPING) {
    parseError(
        "Only small component algorith supports shrink wrapping", lineNum);
    exit(1);
	}
	continuousProblem.algorithm = ALGORITHM_SIMPLE_COMP;
}
#line 4779 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 120:
#line 2540 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ALG_SMALL_COMP");
	continuousProblem.algorithm = ALGORITHM_SMALL_COMP;
}
#line 4788 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 121:
#line 2546 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log(sbuilder() << "type: " << typeid(continuousProblem).name());
	logger.log("ALG_DEF");
	continuousProblem.algorithm = ALGORITHM_DEFAULT;
}
#line 4798 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 122:
#line 2551 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log(sbuilder() << "type: " << typeid(continuousProblem).name());
	logger.log("ALG_DEF");
	continuousProblem.algorithm = ALGORITHM_DEFAULT;
}
#line 4808 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 123:
#line 2559 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("DECOMPOSITION");
}
#line 4816 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 124:
#line 2563 "modelParser.y" /* yacc.c:1646  */
    {
  logger.force("none");
  int varNumber = continuousProblem.stateVarNames.size();
  vector<int> vs;
  for(int i = 0; i < varNumber; i++) {
    vs.push_back(i);
  }
  continuousProblem.components.push_back(vs);
}
#line 4830 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 125:
#line 2572 "modelParser.y" /* yacc.c:1646  */
    {
  logger.force("NO DECOMPOSITION");
  if(continuousProblem.algorithm == ALGORITHM_SMALL_COMP) {
    parseError("No decompostion given for composition algorithm.", lineNum);
    exit(1);
  }
}
#line 4842 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 126:
#line 2580 "modelParser.y" /* yacc.c:1646  */
    { }
#line 4848 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 127:
#line 2582 "modelParser.y" /* yacc.c:1646  */
    { }
#line 4854 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 128:
#line 2585 "modelParser.y" /* yacc.c:1646  */
    {
  for(vector<int>::iterator varIt = (yyvsp[-1].iVec)->begin(); varIt < (yyvsp[-1].iVec)->end(); varIt++) {
    int varId = *varIt;
    
    //check if variable is in any of the other components
    for(vector< vector<int> >::iterator it = continuousProblem.components.begin(); 
        it < continuousProblem.components.end(); it++) {
        
      //check if variable is in the component
      if(it->end() != find(it->begin(), it->end(), varId)) {
        char errMsg[MSG_SIZE];
        string vName;
        continuousProblem.getStateVarName(vName, varId);
        sprintf(errMsg, "Variable %s is already in the component.", vName.c_str());
        parseError(errMsg, lineNum);
        exit(1);
      }
    }
  }

  //add this component to all components
  continuousProblem.components.push_back(*((yyvsp[-1].iVec)));
  delete (yyvsp[-1].iVec);
}
#line 4883 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 129:
#line 2611 "modelParser.y" /* yacc.c:1646  */
    {
	//logger.log("compVarIds1");
	//logger.log(sbuilder() << "*$3: " << *$3);
  
	(yyval.iVec) = (yyvsp[-2].iVec);
  int varId = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));
  if(varId < 0) {
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
  
  //check if variable is in the current component
  if((yyval.iVec)->end() != find((yyval.iVec)->begin(), (yyval.iVec)->end(), varId)) {
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Variable %s is already in the component.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
  (yyval.iVec)->push_back(varId);
	delete (yyvsp[0].identifier);
}
#line 4911 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 130:
#line 2636 "modelParser.y" /* yacc.c:1646  */
    {
	//logger.log("compVarIds2");
	//logger.log(sbuilder() << "*$1: " << *$1);
  int varId = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));
  
  if(varId < 0) {
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}
  
	(yyval.iVec) = new vector<int>();
  (yyval.iVec)->push_back(varId);
	delete (yyvsp[0].identifier);
}
#line 4932 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 131:
#line 2656 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log("plotting1");
	int x = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));
	int y = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_INTERVAL;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_INTERVAL;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 4971 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 132:
#line 2692 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log("plotting2");
	int x = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));
	int y = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State Variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_OCTAGON;
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_OCTAGON;
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 5010 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 133:
#line 2728 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log("plotting3");
	int x = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));
	int y = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_GRID;
	continuousProblem.numSections = (int)(yyvsp[-3].dblVal);
	continuousProblem.plotFormat = PLOT_GNUPLOT;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_GRID;
	hybridProblem.numSections = (int)(yyvsp[-3].dblVal);
	hybridProblem.plotFormat = PLOT_GNUPLOT;

	delete (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 5051 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 134:
#line 2766 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log("plotting4");
	int x = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));
	int y = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_INTERVAL;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_INTERVAL;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 5090 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 135:
#line 2802 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log("plotting5");
	int x = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));
	int y = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State Variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_OCTAGON;
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_OCTAGON;
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 5129 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 136:
#line 2838 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log("plotting6");
	int x = continuousProblem.getIDForStateVar(*(yyvsp[-2].identifier));
	int y = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-2].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	continuousProblem.outputAxes.push_back(x);
	continuousProblem.outputAxes.push_back(y);
	continuousProblem.plotSetting = PLOT_GRID;
	continuousProblem.numSections = (int)(yyvsp[-3].dblVal);
	continuousProblem.plotFormat = PLOT_MATLAB;

	hybridProblem.outputAxes.push_back(x);
	hybridProblem.outputAxes.push_back(y);
	hybridProblem.plotSetting = PLOT_GRID;
	hybridProblem.numSections = (int)(yyvsp[-3].dblVal);
	hybridProblem.plotFormat = PLOT_MATLAB;

	delete (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 5170 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 137:
#line 2877 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("INIT1");
	(yyval.pFlowpipe) = new Flowpipe(*(yyvsp[-1].tmVec), *(yyvsp[0].intVec));

	delete (yyvsp[-1].tmVec);
	delete (yyvsp[0].intVec);
}
#line 5182 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 138:
#line 2886 "modelParser.y" /* yacc.c:1646  */
    {
  logger.log("INIT2");
	Interval intZero;
  logger.logVI(*(yyvsp[0].intVec));
	logger.log("INIT2e");
	(yyval.pFlowpipe) = new Flowpipe(*(yyvsp[0].intVec), intZero);

	delete (yyvsp[0].intVec);
}
#line 5196 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 139:
#line 2897 "modelParser.y" /* yacc.c:1646  */
    {
}
#line 5203 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 140:
#line 2902 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	
	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a parameter.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(!continuousProblem.declareTMVar(*(yyvsp[0].identifier)))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareTMVar(*(yyvsp[0].identifier));
	delete (yyvsp[0].identifier);
}
#line 5240 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 141:
#line 2936 "modelParser.y" /* yacc.c:1646  */
    {
	string tVar("local_t");
	continuousProblem.declareTMVar(tVar);
	hybridProblem.declareTMVar(tVar);

	int id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a state variable.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	
	if(id != -1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "%s has already been declared as a parameter.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if(!continuousProblem.declareTMVar(*(yyvsp[0].identifier)))
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	hybridProblem.declareTMVar(*(yyvsp[0].identifier));
	delete (yyvsp[0].identifier);
}
#line 5281 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 142:
#line 2975 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("taylor_model");
	int id = continuousProblem.getIDForStateVar(*(yyvsp[-8].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-8].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	TaylorModel tmTemp(*(yyvsp[-6].poly), I);
	(yyval.tmVec) = (yyvsp[-9].tmVec);
	(yyval.tmVec)->tms[id] = tmTemp;

	delete (yyvsp[-8].identifier);
	delete (yyvsp[-6].poly);
}
#line 5312 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 143:
#line 3002 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmEmpty;
	(yyval.tmVec) = new TaylorModelVec;

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		(yyval.tmVec)->tms.push_back(tmEmpty);
	}
}
#line 5326 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 144:
#line 3014 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForTMVar(*(yyvsp[-6].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*(yyvsp[-6].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	(yyval.intVec) = (yyvsp[-7].intVec);
	(*(yyval.intVec))[id] = I;

	delete (yyvsp[-6].identifier);
}
#line 5354 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 145:
#line 3039 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.intVec) = new vector<Interval>( continuousProblem.tmVarNames.size() );

	Interval intZero;
	(*(yyval.intVec))[0] = intZero;

	int id = continuousProblem.getIDForTMVar(*(yyvsp[-6].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*(yyvsp[-6].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	(*(yyval.intVec))[id] = I;

	delete (yyvsp[-6].identifier);
}
#line 5386 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 146:
#line 3069 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("intervals1");
	int id = continuousProblem.getIDForStateVar(*(yyvsp[-6].identifier));
	logger.log(*(yyvsp[-6].identifier));
	
	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-6].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	(yyval.intVec) = (yyvsp[-7].intVec);
	(*(yyval.intVec))[id] = I;
	logger.logVI(*(yyval.intVec));

	delete (yyvsp[-6].identifier);
}
#line 5417 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 147:
#line 3096 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("intervals2");
	//exit(1);
	int numVars = continuousProblem.stateVarNames.size();
	logger.log(sbuilder() << "stateVarNames.size = " << continuousProblem.stateVarNames.size());
	(yyval.intVec) = new vector<Interval>(numVars);

	string tVar("local_t");
	continuousProblem.declareTMVar(tVar);

	char name[NAME_SIZE];

	logger.log(local_var_name);
	for(int i=0; i<numVars; ++i)
	{
		sprintf(name, "%s%d", local_var_name, i+1);
		string tmVarName(name);
		logger.log(sbuilder() << "name: " << name);
		continuousProblem.declareTMVar(tmVarName);
	}
}
#line 5443 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 148:
#line 3120 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ODE1");
	(yyval.tmVec) = (yyvsp[-4].tmVec);
	int id = continuousProblem.getIDForStateVar(*(yyvsp[-3].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*(yyvsp[0].poly), intZero);
	(yyval.tmVec)->tms[id] = tmTemp;

	delete (yyvsp[-3].identifier);
	delete (yyvsp[0].poly);
}
#line 5468 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 149:
#line 3141 "modelParser.y" /* yacc.c:1646  */
    {	
	logger.log("ODE2");
	int numVars = continuousProblem.stateVarNames.size();
	logger.log(sbuilder() << "numVars: " << numVars);
	(yyval.tmVec) = new TaylorModelVec;
	TaylorModel tmTemp;
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		(yyval.tmVec)->tms.push_back(tmTemp);
	}
}
#line 5486 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 150:
#line 3157 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.strVec) = (yyvsp[-4].strVec);

	int id = continuousProblem.getIDForStateVar(*(yyvsp[-3].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*(yyval.strVec))[id] = (*(yyvsp[0].identifier));

	delete (yyvsp[-3].identifier);
	delete (yyvsp[0].identifier);
}
#line 5509 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 151:
#line 3176 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size();
	(yyval.strVec) = new vector<string>;

	string empty;
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		(yyval.strVec)->push_back(empty);
	}
}
#line 5526 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 152:
#line 3194 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.tmVec) = (yyvsp[-4].tmVec);

	int id = continuousProblem.getIDForStateVar(*(yyvsp[-3].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-3].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	Interval intZero;
	TaylorModel tmTemp(*(yyvsp[0].poly), intZero);
	(yyval.tmVec)->tms[id] = tmTemp;

	delete (yyvsp[-3].identifier);
	delete (yyvsp[0].poly);
}
#line 5551 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 153:
#line 3215 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size();

	(yyval.tmVec) = new TaylorModelVec;
	TaylorModel tmTemp;
	Interval intZero;

	for(int i=0; i<numVars; ++i)
	{
		(yyval.tmVec)->tms.push_back(tmTemp);
	}
}
#line 5568 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 154:
#line 3236 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) += (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5579 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 155:
#line 3244 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) -= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5590 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 156:
#line 3252 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-1].poly); 
}
#line 5598 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 157:
#line 3257 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) *= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5609 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 158:
#line 3265 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int) (yyvsp[0].dblVal);

	if(exp == 0)
	{
		Interval I(1);
		(yyval.poly) = new Polynomial(I, continuousProblem.tmVarNames.size());
	}
	else
	{
		(yyval.poly) = new Polynomial;
		(*(yyvsp[-2].poly)).pow(*(yyval.poly), exp);
	}

	delete (yyvsp[-2].poly);
}
#line 5630 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 159:
#line 3283 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I(-1);
	(yyval.poly) = (yyvsp[0].poly);
	(yyval.poly)->mul_assign(I);
}
#line 5640 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 160:
#line 3290 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForTMVar(*(yyvsp[0].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int numVars = continuousProblem.tmVarNames.size();
	Interval I(1);

	vector<int> degrees;
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	degrees[id] = 1;
	Monomial monomial(I, degrees);

	(yyval.poly) = new Polynomial(monomial);
	delete (yyvsp[0].identifier);
}
#line 5671 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 161:
#line 3318 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.tmVarNames.size();
	Interval I((yyvsp[0].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 5681 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 162:
#line 3336 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ODEpolynomial+");
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) += (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5693 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 163:
#line 3345 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) -= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5704 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 164:
#line 3353 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-1].poly); 
}
#line 5712 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 165:
#line 3358 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ODEpolynomial*");
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) *= (*(yyvsp[0].poly));
	delete (yyvsp[0].poly);
}
#line 5723 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 166:
#line 3366 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int) (yyvsp[0].dblVal);

	if(exp == 0)
	{
		Interval I(1);
		(yyval.poly) = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		(yyval.poly) = new Polynomial;
		(*(yyvsp[-2].poly)).pow(*(yyval.poly), exp);
	}

	delete (yyvsp[-2].poly);
}
#line 5744 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 167:
#line 3384 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I(-1);
	(yyval.poly) = (yyvsp[0].poly);
	(yyval.poly)->mul_assign(I);
}
#line 5754 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 168:
#line 3391 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ODEpolynomialIDENT");
	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));
	logger.log(sbuilder() << "idPar: " << id);
	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));

		int numVars = continuousProblem.stateVarNames.size()+1;
		(yyval.poly) = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));
		logger.log(sbuilder() << "idSV: " << id);

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		(yyval.poly) = new Polynomial(monomial);
	}

	delete (yyvsp[0].identifier);
}
#line 5801 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 169:
#line 3435 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ODEpolynomial[,]");
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I((yyvsp[-3].dblVal), (yyvsp[-1].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 5812 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 170:
#line 3443 "modelParser.y" /* yacc.c:1646  */
    {
	logger.log("ODEpolynomialNUM");
	logger.log(sbuilder() << "$1: " << (yyvsp[0].dblVal));
	int numVars = continuousProblem.stateVarNames.size()+1;
	logger.log(sbuilder() << "numVars+1:" << numVars);
	Interval I((yyvsp[0].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
	logger.logPoly((yyval.poly));
}
#line 5826 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 171:
#line 3457 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) += (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5837 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 172:
#line 3465 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) -= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5848 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 173:
#line 3473 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-1].poly); 
}
#line 5856 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 174:
#line 3478 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) *= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5867 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 175:
#line 3486 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int) (yyvsp[0].dblVal);

	if(exp == 0)
	{
		Interval I(1);
		(yyval.poly) = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		(yyval.poly) = new Polynomial;
		(*(yyvsp[-2].poly)).pow(*(yyval.poly), exp);
	}

	delete (yyvsp[-2].poly);
}
#line 5888 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 176:
#line 3504 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I(-1);
	(yyval.poly) = (yyvsp[0].poly);
	(yyval.poly)->mul_assign(I);
}
#line 5898 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 177:
#line 3511 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));

		int numVars = continuousProblem.stateVarNames.size()+1;
		(yyval.poly) = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		(yyval.poly) = new Polynomial(monomial);
	}

	delete (yyvsp[0].identifier);
}
#line 5943 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 178:
#line 3553 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I((yyvsp[0].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 5953 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 179:
#line 3564 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) += (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5964 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 180:
#line 3572 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) -= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5975 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 181:
#line 3580 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-1].poly); 
}
#line 5983 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 182:
#line 3585 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) *= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 5994 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 183:
#line 3593 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int) (yyvsp[0].dblVal);

	if(exp == 0)
	{
		Interval I(1);
		(yyval.poly) = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		(yyval.poly) = new Polynomial;
		(*(yyvsp[-2].poly)).pow(*(yyval.poly), exp);
	}

	delete (yyvsp[-2].poly);
}
#line 6015 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 184:
#line 3611 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I(-1);
	(yyval.poly) = (yyvsp[0].poly);
	(yyval.poly)->mul_assign(I);
}
#line 6025 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 185:
#line 3618 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));

		int numVars = continuousProblem.stateVarNames.size()+1;
		(yyval.poly) = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		(yyval.poly) = new Polynomial(monomial);
	}

	delete (yyvsp[0].identifier);
}
#line 6070 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 186:
#line 3660 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I((yyvsp[0].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 6080 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 187:
#line 3672 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForStateVar(*(yyvsp[-8].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-8].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	TaylorModel tmTemp(*(yyvsp[-6].poly), I);
	(yyval.tmVec) = (yyvsp[-9].tmVec);
	(yyval.tmVec)->tms[id] = tmTemp;

	delete (yyvsp[-8].identifier);
	delete (yyvsp[-6].poly);
}
#line 6110 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 188:
#line 3699 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmEmpty;
	(yyval.tmVec) = new TaylorModelVec;

	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i)
	{
		(yyval.tmVec)->tms.push_back(tmEmpty);
	}

	int id = continuousProblem.getIDForStateVar(*(yyvsp[-8].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[-8].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if((yyvsp[-3].dblVal) > (yyvsp[-1].dblVal))
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	Interval I((yyvsp[-3].dblVal),(yyvsp[-1].dblVal));
	TaylorModel tmTemp(*(yyvsp[-6].poly), I);

	(yyval.tmVec)->tms[id] = tmTemp;

	delete (yyvsp[-8].identifier);
	delete (yyvsp[-6].poly);
}
#line 6148 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 189:
#line 3744 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) += (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 6159 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 190:
#line 3752 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) -= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 6170 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 191:
#line 3760 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-1].poly); 
}
#line 6178 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 192:
#line 3765 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) *= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 6189 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 193:
#line 3773 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int) (yyvsp[0].dblVal);

	if(exp == 0)
	{
		Interval I(1);
		(yyval.poly) = new Polynomial(I, continuousProblem.tmVarNames.size());
	}
	else
	{
		(yyval.poly) = new Polynomial;
		(*(yyvsp[-2].poly)).pow(*(yyval.poly), exp);
	}

	delete (yyvsp[-2].poly);
}
#line 6210 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 194:
#line 3791 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I(-1);
	(yyval.poly) = (yyvsp[0].poly);
	(yyval.poly)->mul_assign(I);
}
#line 6220 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 195:
#line 3798 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForTMVar(*(yyvsp[0].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	int numVars = continuousProblem.tmVarNames.size();
	Interval I(1);

	vector<int> degrees;
	for(int i=0; i<numVars; ++i)
	{
		degrees.push_back(0);
	}

	degrees[id] = 1;
	Monomial monomial(I, degrees);

	(yyval.poly) = new Polynomial(monomial);
	delete (yyvsp[0].identifier);
}
#line 6251 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 196:
#line 3826 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.tmVarNames.size();
	Interval I((yyvsp[-3].dblVal), (yyvsp[-1].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 6261 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 197:
#line 3833 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.tmVarNames.size();
	Interval I((yyvsp[0].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 6271 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 198:
#line 3860 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.ptm) = (yyvsp[-2].ptm);
	(yyvsp[-2].ptm)->add_assign(*(yyvsp[0].ptm));
	delete (yyvsp[0].ptm);
}
#line 6281 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 199:
#line 3867 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.ptm) = (yyvsp[-2].ptm);
	(yyvsp[-2].ptm)->sub_assign(*(yyvsp[0].ptm));
	delete (yyvsp[0].ptm);
}
#line 6291 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 200:
#line 3874 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.ptm) = (yyvsp[-2].ptm);

	Interval intPoly1, intPoly2, intTrunc;

	(yyvsp[0].ptm)->polyRangeNormal(intPoly2, parseSetting.step_exp_table);
	(yyvsp[-2].ptm)->mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, *(yyvsp[0].ptm), intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);

	parseSetting.ranges.push_back(intPoly1);
	parseSetting.ranges.push_back(intPoly2);
	parseSetting.ranges.push_back(intTrunc);

	delete (yyvsp[0].ptm);
}
#line 6310 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 201:
#line 3890 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.ptm) = (yyvsp[-1].ptm);
}
#line 6318 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 202:
#line 3895 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmTemp;
	(yyvsp[0].ptm)->rec_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(yyval.ptm) = (yyvsp[-2].ptm);

	Interval intPoly1, intPoly2, intTrunc;

	tmTemp.polyRangeNormal(intPoly2, parseSetting.step_exp_table);
	(yyvsp[-2].ptm)->mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, tmTemp, intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);

	parseSetting.ranges.push_back(intPoly1);
	parseSetting.ranges.push_back(intPoly2);
	parseSetting.ranges.push_back(intTrunc);

	delete (yyvsp[0].ptm);
}
#line 6340 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 203:
#line 3914 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmTemp;
	(yyvsp[-1].ptm)->exp_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*(yyvsp[-1].ptm) = tmTemp;
	(yyval.ptm) = (yyvsp[-1].ptm);
}
#line 6352 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 204:
#line 3923 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmTemp;
	(yyvsp[-1].ptm)->sin_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*(yyvsp[-1].ptm) = tmTemp;
	(yyval.ptm) = (yyvsp[-1].ptm);
}
#line 6364 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 205:
#line 3932 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmTemp;
	(yyvsp[-1].ptm)->cos_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*(yyvsp[-1].ptm) = tmTemp;
	(yyval.ptm) = (yyvsp[-1].ptm);
}
#line 6376 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 206:
#line 3941 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmTemp;
	(yyvsp[-1].ptm)->log_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*(yyvsp[-1].ptm) = tmTemp;
	(yyval.ptm) = (yyvsp[-1].ptm);
}
#line 6388 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 207:
#line 3950 "modelParser.y" /* yacc.c:1646  */
    {
	TaylorModel tmTemp;
	(yyvsp[-1].ptm)->sqrt_taylor(tmTemp, parseSetting.ranges, parseSetting.step_exp_table, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	*(yyvsp[-1].ptm) = tmTemp;
	(yyval.ptm) = (yyvsp[-1].ptm);
}
#line 6400 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 208:
#line 3959 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int)(yyvsp[0].dblVal);

	if(exp == 0)
	{
		Interval I(1);
		(yyval.ptm) = new TaylorModel(I, continuousProblem.stateVarNames.size()+1);
	}
	else
	{
		TaylorModel res = *(yyvsp[-2].ptm);
		TaylorModel pow = *(yyvsp[-2].ptm);
		
		Interval intPoly1, intPoly2, intTrunc;
		
		for(int degree = exp - 1; degree > 0;)
		{
			pow.polyRangeNormal(intPoly2, parseSetting.step_exp_table);
		
			if(degree & 1)
			{
				res.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, pow, intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);
				
				parseSetting.ranges.push_back(intPoly1);
				parseSetting.ranges.push_back(intPoly2);
				parseSetting.ranges.push_back(intTrunc);
			}
			
			degree >>= 1;
			
			if(degree > 0)
			{
				pow.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, pow, intPoly2, parseSetting.step_exp_table, parseSetting.order, parseSetting.cutoff_threshold);
				
				parseSetting.ranges.push_back(intPoly1);
				parseSetting.ranges.push_back(intPoly2);
				parseSetting.ranges.push_back(intTrunc);
			}
		}
		
		(yyval.ptm) = new TaylorModel(res);
	}

	delete (yyvsp[-2].ptm);
}
#line 6450 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 209:
#line 4006 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I(-1);
	(yyval.ptm) = (yyvsp[0].ptm);
	(yyval.ptm)->mul_assign(I);
}
#line 6460 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 210:
#line 4013 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.ptm) = new TaylorModel;
	(*(yyval.ptm)) = parseSetting.flowpipe.tms[id];

	Interval intTemp;
	(*(yyval.ptm)).expansion.ctrunc_normal(intTemp, parseSetting.step_exp_table, parseSetting.order);
	parseSetting.ranges.push_back(intTemp);
	(*(yyval.ptm)).remainder += intTemp;

	delete (yyvsp[0].identifier);
}
#line 6486 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 211:
#line 4036 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I((yyvsp[-3].dblVal), (yyvsp[-1].dblVal));
	(yyval.ptm) = new TaylorModel(I, numVars);
}
#line 6496 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 212:
#line 4059 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.pint) = (yyvsp[-2].pint);
	(*(yyval.pint)) += (*(yyvsp[0].pint));
	delete (yyvsp[0].pint);
}
#line 6506 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 213:
#line 4066 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.pint) = (yyvsp[-2].pint);
	(*(yyval.pint)) -= (*(yyvsp[0].pint));
	delete (yyvsp[0].pint);
}
#line 6516 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 214:
#line 4073 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.pint) = new Interval;

	(*(yyval.pint)) = (*parseSetting.iterRange) * (*(yyvsp[0].pint));
	++parseSetting.iterRange;
	(*(yyval.pint)) += (*parseSetting.iterRange) * (*(yyvsp[-2].pint));
	(*(yyval.pint)) += (*(yyvsp[-2].pint)) * (*(yyvsp[0].pint));
	++parseSetting.iterRange;
	(*(yyval.pint)) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete (yyvsp[-2].pint);
	delete (yyvsp[0].pint);
}
#line 6535 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 215:
#line 4089 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.pint) = (yyvsp[-1].pint);
}
#line 6543 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 216:
#line 4094 "modelParser.y" /* yacc.c:1646  */
    {
	Interval intTemp;
	rec_taylor_only_remainder(intTemp, *(yyvsp[0].pint), parseSetting.iterRange, parseSetting.order);

	(yyval.pint) = new Interval;

	(*(yyval.pint)) = (*parseSetting.iterRange) * intTemp;
	++parseSetting.iterRange;
	(*(yyval.pint)) += (*parseSetting.iterRange) * (*(yyvsp[-2].pint));
	(*(yyval.pint)) += (*(yyvsp[-2].pint)) * intTemp;
	++parseSetting.iterRange;
	(*(yyval.pint)) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete (yyvsp[-2].pint);
	delete (yyvsp[0].pint);
}
#line 6565 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 217:
#line 4113 "modelParser.y" /* yacc.c:1646  */
    {
	Interval intTemp;
	exp_taylor_only_remainder(intTemp, *(yyvsp[-1].pint), parseSetting.iterRange, parseSetting.order);

	(*(yyvsp[-1].pint)) = intTemp;
	(yyval.pint) = (yyvsp[-1].pint);
}
#line 6577 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 218:
#line 4122 "modelParser.y" /* yacc.c:1646  */
    {
	Interval intTemp;
	sin_taylor_only_remainder(intTemp, *(yyvsp[-1].pint), parseSetting.iterRange, parseSetting.order);

	(*(yyvsp[-1].pint)) = intTemp;
	(yyval.pint) = (yyvsp[-1].pint);
}
#line 6589 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 219:
#line 4131 "modelParser.y" /* yacc.c:1646  */
    {
	Interval intTemp;
	cos_taylor_only_remainder(intTemp, *(yyvsp[-1].pint), parseSetting.iterRange, parseSetting.order);

	(*(yyvsp[-1].pint)) = intTemp;
	(yyval.pint) = (yyvsp[-1].pint);
}
#line 6601 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 220:
#line 4140 "modelParser.y" /* yacc.c:1646  */
    {
	Interval intTemp;
	log_taylor_only_remainder(intTemp, *(yyvsp[-1].pint), parseSetting.iterRange, parseSetting.order);

	(*(yyvsp[-1].pint)) = intTemp;
	(yyval.pint) = (yyvsp[-1].pint);
}
#line 6613 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 221:
#line 4149 "modelParser.y" /* yacc.c:1646  */
    {
	Interval intTemp;
	sqrt_taylor_only_remainder(intTemp, *(yyvsp[-1].pint), parseSetting.iterRange, parseSetting.order);

	(*(yyvsp[-1].pint)) = intTemp;
	(yyval.pint) = (yyvsp[-1].pint);
}
#line 6625 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 222:
#line 4158 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int)(yyvsp[0].dblVal);

	if(exp == 0)
	{
		Interval intZero;
		(*(yyvsp[-2].pint)) = intZero;
		(yyval.pint) = (yyvsp[-2].pint);
	}
	else
	{
		Interval res(*(yyvsp[-2].pint));
		Interval pow(*(yyvsp[-2].pint));
		
		Interval intPoly1, intPoly2, intTrunc;
		
		for(int degree = exp - 1; degree > 0;)
		{
			if(degree & 1)
			{
				Interval intTemp;
				intTemp = (*parseSetting.iterRange) * pow;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange) * res;
				intTemp += pow * res;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange);
				++parseSetting.iterRange;
			
				res = intTemp;
			}
			
			degree >>= 1;
			
			if(degree > 0)
			{
				Interval intTemp;
				intTemp = (*parseSetting.iterRange) * pow;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange) * pow;
				intTemp += pow * pow;
				++parseSetting.iterRange;
				intTemp += (*parseSetting.iterRange);
				++parseSetting.iterRange;
			
				pow = intTemp;
			}
		}
		
		(yyval.pint) = new Interval(res);
	}

	delete (yyvsp[-2].pint);
}
#line 6684 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 223:
#line 4214 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.pint) = (yyvsp[0].pint);
	(yyval.pint)->inv_assign();
}
#line 6693 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 224:
#line 4220 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.pint) = new Interval;
	(*(yyval.pint)) = parseSetting.flowpipe.tms[id].getRemainder();
	
	(*(yyval.pint)) += (*parseSetting.iterRange);
	++parseSetting.iterRange;

	delete (yyvsp[0].identifier);
}
#line 6717 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 225:
#line 4241 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.pint) = new Interval;
}
#line 6725 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 226:
#line 4263 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) += (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 6736 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 227:
#line 4271 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) -= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 6747 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 228:
#line 4279 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) *= (*(yyvsp[0].poly));
	(yyval.poly)->nctrunc(parseSetting.order);
	(yyval.poly)->cutoff(parseSetting.cutoff_threshold);

	delete (yyvsp[0].poly);
}
#line 6760 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 229:
#line 4289 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-1].poly);
}
#line 6768 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 230:
#line 4294 "modelParser.y" /* yacc.c:1646  */
    {
	Polynomial polyTemp;
	(yyvsp[0].poly)->rec_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*(yyvsp[-2].poly)) *= polyTemp;
	(yyvsp[-2].poly)->nctrunc(parseSetting.order);
	(yyval.poly) = (yyvsp[-2].poly);
	(yyval.poly)->cutoff(parseSetting.cutoff_threshold);

	delete (yyvsp[0].poly);
}
#line 6784 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 231:
#line 4307 "modelParser.y" /* yacc.c:1646  */
    {
	Polynomial polyTemp;
	(yyvsp[-1].poly)->exp_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*(yyvsp[-1].poly)) = polyTemp;
	(yyval.poly) = (yyvsp[-1].poly);
}
#line 6796 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 232:
#line 4316 "modelParser.y" /* yacc.c:1646  */
    {
	Polynomial polyTemp;
	(yyvsp[-1].poly)->sin_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*(yyvsp[-1].poly)) = polyTemp;
	(yyval.poly) = (yyvsp[-1].poly);
}
#line 6808 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 233:
#line 4325 "modelParser.y" /* yacc.c:1646  */
    {
	Polynomial polyTemp;
	(yyvsp[-1].poly)->cos_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*(yyvsp[-1].poly)) = polyTemp;
	(yyval.poly) = (yyvsp[-1].poly);
}
#line 6820 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 234:
#line 4334 "modelParser.y" /* yacc.c:1646  */
    {
	Polynomial polyTemp;
	(yyvsp[-1].poly)->log_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*(yyvsp[-1].poly)) = polyTemp;
	(yyval.poly) = (yyvsp[-1].poly);
}
#line 6832 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 235:
#line 4343 "modelParser.y" /* yacc.c:1646  */
    {
	Polynomial polyTemp;
	(yyvsp[-1].poly)->sqrt_taylor(polyTemp, continuousProblem.stateVarNames.size()+1, parseSetting.order, parseSetting.cutoff_threshold);

	(*(yyvsp[-1].poly)) = polyTemp;
	(yyval.poly) = (yyvsp[-1].poly);
}
#line 6844 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 236:
#line 4352 "modelParser.y" /* yacc.c:1646  */
    {
	int exp = (int) (yyvsp[0].dblVal);
	
	(yyval.poly) = new Polynomial;
	
	(*(yyvsp[-2].poly)).pow(*(yyval.poly), exp, parseSetting.order);
	(yyval.poly)->cutoff(parseSetting.cutoff_threshold);
}
#line 6857 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 237:
#line 4362 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[0].poly);
	(yyval.poly)->inv_assign();
}
#line 6866 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 238:
#line 4368 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.poly) = new Polynomial;
	parseSetting.flowpipe.tms[id].getExpansion(*(yyval.poly));
	
	(*(yyval.poly)).nctrunc(parseSetting.order);
}
#line 6887 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 239:
#line 4386 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I((yyvsp[-3].dblVal), (yyvsp[-1].dblVal));
	(yyval.poly) = new Polynomial(I, continuousProblem.stateVarNames.size()+1);
}
#line 6896 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 240:
#line 4406 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '+';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 6910 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 241:
#line 4417 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '-';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 6924 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 242:
#line 4428 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '*';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 6938 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 243:
#line 4439 "modelParser.y" /* yacc.c:1646  */
    {
	string str;
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 6952 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 244:
#line 4450 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '/';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 6966 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 245:
#line 4461 "modelParser.y" /* yacc.c:1646  */
    {
	string str("exp");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 6980 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 246:
#line 4472 "modelParser.y" /* yacc.c:1646  */
    {
	string str("sin");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 6994 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 247:
#line 4483 "modelParser.y" /* yacc.c:1646  */
    {
	string str("cos");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7008 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 248:
#line 4494 "modelParser.y" /* yacc.c:1646  */
    {
	string str("log");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7022 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 249:
#line 4505 "modelParser.y" /* yacc.c:1646  */
    {
	string str("sqrt");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7036 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 250:
#line 4516 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += '^';

	char strNum[NUM_LENGTH];
	sprintf(strNum, "%d", (int)(yyvsp[0].dblVal));
	string num(strNum);
	(*(yyvsp[-2].identifier)) += num;

	(yyval.identifier) = (yyvsp[-2].identifier);
}
#line 7051 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 251:
#line 4528 "modelParser.y" /* yacc.c:1646  */
    {
	string str;
	str += '-';
	str += (*(yyvsp[0].identifier));
	(*(yyvsp[0].identifier)) = str;

	(yyval.identifier) = (yyvsp[0].identifier);
}
#line 7064 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 252:
#line 4538 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));

		(yyval.identifier) = new string;
		range.toString(*(yyval.identifier));
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		(yyval.identifier) = (yyvsp[0].identifier);
	}
}
#line 7095 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 253:
#line 4566 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.identifier) = new string;
	char strNum_lo[NUM_LENGTH], strNum_up[NUM_LENGTH];
	sprintf(strNum_lo, "%.20e", (yyvsp[-3].dblVal));
	sprintf(strNum_up, "%.20e", (yyvsp[-1].dblVal));

	string num_lo(strNum_lo);
	string num_up(strNum_up);

	(*(yyval.identifier)) += '[';
	(*(yyval.identifier)) += num_lo;
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += ',';
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += num_up;
	(*(yyval.identifier)) += ']';
}
#line 7117 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 254:
#line 4585 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.identifier) = new string;
	char strNum[NUM_LENGTH];
	sprintf(strNum, "%.20e", (yyvsp[0].dblVal));
	string num(strNum);

	(*(yyval.identifier)) += '[';
	(*(yyval.identifier)) += num;
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += ',';
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += num;
	(*(yyval.identifier)) += ']';
}
#line 7136 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 255:
#line 4619 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '+';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 7150 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 256:
#line 4630 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '-';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 7164 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 257:
#line 4641 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '*';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 7178 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 258:
#line 4652 "modelParser.y" /* yacc.c:1646  */
    {
	string str;
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7192 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 259:
#line 4663 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += '/';
	(*(yyvsp[-2].identifier)) += ' ';
	(*(yyvsp[-2].identifier)) += (*(yyvsp[0].identifier));

	(yyval.identifier) = (yyvsp[-2].identifier);
	delete (yyvsp[0].identifier);
}
#line 7206 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 260:
#line 4674 "modelParser.y" /* yacc.c:1646  */
    {
	string str("exp");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7220 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 261:
#line 4685 "modelParser.y" /* yacc.c:1646  */
    {
	string str("sin");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7234 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 262:
#line 4696 "modelParser.y" /* yacc.c:1646  */
    {
	string str("cos");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7248 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 263:
#line 4707 "modelParser.y" /* yacc.c:1646  */
    {
	string str("log");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7262 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 264:
#line 4718 "modelParser.y" /* yacc.c:1646  */
    {
	string str("sqrt");
	str += '(';
	str += (*(yyvsp[-1].identifier));
	str += ')';
	(*(yyvsp[-1].identifier)) = str;

	(yyval.identifier) = (yyvsp[-1].identifier);
}
#line 7276 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 265:
#line 4729 "modelParser.y" /* yacc.c:1646  */
    {
	(*(yyvsp[-2].identifier)) += '^';

	char strNum[NUM_LENGTH];
	sprintf(strNum, "%d", (int)(yyvsp[0].dblVal));
	string num(strNum);
	(*(yyvsp[-2].identifier)) += num;

	(yyval.identifier) = (yyvsp[-2].identifier);
}
#line 7291 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 266:
#line 4741 "modelParser.y" /* yacc.c:1646  */
    {
	string str;
	str += '-';
	str += (*(yyvsp[0].identifier));
	(*(yyvsp[0].identifier)) = str;

	(yyval.identifier) = (yyvsp[0].identifier);
}
#line 7304 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 267:
#line 4751 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(yyval.identifier) = (yyvsp[0].identifier);
}
#line 7322 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 268:
#line 4766 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.identifier) = new string;
	char strNum[NUM_LENGTH];
	sprintf(strNum, "%.20e", ((yyvsp[-3].dblVal)+(yyvsp[-1].dblVal))/2);

	string num(strNum);

	(*(yyval.identifier)) += '[';
	(*(yyval.identifier)) += num;
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += ',';
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += num;
	(*(yyval.identifier)) += ']';
}
#line 7342 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 269:
#line 4783 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.identifier) = new string;
	char strNum[NUM_LENGTH];
	sprintf(strNum, "%.20e", (yyvsp[0].dblVal));
	string num(strNum);

	(*(yyval.identifier)) += '[';
	(*(yyval.identifier)) += num;
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += ',';
	(*(yyval.identifier)) += ' ';
	(*(yyval.identifier)) += num;
	(*(yyval.identifier)) += ']';
}
#line 7361 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 270:
#line 4814 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) += (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 7372 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 271:
#line 4822 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-2].poly);
	(*(yyval.poly)) -= (*(yyvsp[0].poly));

	delete (yyvsp[0].poly);
}
#line 7383 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 272:
#line 4830 "modelParser.y" /* yacc.c:1646  */
    {
	if((yyvsp[-2].poly)->degree() + (yyvsp[0].poly)->degree() > 1)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Polynomial is not linear.");
		parseError(errMsg, lineNum);
		exit(1);
	}
	else
	{
		(yyval.poly) = (yyvsp[-2].poly);
		(*(yyval.poly)) *= (*(yyvsp[0].poly));
		delete (yyvsp[0].poly);
	}
}
#line 7403 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 273:
#line 4847 "modelParser.y" /* yacc.c:1646  */
    {
	(yyval.poly) = (yyvsp[-1].poly); 
}
#line 7411 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 274:
#line 4852 "modelParser.y" /* yacc.c:1646  */
    {
	Interval I(-1);
	(yyval.poly) = (yyvsp[0].poly);
	(yyval.poly)->mul_assign(I);
}
#line 7421 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 275:
#line 4859 "modelParser.y" /* yacc.c:1646  */
    {
	int id = continuousProblem.getIDForPar(*(yyvsp[0].identifier));

	if(id >= 0)
	{
		Interval range;
		continuousProblem.getRangeForPar(range, *(yyvsp[0].identifier));

		int numVars = continuousProblem.stateVarNames.size()+1;
		(yyval.poly) = new Polynomial(range, numVars);
	}
	else
	{
		id = continuousProblem.getIDForStateVar(*(yyvsp[0].identifier));

		if(id < 0)
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "State variable %s is not declared.", (*(yyvsp[0].identifier)).c_str());
			parseError(errMsg, lineNum);
			exit(1);
		}

		int numVars = continuousProblem.stateVarNames.size()+1;
		Interval I(1);

		vector<int> degrees;
		for(int i=0; i<numVars; ++i)
		{
			degrees.push_back(0);
		}

		degrees[id+1] = 1;
		Monomial monomial(I, degrees);

		(yyval.poly) = new Polynomial(monomial);
	}

	delete (yyvsp[0].identifier);
}
#line 7466 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 276:
#line 4901 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I((yyvsp[0].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 7476 "modelParser.tab.c" /* yacc.c:1646  */
    break;

  case 277:
#line 4908 "modelParser.y" /* yacc.c:1646  */
    {
	int numVars = continuousProblem.stateVarNames.size()+1;
	Interval I((yyvsp[-3].dblVal), (yyvsp[-1].dblVal));
	(yyval.poly) = new Polynomial(I, numVars);
}
#line 7486 "modelParser.tab.c" /* yacc.c:1646  */
    break;


#line 7490 "modelParser.tab.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 4935 "modelParser.y" /* yacc.c:1906  */


int yyerror(const char * what)
{
	fprintf(stderr, "Error line %d: %s\n", lineNum,what);
	err=true;
	return 1;
}

int yyerror(string what)
{
	cerr << "Error line "<<lineNum<<" "<<what<<endl;
	err=true;
	return 1;
}
