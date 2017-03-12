/* A Bison parser, made by GNU Bison 3.0.  */

/* Bison interface for Yacc-like parsers in C

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
#line 22 "modelParser.y" /* yacc.c:1909  */

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

#line 155 "modelParser.tab.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_MODELPARSER_TAB_H_INCLUDED  */
