%{
/*
 * prototye for preparing a Flamola model for Omuses
 *
 * use OmSim 'Flat model' as input stream
 * rf, 6/28/95
 *
 * output C code also for system equations
 * rf, 4/21/96
 *
 * convert {variable}^{number} to pow(...)
 * rf, 5/13/96
 *
 * rf, 11/4/96
 *  -- treatment of conditional expressions
 *  -- functions must not have spaces between fname and "("
 *     as this is now reserved for conditional expressions
 *
 * rf, 1/18/97
 *  -- offset to consider discrete states
 *  -- index kk instead of k for external inputs
 *
 * rf, 1/18/97
 *  -- new syntax of Omuses 1.2
 *
 * rf, 8/24/98
 *  -- bugfix (indentifiers may contain underscores)
 *
 /

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#ifndef NULL
#define NULL 0
#endif

/* declare symbol table */
typedef enum {
  FUNCTION,
  VARIABLE,
  GLOBAL,
  INPUT,
  NUMBER
} SymType;

typedef enum {
  STATE = 1,	/* dynamic state variable */
  CONTROL = 2,	/* control input */
  TEMPORAL = 4,	/* intermediate variable of dynamic part */
  RESOLVED = 8, /* varariable is defined by an equation */
  USED = 16	/* variable is used in other equations */
} SymAttr;  

typedef struct SymEntry {
  char *name;
  char *value;
  int number;
  SymType type;
  SymAttr attr;
  struct SymEntry *next;
  struct SymEntry *prev;
} SymEntry;

static void mname_store(char *name);
static void symexp_store(char *name);
static SymEntry *symlook(char *name, SymType type, SymAttr attr);
static SymEntry *symentry;

static void print_maple_defs();

%}

blanc		[ \t]
digit		[0-9]
letter		[A-Za-z_]
id		{letter}({letter}|{digit})*
variable	{id}(\.{id})*
expansion	"{"{variable}"}"|"{}"
global		{id}"::"{id}
number		({digit}+|({digit}*\.{digit}+))([Ee][+-]?{digit}+)?
input		"*continuousinput*"\([^)]*"time()"\)

/*
 * states:
 *  INITIAL	Lex's initial state
 * Snapshot
 *  DECLARE	default and initial values for parameters variables, resp.
 * Parameter part, implicit discrete part
 *  DISCRETE
 * Continuous manipulated
 *  DYNAMIC	dynamic part
 */

%s DECLARE
%s DISCRETE
%s DYNAMIC
%%


"%" {
  while (input() != '\n'); /* eat up comments */
}

<INITIAL>"model "{id} {
  mname_store(yytext+6);
  input(); /* following semicolon */
}

"declare:"|"% Parameters:"|"% Continuous State Variables:" {
  BEGIN(DECLARE);
}
<DECLARE>{variable}/".default" {
  symentry = symlook(yytext, VARIABLE, 0);
}
<DECLARE>{variable}/".initial" {
  symentry = symlook(yytext, VARIABLE, 0);
}
<DECLARE>{number} {
  if (symentry) {
    symentry->value = (char *)malloc(strlen(yytext) + 1);
    strcpy(symentry->value, yytext);
  }
}

"discrete:" {
  printf("    // bound parameters and implicit discrete part\n");
  BEGIN(DISCRETE);
}
<DISCRETE>{variable}/{blanc}*:= {
  printf("    %s", symlook(yytext, VARIABLE, RESOLVED)->name);
}
<DISCRETE,DYNAMIC>{variable}{blanc}*"^"{blanc}*{number} {
  int i;
  for (i = 0; yytext[i] != '^'; i++);
  yytext[i] = '\0';
  printf("pow(%s, %s)", symlook(yytext, VARIABLE, RESOLVED)->name,
                        yytext + i + 1);
  yytext[i] = '^';
}
<DISCRETE>{variable} {
  printf("%s", symlook(yytext, VARIABLE, USED)->name);
}
<DISCRETE>{global} {
  printf("%s", symlook(yytext, GLOBAL, USED)->name);
}
<DISCRETE>{id}/"(" {
  printf("%s", symlook(yytext, FUNCTION, RESOLVED | USED)->name);
}
<DISCRETE>:= {
  printf("=");
}
<DISCRETE>";" {
  printf(";\n");
}
<DISCRETE>. {
  printf("%s", yytext);
}

"dynamic:" {
  printf("\n  // dynamic model equations\n");
  BEGIN(DYNAMIC);
}
<DYNAMIC>^{variable}/{blanc}*"," {
  symlook(yytext, VARIABLE, TEMPORAL);
}
<DYNAMIC>{variable}/"'" {
  printf("  xp[offs+%d]", symlook(yytext, VARIABLE, STATE)->number);
  input(); /* eat up following quote */
}
<DYNAMIC>{variable}/{blanc}*:= {
  /* backing store symentry for detecting control inputs */
  if (strcmp(yytext + strlen(yytext) - 8, ".control") == 0)
    symentry = symlook(yytext, VARIABLE, CONTROL | TEMPORAL);
  else
    symentry = symlook(yytext, VARIABLE, TEMPORAL);
  printf("  %s", symentry->name);
}
<DYNAMIC>{variable} {
  if (strcmp(yytext, "then") == 0)
    printf("? (adouble)(");
  else if (strcmp(yytext, "else") == 0)
    printf("): (adouble)");
  else if (strcmp(yytext, "if") != 0)
    printf("%s", symlook(yytext, VARIABLE, USED)->name);
}
<DYNAMIC>{global} {
  printf("%s", symlook(yytext, GLOBAL, USED)->name);
}
<DYNAMIC>{input} {
  if (symentry->attr & CONTROL) { /* last seen symbol was a control */
    symentry = symlook(yytext, INPUT, CONTROL | USED);
    printf("u[%d]", symentry->number);
  }
  else
    printf("%s[kk]", symlook(yytext, INPUT, USED)->name);
}
<DYNAMIC>{id}/"(" {
  printf("%s", symlook(yytext, FUNCTION, RESOLVED | USED)->name);
}
<DYNAMIC>{number} {
  /* introduce symbols, as Maple can't handle float numbers */
/*  printf("%s", symlook(yytext, NUMBER, RESOLVED | USED)->name); */
  printf("%s", symlook(yytext, NUMBER, RESOLVED | USED)->value);
}
<DYNAMIC>:= {
  printf("=");
}
<DYNAMIC>";" {
  printf(";\n");
}
<DYNAMIC>. {
  printf("%s", yytext);
}
<DYNAMIC>"output:" {
  printf("\n");
  BEGIN(INITIAL);
}

.|\n
%%

/* bottom of list for symbol table */
static SymEntry *symbot;

/*
 * -- store symbol expansion in "symexp"
 * -- store length of name of a model in "mnamelen"
 * -- process current symbol name in "symbol"
 *
 */
static char *symexp;
static int symexpsize;	/* allocated memory */
static int symexplen;	/* currently stored string */
static char *mname;
static int mnamelen;
static char *symbol;
static int symbolsize;
static int nnumbers;	/* count symbols introduced for numbers */
static int nstates;	/* count symbols introduced for states */
static int ncontrols;	/* count symbols introduced for controls */

int main()
{
  SymEntry *entry;
  SymEntry **states;
  int after_entry;
  int i;

  symbot = NULL;

  mname = (char *)0;
  mnamelen = 0;

  symexpsize = 1;
  symexp = (char *)malloc(symexpsize);
  symexp[0] = 0;

  symbolsize = 10;
  symbol = (char *)malloc(symbolsize);

  nnumbers = 0;
  nstates = 0;

  yylex();

  /* sort states according to flat model */
  states = (SymEntry **)malloc(nstates * sizeof(SymEntry *));
  entry = symbot;
  while (entry != NULL) {
    if (entry->attr & STATE) {
      states[entry->number] = entry;
    }
    entry = entry->next;
  }

  printf("  // state assignments\n");
  for (i = 0; i < nstates; i++) {
    entry = states[i];
    printf("  %s = x[offs+%d];\n", entry->name, entry->number);
  }

  printf("\n    // initial state constraints\n");
  for (i = 0; i < nstates; i++) {
    entry = states[i];
    if (entry->type == VARIABLE && entry->value != NULL)
      printf("    x.min[offs+%d] = x.max[offs+%d] = %s;\t// %s\n",
	     entry->number, entry->number, entry->value, entry->name); 
  }

  printf("\n    // initial states\n");
  for (i = 0; i < nstates; i++) {
    entry = states[i];
    if (entry->type == VARIABLE && entry->value != NULL)
      printf("    x.initial[offs+%d] = %s;\t// %s\n",
	     entry->number, entry->value, entry->name); 
  }

  printf("\n  // default values for parameters\n");
  entry = symbot;
  while (entry != NULL) {
    if (!(entry->attr & (STATE | CONTROL | TEMPORAL | RESOLVED))) {
      if (entry->type == VARIABLE && entry->value != NULL)
      	printf("  %s = %s;\n", entry->name, entry->value);
    }
    entry = entry->next;
  }

  printf("\n  // interface elements for unbound variables\n");
  entry = symbot;
  while (entry != NULL) {
    if (entry->attr == USED) {
      switch (entry->type) {
      case INPUT:
	printf("  _ifList.append(new If_RealVec(\"prg_%s\", &%s));\n",
	       entry->name, entry->name);
	break;
      default:
	printf("  _ifList.append(new If_Real(\"prg_%s\", &%s));\n",
	       entry->name, entry->name);
	break;
      }
    }
    entry = entry->next;
  }

  printf("\n  // model inputs and parameters\n");
  entry = symbot;
  while (entry != NULL) {
    if (!(entry->attr & (STATE | CONTROL | TEMPORAL))
	&& (entry->attr & (USED | RESOLVED))) {
      switch (entry->type) {
      case VARIABLE:
      case GLOBAL:
	printf("  double %s;\n", entry->name);
	break;
      case INPUT:
	printf("  VECP %s;\n", entry->name);
	break;
      default:
	break;
      }
    }
    entry = entry->next;
  }

  printf("\n  // dynamic model variables\n");
  entry = symbot;
  after_entry = 0;
  while (entry != NULL) {
    if (entry->attr & (TEMPORAL | STATE)) {
      if (after_entry)
	printf(", ");
      else
        printf("  adouble ");
      printf("%s", entry->name);
      after_entry = 1;
    }
    entry = entry->next;
  }
  printf(";\n");

  printf("\n  # interface elements for unresolved variables\n");
  entry = symbot;
  while (entry != NULL) {
    if (entry->attr == USED && entry->value == NULL) {
      printf("  prg_%s\t\n", entry->name);
    }
    entry = entry->next;
  }

  return 0;
}


static void print_maple_defs()
{
  SymEntry *entry;
  int after_entry;

  /* define set of float constants */
  printf("floats := {\n");
  entry = symbot;
  after_entry = 0;
  while (entry != NULL) {
    if (entry->type == NUMBER) {
      if (after_entry)
	printf(",\n");
      printf("%s = %s", entry->name, entry->value);
      after_entry = 1;
    }
    entry = entry->next;
  }
  printf("\n}:\n");

  /* define set of temporals */
  printf("temporals := {\n");
  entry = symbot;
  after_entry = 0;
  while (entry != NULL) {
    if (entry->attr & TEMPORAL) {
      if (after_entry)
	printf(",\n");
      printf("%s", entry->name);
      after_entry = 1;
    }
    entry = entry->next;
  }
  printf("\n}:\n");

  /* define set of states */
  printf("states := [\n");
  entry = symbot;
  after_entry = 0;
  while (entry != NULL) {
    if (entry->attr & STATE) {
      if (after_entry)
	printf(",\n");
      printf("%s", entry->name);
      after_entry = 1;
    }
    entry = entry->next;
  }
  printf("\n]:\n");

  /* controls must be added by hand */
  printf("controls := [\n]:\n");

  printf("\n");
}

static void dot2score(char *name)
{
  while (*name != 0) {
    if (*name == '.' || *name == ':')
      *name = '_';
    name ++;
  }
}

static void mname_store(char *name)
{
  /*
   * -- determine mnamelen, including '.'
   * -- store model name
   */

  mnamelen = strlen(name) + 1;
  mname = (char *)malloc(mnamelen + 1);
  strcpy(mname, name);
  mname[mnamelen - 1] = '.';
  mname[mnamelen] = 0;
}

static void symexp_store(char *name)
{
  /*
   * -- store expansion for OHM names, add a '.'
   */

  symexplen = strlen(name) - 2 - mnamelen; /* without "{}" */
  if (symexplen < 0) {
    symexp[0] = 0;
    symexplen = 0;
    return;
  }
  name += 1 + mnamelen;
  if (symexplen + 2 > symexpsize) {
    symexpsize = symexplen + 2;
    symexp = (char *)realloc(symexp, symexpsize);
  }
  strcpy(symexp, name);
  symexp[symexplen++] = '.';
  symexp[symexplen] = 0;
}

static SymEntry *symlook(char *name, SymType type, SymAttr attr)
{
  SymEntry *entry, *insert;
  int maxsize, len;
  int input, cmp;

  /*
   * -- convert name into a valid symbol
   * -- look symbol up in symbol table
   * -- add a new entry, if not found
   */

  maxsize = symexplen + strlen(name) + 1;
  if (maxsize > symbolsize) {
    symbolsize = maxsize;
    symbol = (char *)realloc(symbol, symbolsize);
  }

  switch (type) {

  case VARIABLE:
    /* manage symbol expansion */
    if (symexplen != 0 && strncmp(name, mname, mnamelen) != 0) {
      /* add expansion */
      strcpy(symbol, symexp);
      strcpy(symbol + symexplen, name);
    }
    else if (strncmp(name, mname, mnamelen) == 0) {
      /* skip model name */
      name += mnamelen;
      strcpy(symbol, name);
    }
    else
      strcpy(symbol, name);
    break;

  case INPUT:
    len = strlen(symentry->name);
    if (strcmp(symentry->name + len - 8, "_control") == 0) {
      sprintf(symbol, "%s", symentry->name);
      symbol[len - 8] = '\0';
    }
    else if (strcmp(symentry->name + len - 6, "_input") == 0) {
      sprintf(symbol, "%s", symentry->name);
      symbol[len - 6] = '\0';
    }
    else {
      /* convert "*continuousinput*(.., ..., <nr>, time())" into "_i<nr>" */
      while (*name ++ != ',');
      while (*name ++ != ',');
      sscanf(name, "%d,", &input);
      sprintf(symbol, "_i%d", input);
    }
    break;

  case NUMBER:
    /* convert float constant into _f<nr>" */
    entry = symbot;
    while (entry != NULL) {
      if (entry->type == NUMBER && strcmp(name, entry->value) == 0) {
	strcpy(symbol, entry->name);
	break;
      }
      entry = entry->next;
    }
    if (entry == NULL) {
      sprintf(symbol, "_f%d", ++nnumbers);
    }
    break;

  default:
    strcpy(symbol, name);
    break;
  }

  if (type != NUMBER)
    dot2score(symbol);

  entry = symbot;
  insert = NULL;
  while (entry != NULL) {
    cmp = strcmp(symbol, entry->name);
    if (cmp == 0 && entry->type == type) {
      if (attr & STATE && !(entry->attr & STATE)) {
        entry->number = nstates++;
      }
      entry->attr |= attr;
      return entry;
    }
    if (cmp > 0)
      insert = entry;
    entry = entry->next;
  }

  /*
   * allocate and initialize a new entry
   */
  entry = (SymEntry *)malloc(sizeof(SymEntry));
  if (insert == NULL) {
    /* go on the bottom */
    entry->prev = NULL;
    entry->next = symbot;
    symbot = entry;
  }
  else {
    /* go behind insert */
    entry->prev = insert;
    entry->next = insert->next;
    insert->next = entry;
  }
  entry->name = (char *)malloc(strlen(symbol) + 1);
  strcpy(entry->name, symbol);
  entry->type = type;
  entry->attr = attr;
  entry->number = -1;
  if (type == NUMBER) {
    entry->value = (char *)malloc(strlen(name) + 1);
    strcpy(entry->value, name);
  }
  else
    entry->value = NULL;

  if (attr & STATE) {
    entry->number = nstates++;
  }
  else if (type == INPUT && attr & CONTROL) {
    entry->number = ncontrols++;
  }

  return entry;
}
