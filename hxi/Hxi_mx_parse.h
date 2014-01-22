/**
 * @file Hxi_mx_parse.h
 *   Parse an argument string passed to a Simulink S-function and
 *   create according mxArrays.
 *
 * (Simulink is a registered trademarks of The MathWorks, Inc.)
 *
 * rf, 10/19/2002
 */

/*
    Copyright (C) 1994--2005  Ruediger Franke

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; 
    version 2 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library (file COPYING.LIB);
    if not, write to the Free Software Foundation, Inc.,
    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#if !defined(Hxi_mx_parse_H)
#define Hxi_mx_parse_H

// include mxArray data type definitions from Matlab or Hxi
#include "simstruc.h"

#include <Meschach.h> // sscan_real

/** Hqp eXternal Interfaces */
namespace Hxi {

/**
 * @name Functions for parsing S-function argument strings
 */
//@{

/**
 * Return a pointer to the first non-whitespace character at or
 * behind start.
 * Note: newline '\n' is considered as an ordinary whitespace character,
 *       i.e. a semicolon ';' is required to delimit rows of a matrix.
 */
static const char *mx_forward_whitespaces(const char *start)
{
  const char *str;
  for (str = start; *str == ' ' || *str == '\t' || *str == '\n'; str++);
  return str;
}

/**
 * Return a pointer to the first non-whitespace character behind
 * the current argument.
 */
static const char *mx_forward_argument(const char *start)
{
  const char *str = mx_forward_whitespaces(start);
  int nparenthesis = 0;
  bool within_string = false;
  bool ready = false;
  for (; !ready; str++) {
    switch (*str) {

      // end of string
    case '\0':
      ready = true;
      break;

      // delimiter between arguments
    case ' ':
    case '\t':
    case '\n':
    case ',':
    case ';':
      if (nparenthesis == 0 && !within_string)
	ready = true;
      break;

      // parenthesis
    case '[':
    case '{':
      ++nparenthesis;
      break;

    case ']':
    case '}':
      if (nparenthesis > 0)
	--nparenthesis;
      else {
	ready = true;
      }
      break;

      // strings possibly containing special characters
    case '\'':
      if (!within_string)
	within_string = true;
      else
	within_string = false;
      break;

    default:
      break;
    }
  }
  return mx_forward_whitespaces(--str);
}

/**
 * Count numbers of columns n of argument row.
 * Return a pointer to the first non-whitespace character behind the row.
 * Note: Comma ',' or whitespaces are valid column delimiters.
 */
static const char *mx_count_columns(const char *row, int &n)
{
  const char *str = mx_forward_whitespaces(row);
  for (n = 0;
       *str != ';' && *str != ']'  && *str != '}' && *str != '\0';
       n++) {
    str = mx_forward_argument(str);
    if (*str == ',')
      str = mx_forward_whitespaces(++str);
  }
  return str;
}

/**
 * Count numbers of rows m and columns n of argument arg.
 * Return a pointer to the first non-whitespace character behind the argument.
 * Note: Brackets '[' and ']' are array delimiters.
 *       Semicolon ';' is required as row delimiter.
 */
static const char *mx_count_dimensions(const char *arg, int &m, int &n)
{
  const char *str = mx_forward_whitespaces(arg);
  if (*str != '[') {
    m = 1; // scalar value
    n = 1;
  } else {
    str = mx_forward_whitespaces(++str); // step into array
    n = 0;
    for (m = 0; *str != ']' && *str != '\0'; m++) {
      str = mx_count_columns(str, n);
      if (*str == ';')
	str = mx_forward_whitespaces(++str);
    }
  }
  return str;
}

/**
 * Parse one argument and return a dynamically allocated mxArray.
 * The function recognizes 
 *   - strings delimited by '\'' or '"'
 *   - scalar doubles
 *   - matrices of doubles
 *   - cell arrays enclosed by '{' '}' are parsed as strings
 */
static mxArray *mx_parse_argument(SimStruct *S, const char *arg)
{
  const char *str = mx_forward_whitespaces(arg);
  const char *str1;
  int m, n, i, j;
  mxArray *mxa = NULL;

  if (*str == '\'' || *str == '"') {
    char *arg_str, c, c_old;
    int idx;
    // obtain end of string argument
    str1 = mx_forward_argument(str);
    str++; // step into string
    while (*str1 != '\'' && *str1 != '"' && str1 > str)
      str1--;
    // create string
    // note: we need to check for quotes on quotes
    n = str1 - str;
    arg_str = (char *)malloc(n + 1);
    c_old = '\0';
    for (idx = 0, j = 0; j < n; j++) {
      c = *(str + j);
      if (!(c == '\'' && c_old == '\'')) {
	arg_str[idx++] = c;
	c_old = c;
      }
      else {
	c_old = '\0';
      }
    }
    arg_str[idx] = '\0';
    mxa = hmxCreateString(S, arg_str);
    free(arg_str);
  }
  else if (*str == '{') {
    // store cell array without parsing in an mxString
    char *arg_str;
    // obtain end of string argument
    str1 = mx_forward_argument(str);
    while (*str1 != '}' && str1 > str)
      str1--;
    str1++;  // step behind cell array
    // create string
    n = str1 - str;
    arg_str = (char *)malloc(n + 1);
    for (j = 0; j < n; j++) {
      arg_str[j] = *(str + j);
    }
    arg_str[j] = '\0';
    mxa = hmxCreateString(S, arg_str);
    free(arg_str);
  }
  else {
    // create mxDoubleMatrix
    mx_count_dimensions(str, m, n);
    mxa = hmxCreateDoubleMatrix(S, m, n, mxREAL);
    if (*str == '[')
      str = mx_forward_whitespaces(++str); // step into array
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        mxGetPr(mxa)[j*m+i] = sscan_real(str);
	str = mx_forward_argument(str);
	if (*str == ',')
	  str = mx_forward_whitespaces(++str); // step over col delimiter
      }
      if (*str == ';')
	str = mx_forward_whitespaces(++str); // step over row delimiter
    }
  }

  return mxa;
}

/** Print mxArray to FILE. */
static void mx_foutput(FILE *fp, mxArray *mx_arg)
{
  if (mxIsDouble(mx_arg)) {
    int i, j;
    int m = mxGetM(mx_arg);
    int n = mxGetN(mx_arg);
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
	fprintf(fp, "%g\t", mxGetPr(mx_arg)[j*m+i]);
      }
      fprintf(fp, "\n");
    }
  }
  else if (mxIsChar(mx_arg)) {
    char *str = mxArrayToString(mx_arg);
    fprintf(fp, "%s\n", str);
    mxFree(str);
  }
  else {
    fprintf(fp, "Error: unknown mxArray data type!\n");
  }
}

//@}

}; // namespace Hxi

#endif
