/* tpc.c --
 *    - build from Tcl scripts a c_file with a global string
 *    - usage: tpc [-o <c_file>] [-ts <max_token_size>] tcl_file ...
 *    - default c_file: tpc.out
 *
 * Copyright (c) 1994 R"udiger Franke
 * All Rights Reserved.
 * 
 * Redistribution and use in any form, with or without modification, 
 * is permitted, provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in other form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *       This product includes software developed by R"udiger Franke.
 * 4. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main (int argc, char *argv[])
{
  char *outputFile;
  int  tsize;
  FILE *in, *out;
  int  argnum;
  int  c, index, written;
  int  linestart;

  /*
   * usage
   */

  if (argc == 1) {
    printf ("usage:");
    printf ("  tpc [-o <c_file>] [-ts <max_token_size>] tcl_file ...\n");
    return 0;
  }

  /*
   *  parse command line options
   *  (don't use getopt() to get rid of unistd.h)
   */

  outputFile = "tpc.out";
  tsize = -1;

  argnum = 1;
  while (argc > argnum && argv[argnum][0] == '-') {
    if (strncmp(argv[argnum], "-o", 2) == 0) {
      if (strlen(argv[argnum]) > 2) {
	outputFile = &argv[argnum][2];
	argnum ++;
      }
      else {
        if (argc > argnum + 1) {
	  outputFile = argv[argnum+1];
	  argnum += 2;
	} else {
	  fprintf (stderr, "tpc: missing file name for -o\n");
	  return -1;
	}
      }
    }
    else if (strncmp(argv[argnum], "-ts", 3) == 0) {
      if (strlen(argv[argnum]) > 3) {
	tsize = atoi(&argv[argnum][3]);
	argnum ++;
      }
      else {
        if (argc > argnum + 1) {
	  tsize = atoi(argv[argnum+1]);
	  argnum += 2;
	} else {
	  fprintf (stderr, "tpc: missing token size for -ts\n");
	  return -1;
	}
      }
    }
    else {
      fprintf (stderr, "tpc: unknown option \"%s\".\n", argv[argnum]);
      return -1;
    }
  }

 if (argc == argnum) {
   fprintf(stderr, "tpc: no input files specified.\n");
   return -1;
 }

 /*
  *  open output file, write variable declaration
  */

  if ((out= fopen (outputFile, "w")) == NULL) {
    fprintf (stderr, "tpc: open output file %s", outputFile);
    perror ("");
    return -1;
  }
  
  if (tsize > 0)
    fprintf (out, "\nstatic char _tpc_array[][%d] = {\n\"", tsize);
  else
    fprintf (out, "\nstatic char _tpc_array[] = {\n\"");
  
  /*
   *  cat all source files onto output file
   */

  written = 0;
  linestart = 1;
  while (argnum < argc) {
    if( (in= fopen( argv[argnum], "r")) == NULL) {
      fprintf( stderr, "tpc: open file %s", argv[argnum]);
      perror( "");
      continue;
    }

    while (1) {

      c = getc (in);

      /* strip comments and leading whitespaces */

      if (linestart == 1) {
	switch (c) {
	  case '#':
	    while ((c = getc (in) != '\n') && c != EOF) {};
            break;

	  case ' ':
	  case '\t':
	  case '\n':
	  case EOF:
	    break;

	  default:
	    linestart = 0;
        }	  
      }

      /* end test */
    
      if (c == EOF)
	break;

      /* output read character */

      if (linestart == 0) {
        switch (c) {
          case '\t':
	    fprintf( out, "\\t");
	    written ++;
	    break;

          case '\n':
	    fprintf( out, "\\n");
	    written ++;
	    linestart = 1;
	    break;

          case '\\':
          case '\"':
	    fputc( '\\', out);
          default:
	    fputc( c, out);
	    written ++;
        }
      }

      /* start new token if current is full */

      if (tsize > 0)
	if (written == tsize) {
	  fprintf (out, "\",\n\"");
	  written = 0;
	}
    }

    fclose (in);
    argnum++;
  }

  fprintf (out, "\"\n};\n\n");

  /*
   * declare global variable for whole string
   */

  fprintf (out, "char *");

  /* get relative file name */
  index = strlen (outputFile) - 1;
  while (outputFile[index] != '/' && index -- > 0);
  index ++;
  written = 0;
  while ((c = outputFile[index++]) != '\0') {
    if( c == '.' || c == '-')
      break;
    fputc (c, out);
    written++;
  }

  if (written == 0) {
    fprintf (stderr, "tpc: can't build variable name from \"%s\".\n",
             outputFile);
    fclose (out);
    return -1;
  }
  if (tsize > 0)
    fprintf (out, " = _tpc_array[0];\n");
  else
    fprintf (out, " = _tpc_array;\n");

  fclose (out);
  return 0;
}

