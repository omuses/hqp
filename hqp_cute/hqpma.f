C
C  HQP driver for problems derived from SIF files.
C
C  R. Franke (based upon A. R. Conn's and Ph. Toint's cobma.f)
C  10/30/1995.
C
      PROGRAM HQPMA
C
      CHARACTER * 256   ARG1
C
      CALL GETARG(1, ARG1)
      CALL TCLMN(ARG1)
C
      END
C
      SUBROUTINE CSIZE(N, M, NNZJ, NNZH)
      INTEGER N, M, NNZJ, NNZH
C
      INTEGER NG, NNZJMX, NNZHMX
CTOY  PARAMETER      ( NNZJMX =   9000, NNZHMX =   9000 )
CMED  PARAMETER      ( NNZJMX =  90000, NNZHMX =  90000 )
CBIG  PARAMETER      ( NNZJMX = 900000, NNZHMX = 900000 )
C
      INTEGER          INPUT, IOUT
C
CSALF CHARACTER * 64   PRBDAT
CSALF PARAMETER      ( INPUT = 55, IOUT = 6 )
CUNIX CHARACTER * 64   PRBDAT
CUNIX PARAMETER      ( INPUT = 55, IOUT = 6 )
CVMS  CHARACTER * 64   PRBDAT
CVMS  PARAMETER      ( INPUT = 55, IOUT = 6 )
CWFC  CHARACTER * 64   PRBDAT
CWFC  PARAMETER      ( INPUT = 55, IOUT = 6 )
C
C  Build data input file name
C
CSALF PRBDAT = 'OUTSDIF.DAT'
CUNIX PRBDAT = 'OUTSDIF.d'
CVMS  PRBDAT = 'OUTSDIF.DAT'
CWFC  PRBDAT = 'OUTSDIF.DAT'
C
C  Open the input file.
C
CSALF OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED' )
CUNIX OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED', STATUS = 'OLD' )
CUNIX REWIND INPUT
CVMS  OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED',
CVMS *       STATUS = 'OLD', CARRIAGECONTROL = 'LIST' )
CWFC  OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED' )
C
C  Find out problem size.
C
      READ( INPUT, 1000 ) N , NG
      CLOSE ( INPUT )
C
      M = NG
      NNZJ = NNZJMX
      NNZH = NNZHMX
C
C  Non-executable statements.
C
 1000 FORMAT( 10I6 )
C
C     end of CSIZE
C
      END
C
C subroutine CINIT
C
      SUBROUTINE CINIT(N, M, X0, BL, BU, BINF,
     *                 EQUATN, LINEAR, V0, CL, CU,
     *                 EFIRST, LFIRST, NVFRST)
      INTEGER          N, M
CS    REAL             X0    ( N ), BL    ( N ), BU    ( N )
CS    REAL             V0    ( M ), CL    ( M ), CU    ( M )
CS    REAL             BINF
CD    DOUBLE PRECISION X0    ( N ), BL    ( N ), BU    ( N )
CD    DOUBLE PRECISION V0    ( M ), CL    ( M ), CU    ( M )
CD    DOUBLE PRECISION BINF
      LOGICAL          EQUATN( M ), LINEAR( M )
      LOGICAL          EFIRST, LFIRST, NVFRST
C
      INTEGER          INPUT, IOUT, NMAX, MMAX
C
CSALF CHARACTER * 64   PRBDAT
CSALF PARAMETER      ( INPUT = 55, IOUT = 6 )
CUNIX CHARACTER * 64   PRBDAT
CUNIX PARAMETER      ( INPUT = 55, IOUT = 6 )
CVMS  CHARACTER * 64   PRBDAT
CVMS  PARAMETER      ( INPUT = 55, IOUT = 6 )
CWFC  CHARACTER * 64   PRBDAT
CWFC  PARAMETER      ( INPUT = 55, IOUT = 6 )
CS    REAL             BIGINF
CD    DOUBLE PRECISION BIGINF
CS    PARAMETER      ( BIGINF = 9.0E+19 )
CD    PARAMETER      ( BIGINF = 9.0D+19 )
C
C  Build data input file name
C
CSALF PRBDAT = 'OUTSDIF.DAT'
CUNIX PRBDAT = 'OUTSDIF.d'
CVMS  PRBDAT = 'OUTSDIF.DAT'
CWFC  PRBDAT = 'OUTSDIF.DAT'
C
C  Open the input file.
C
CSALF OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED' )
CUNIX OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED', STATUS = 'OLD' )
CUNIX REWIND INPUT
CVMS  OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED',
CVMS *       STATUS = 'OLD', CARRIAGECONTROL = 'LIST' )
CWFC  OPEN ( INPUT, FILE = PRBDAT, FORM = 'FORMATTED' )
C
C     Call CSETUP
C
      NMAX = N
      MMAX = M
      CALL CSETUP(INPUT, IOUT, N, M, X0, BL, BU, NMAX,
     *            EQUATN, LINEAR, V0, CL, CU, MMAX, 
     *            EFIRST, LFIRST, NVFRST)
C
      CLOSE ( INPUT )
C      
      BINF = BIGINF
C
C End of CINIT
C
      END
C
C SUBROUTINE CWRTSN
C
      SUBROUTINE CWRTSN(N, M, HEADER, F, X, V)
      INTEGER N, M
      CHARACTER * 60 HEADER
CS    REAL             F, X ( N ), V ( M )
CD    DOUBLE PRECISION F, X ( N ), V ( M )
C
      INTEGER NMAX, MMAX
CTOY  PARAMETER      ( NMAX =   500, MMAX =   500 )
CMED  PARAMETER      ( NMAX =  5000, MMAX =  5000 )
CBIG  PARAMETER      ( NMAX = 50000, MMAX = 50000 )
C
      CHARACTER * 10 PNAME
      CHARACTER * 10 XNAMES ( NMAX )
      CHARACTER * 10 GNAMES ( MMAX )
C
      INTEGER          OUTPUT, IOUT, I
C
CSALF CHARACTER * 64   SLNDAT
CSALF PARAMETER      ( OUTPUT = 56, IOUT = 6 )
CUNIX CHARACTER * 64   SLNDAT
CUNIX PARAMETER      ( OUTPUT = 56, IOUT = 6 )
CVMS  CHARACTER * 64   SLNDAT
CVMS  PARAMETER      ( OUTPUT = 56, IOUT = 6 )
CWFC  CHARACTER * 64   SLNDAT
CWFC  PARAMETER      ( OUTPUT = 56, IOUT = 6 )
C
C  Build data output file name
C
CSALF SLNDAT = 'SOLUTION.DAT'
CUNIX SLNDAT = 'SOLUTION.d'
CVMS  SLNDAT = 'SOLUTION.DAT'
CWFC  SLNDAT = 'SOLUTION.DAT'
C
C Check memory size
C
      IF ( N .GT. NMAX ) THEN
         IF ( IOUT .GT. 0 ) THEN
            WRITE( IOUT, 2100 ) 'XNAMES', 'NMAX  ', N - NMAX
         END IF
         STOP
      END IF
      IF ( M .GT. MMAX ) THEN
         IF ( IOUT .GT. 0 ) THEN
            WRITE( IOUT, 2100 ) 'GNAMES', 'MMAX  ', M - MMAX
         END IF
         STOP
      END IF
C
C  Open the output file.
C
CSALF OPEN ( OUTPUT, FILE = SLNDAT, FORM = 'FORMATTED' )
CUNIX OPEN ( OUTPUT, FILE = SLNDAT, FORM = 'FORMATTED',
CUNIX*       STATUS = 'UNKNOWN' )
CUNIX REWIND OUTPUT
CVMS  OPEN ( OUTPUT, FILE = SLNDAT, STATUS = 'UNKNOWN',
CVMS *       CARRIAGECONTROL = 'LIST' )
CWFC  OPEN ( OUTPUT, FILE = SLNDAT, FORM = 'FORMATTED',
CWFC *       STATUS = 'UNKNOWN')
CWFC  REWIND OUTPUT
C
C  Get names
C
      CALL CNAMES(N, M, PNAME, XNAMES, GNAMES)
C
C  Write solution
C
      WRITE( OUTPUT, 2000 ) HEADER
      WRITE( OUTPUT, 2010 ) PNAME
      IF ( N .GT. 0 ) THEN
         WRITE( OUTPUT, 2020 )
         DO 1500 I = 1, N
            WRITE( OUTPUT, 2030 ) XNAMES( I ), X( I )
 1500    CONTINUE
         WRITE( OUTPUT, 2040 )
         DO 1510 I = 1, M
            WRITE( OUTPUT, 2030 ) GNAMES( I ), V( I )
 1510    CONTINUE
         WRITE( OUTPUT, 2060 )
         WRITE( OUTPUT, 2070 ) F
      ELSE
         WRITE( OUTPUT, 2090 )
      END IF
C
C  Close output file
C
      CLOSE ( OUTPUT )
C
C  Non-executable statements.
C
 2000 FORMAT( '*   ', A60)
 2010 FORMAT( '*   problem name: ', A10)
 2020 FORMAT( /, '*   variables', /)
 2030 FORMAT( '    SOLUTION  ', A10, 1P, D12.5 )
 2040 FORMAT( /, '*   Lagrange multipliers', /)
 2060 FORMAT( /, '*   objective function value', /)
 2070 FORMAT( ' XL SOLUTION  ', 10X, 1P, D12.5 )
 2090 FORMAT( /, '*   no solution')
 2100 FORMAT( /, ' ** SUBROUTINE CWRTSN: array length ', A6,
     *        ' too small.', /, ' -- Miminimization abandoned.',
     *        /, ' -- Increase the parameter ', A6, ' by at least ', I8,
     *           ' and restart.'  )
C
C     end of CWRTSN
C
      END
