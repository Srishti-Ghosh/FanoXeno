      program test
c    Test program for CG coef.
    1 continue
      read (5,*) AJ1,AJ2,AJ3,AM1,AM2,AM3
      R = CF(AJ1,AJ2,AJ3,AM1,AM2 ,AM3)
      write ( 6,*) AJ1,AJ2,AJ3,AM1,AM2,AM3,R
      go to 1
    2 continue
      stop
      end  
      FUNCTION CF(AJ1,AJ2,AJ3,AM1,AM2,AM3)                              ACYK0245
C    FUNCTION PROGRAM FOR THE C.G. COEFFICIENT USING THE HYPERGEOMETRIC ACYK0246
C    FUNCTION FORMULA GIVEN BY (2.8) OF TEXT.                           ACYK0247
      dimension FCT(500)                                                ACYK0248
C    ALL THE 500 FACTORIALS ARE COMPUTED FIRST AND SET UP AS AN ARRAY.  ACYK0019
      FCT(1) = 0.0                                                      ACYK0020
      FCT(2) = 0.0                                                      ACYK0021
      DO 10 N = 3,500                                                   ACYK0022
      FF = N-1                                                          ACYK0023
      FCT(N) = FCT(N-1) + ALOG(FF)                                      ACYK0024
   10 continue   
      CF = 0.0                                                          ACYK0249
      J = AJ1 + AJ2 + AJ3                                               ACYK0250
C    THE TRIANGULAR INEQUALITY IS CHECKED FIRST.                        ACYK0251
      CHK = TRIA(AJ1,AJ2,AJ3)                                           ACYK0252
      IF(CHK.EQ.0.0.OR.((AM1+AM2).NE.AM3)) RETURN                       ACYK0253
      IF((AJ1.EQ.0.0.AND.AM1.EQ.0.0).OR.(AJ2.EQ.0.0.AND.AM2.EQ.0.0))    ACYK0254
     1 GO TO 120                                                        ACYK0255
      IF(AJ3.EQ.0.0.AND.AM3.EQ.0.0) GO TO 130                           ACYK0256
      JX = -AJ1 + AJ2 + AJ3                                             ACYK0257
      JY = AJ1 - AJ2 + AJ3                                              ACYK0258
      JZ = AJ1 + AJ2 - AJ3                                              ACYK0259
      IF(AM1.EQ.0.0.AND.AM2.EQ.0.0) GO TO 110                           ACYK0260
      JJM = AJ1 - AJ2 + AM3                                             ACYK0261
      N1 = AJ1 + AM1                                                    ACYK0262
      N2 = AJ1 - AM1                                                    ACYK0263
      N3 = AJ2 + AM2                                                    ACYK0264
      N4 = AJ2 - AM2                                                    ACYK0265
      N5 = AJ3 - AM3                                                    ACYK0266
      N6 = AJ3 + AM3                                                    ACYK0267
C    THE POSITIVE NATURE OF THE DENOMINATOR PARAMETERS IS CHECKED FOR   ACYK0268
C    THE SELECTION OF THE VALID F32 FUNCTION.                           ACYK0269
      IF(N5.GE.N2.AND.N6.GE.N3) GO TO 20                                ACYK0270
      IF(N1.GE.N4.AND.N2.GE.N5) GO TO 30                                ACYK0271
      IF(N3.GE.N6.AND.N4.GE.N1) GO TO 40                                ACYK0272
      IF(N5.GE.N4.AND.N6.GE.N1) GO TO 50                                ACYK0273
      IF(N1.GE.N6.AND.N2.GE.N3) GO TO 60                                ACYK0274
      IF(N5.GE.N2.AND.N4.GE.N5) GO TO 70                                ACYK0275
   20 IA = N2                                                           ACYK0276
      IB = N3                                                           ACYK0277
      IC = JZ                                                           ACYK0278
      ID = 1 + N5 - N2                                                  ACYK0279
      IE = 1 + N6 - N3                                                  ACYK0280
      M = N1 - N4                                                       ACYK0281
      GO TO 80                                                          ACYK0282
   30 IA = N4                                                           ACYK0283
      IB = N5                                                           ACYK0284
      IC = JX                                                           ACYK0285
      ID = 1 + N1 - N4                                                  ACYK0286
      IE = 1 + N2 - N5                                                  ACYK0287
      M = N3 - N6                                                       ACYK0288
      GO TO 80                                                          ACYK0289
   40 IA = N6                                                           ACYK0290
      IB = N1                                                           ACYK0291
      IC = JY                                                           ACYK0292
      ID = 1 + N3 - N6                                                  ACYK0293
      IE = 1 + N4 - N1                                                  ACYK0294
      M = N5 - N2                                                       ACYK0295
      GO TO 80                                                          ACYK0296
   50 IA = N4                                                           ACYK0297
      IB = N1                                                           ACYK0298
      IC = JZ                                                           ACYK0299
      ID = 1 + N5 - N4                                                  ACYK0300
      IE = 1 + N6 - N1                                                  ACYK0301
      M = J + N3 - N2                                                   ACYK0302
      GO TO 80                                                          ACYK0303
   60 IA = N6                                                           ACYK0304
      IB = N3                                                           ACYK0305
      IC = JX                                                           ACYK0306
      ID = 1 + N1 - N6                                                  ACYK0307
      IE = 1 + N2 - N3                                                  ACYK0308
      M = J + N5 - N4                                                   ACYK0309
      GO TO 80                                                          ACYK0310
   70 IA = N2                                                           ACYK0311
      IB = N5                                                           ACYK0312
      IC = JY                                                           ACYK0313
      ID = 1 + N5 - N2                                                  ACYK0314
      IE = 1 + N4 - N5                                                  ACYK0315
      M = J + N1 - N6                                                   ACYK0316
   80 F32 = 1.0                                                         ACYK0317
      CONST =FCT(N1+1)+FCT(N2+1)+FCT(N3+1)+FCT(N4+1)+FCT(N5+1)+FCT(N6+1)ACYK0318
     1     +FCT(JX+1)+FCT(JY+1)+FCT(JZ+1)-FCT(J+2)-2.0*(FCT(IA+1)       ACYK0319
     2     +FCT(IB+1)+FCT(IC+1)+FCT(ID)+FCT(IE))                        ACYK0320
C    IF ANY ONE OF THE NUMERATOR PARAMETERS IS ZERO, THEN THE VALUE OF  ACYK0321
C    THE F32 IS SET EQUAL TO 1.0 AND THE PROGRAM SEGMENT FOR F32 IS SKIPACYK0322
      IF ( IA.EQ.0.OR.IB.EQ.0.OR.IC.EQ.0) GO TO 100                     ACYK0323
      N = MIN0(IA,IB,IC)                                                ACYK0324
      IX = N - 1                                                        ACYK0325
   90 RN = (IX-IA)*(IX-IB)*(IX-IC)                                      ACYK0326
      RND = (IX+ID)*(IX+IE)*(IX+1)                                      ACYK0327
      F32 = 1.0 + RN*F32/RND                                            ACYK0328
      IX = IX - 1                                                       ACYK0329
      IF(IX.GE.0) GO TO 90                                              ACYK0330
  100 CF = PHASE(JJM+M)*SQRT((2.0*AJ3+1.0)*EXP(CONST))*F32              ACYK0331
      RETURN                                                            ACYK0332
C    SPECIAL VALUES OF THE CLEBSCH-GORDAN COEFFICIENT GIVEN BY (2.14) OFACYK0333
  110 IF(PHASE(J).EQ.(-1.0)) RETURN                                     ACYK0334
      IF(AJ1.EQ.0.0.OR.AJ2.EQ.0.0) GO TO 120                            ACYK0335
      IF(AJ3.EQ.0.0) GO TO 130                                          ACYK0336
      N10 = J/2                                                         ACYK0337
      L1 = AJ1                                                          ACYK0338
      L2 = AJ2                                                          ACYK0339
      L3 = AJ3                                                          ACYK0340
      N11 = (L1 + L2 - L3)/2                                            ACYK0341
      CF = PHASE(N11)*SQRT((2.0*AJ3+1.0)*EXP(FCT(JX+1)+FCT(JY+1)        ACYK0342
     1     +FCT(JZ+1)-FCT(J+2)+2.0*(FCT(N10+1)-FCT(N10-L1+1)            ACYK0343
     2     -FCT(N10-L2+1)-FCT(N10-L3+1))))                              ACYK0344
      RETURN                                                            ACYK0345
  120 CF = 1.0                                                          ACYK0346
      RETURN                                                            ACYK0347
  130 N = AJ1 - AM1                                                     ACYK0348
      CF = PHASE(N)/SQRT(2.0*AJ1+1.0)                                   ACYK0349
      RETURN                                                            ACYK0350
      END                                                               ACYK0351
      FUNCTION PHASE(N)                                                 ACYK0115
C    FUNCTION PHASE(N) IS FOR COMPUTING (-1)**N.                        ACYK0116
      PHASE = 1.0                                                       ACYK0117
      M = (N/2)*2                                                       ACYK0118
      IF(M.NE.N) PHASE = -1.0                                           ACYK0119
      RETURN                                                            ACYK0120
      END                                                               ACYK0121
      FUNCTION TRIA(X,Y,Z)                                              ACYK0122
C    FUNCTION TRIA(X,Y,Z) IS FOR CHECKING THE TRIANGULAR INEQUALITY     ACYK0123
C    CONDITION ON ANY THREE ANGULAR MOMENTA.                            ACYK0124
      TRIA = 0.0                                                        ACYK0125
      IF((X+Y-Z).GE.0..AND.(X-Y+Z).GE.0..AND.(-X+Y+Z).GE.0.) TRIA = 1.0 ACYK0126
      RETURN                                                            ACYK0127
      END                                                               ACYK0128

