Subroutine muller(kn, n, rts, maxit, ep1, ep2, fn, fnreal)
    Implicit Complex (A-H, O-Z)
    Complex num, lambda
    Real ep1, ep2, eps1, eps2
    Real abso
    External fn
    Real aimag
    Logical fnreal
    Dimension rts(6)
    aimag(x) = (0.E0, -1.E0)*x
  
    ! initialization.
    eps1 = amax1(ep1, 1.E-12)
    eps2 = amax1(ep2, 1.E-20)
  
    Do i = kn+1, kn+n
      kount = 0
    ! compute first three estimates for root as,
    !     rts(i) + 0.5 , rts(i) - 0.5 , rts(i).
      abso = cabs(rts(i))
      If (abso) 11, 12, 11
      11 firsss = rts(i)/100.0E0
      Goto 13
      12 firsss = cmplx(1.0E0, 0.0E0)
      13 Continue
      secdd = firsss/100.0E0
      1 h = .5E0*firsss
      rt = rts(i) + h
      Assign 10 To nn
      Goto 70
      10 delfpr = frtdef
      rt = rts(i) - h
      Assign 20 To nn
      Goto 70
      20 frtprv = frtdef
      delfpr = frtprv - delfpr
      rt = rts(i)
      Assign 30 To nn
      Goto 70
      30 Assign 80 To nn
      lambda = -0.5
    ! compute next estimate for root.
      40 delf = frtdef - frtprv
      dfprlm = delfpr*lambda
      num = -frtdef*(1.0+lambda)*2
      g = (1.0+lambda*2)*delf - lambda*dfprlm
      sqr = g*g + 2.0*num*lambda*(delf-dfprlm)
      If (fnreal .And. real(sqr)<0.0) sqr = 0.0
      sqr = csqrt(sqr)
      den = g + sqr
      If (real(g)*real(sqr)+aimag(g)*aimag(sqr)<0.0) den = g - sqr
      If (cabs(den)==0.0) den = 1.0
      lambda = num/den
      frtprv = frtdef
      delfpr = delf
      h = h*lambda
      rt = rt + h
      If (kount>maxit) Write (*, 1492)
      If (kount>maxit) Goto 100
    !
      70 kount = kount + 1
      frt = fn(rt)
      frtdef = frt
      If (i<2) Goto 75
    !
      Do j = 2, i
        den = rt - rts(j-1)
        If (cabs(den)<eps2) Goto 79
        frtdef = frtdef/den
      End Do
      75 Goto nn, (10, 20, 30, 80)
      79 rts(i) = rt + secdd
      Goto 1
    ! check for convergence.
      80 If (cabs(h)<eps1*cabs(rt)) Goto 100
      If (amax1(cabs(frt),cabs(frtdef))<eps2) Goto 100
    !
    ! check for divergence.
      If (cabs(frtdef)<10.0*cabs(frtprv)) Goto 40
      h = h/2.0
      lambda = lambda/2.0
      rt = rt - h
      Goto 70
      100 rts(i) = rt
    End Do
    Return
    1492 Format (' lack of convergence in muller')
  End Subroutine muller