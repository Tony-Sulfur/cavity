! closs(z) is a factor to account for wall losses of the te(m,n) mode
! at position z. closs(z)=1.0 for zero wall resistivity.
Function closs(z)
    Implicit Real (A, B, D-H, J-Z), Complex (C)
    Common /const/ci, pi, realc
    Common /mode/wr, wi, fr0, fi0, xmn, im
    cw = cmplx(wr, wi)
    delta = sqrt(realc**2*rho(z)/(2.0*pi*wr*9.0E9))
    rw = radius(z)
    wcmn = xmn*realc/rw
    cdum1 = (1.0+ci)*delta/rw
    cdum2 = (float(im**2)/(xmn**2-float(im**2)))*(cw/wcmn)**2
    closs = 1.0 - cdum1*(1.0+cdum2)
    Return
  End Function closs