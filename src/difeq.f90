
Subroutine difeq(y, dy, ieqfst, ieqlst) 
    Implicit Real (A, B, D-H, J-Z), Complex (C)
    Dimension y(20), dy(20)
    Common /const/ci, pi, realc
    Common /mode/wr, wi, fr0, fi0, xmn, im
    cw = cmplx(wr, wi)
    z = y(5)
    rw = radius(z)
    ckz2 = (cw/realc)**2 - xmn**2*closs(z)/rw**2
    kz2r = ckz2
    kz2i = -ci*ckz2
    dy(1) = y(3)
    dy(2) = y(4)
    dy(3) = -kz2r*y(1) + kz2i*y(2)
    dy(4) = -kz2i*y(1) - kz2r*y(2)
    dy(5) = 1.0
    Return
  End Subroutine difeq