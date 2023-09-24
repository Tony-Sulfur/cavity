! Called by muller to check D(w)=0


Function cbc(cw)
    Implicit Real (A, B, D-H, J-Z), Complex (C)
    Dimension y(20), dy(20), q(20)
    Dimension zmark(100), rwl(100), rwr(100), rhol(100), rhor(100)
    Dimension az(10001), arw(10001), anrho(10001), afamp(10001), afphse(10001)
    Common /const/ci, pi, realc
    Common /ckt/zmark, rwl, rwr, rhol, rhor, izmark, izstep
    Common /mode/wr, wi, fr0, fi0, xmn, im
    Common /diag/az, arw, anrho, afamp, afphse, idiag, icont
    External difeq
    wr = cw
    wi = -ci*cw
  
    ! initialize rf field array y at left boundary.
    z1 = zmark(1)
    rw = radius(z1)
    If ((wr/realc)**2 >= xmn**2/rw**2) then ! aobve(equal) the cutoff
      ckz2 = (cw/realc)**2 - xmn**2*closs(z1)/rw**2
      ckz = csqrt(ckz2)
      kzr = ckz
      kzi = -ci*ckz
      y(1) = fr0
      y(2) = fi0
      y(3) = kzr*y(2) + kzi*y(1)
      y(4) = -kzr*y(1) + kzi*y(2)
      y(5) = z1
    else  ! below the cutoff
      ckapa2 = xmn**2*closs(z1)/rw**2 - (cw/realc)**2
      ckapa = csqrt(ckapa2)
      kapar = ckapa
      kapai = -ci*ckapa
      y(1) = fr0
      y(2) = fi0
      y(3) = kapar*y(1) - kapai*y(2)
      y(4) = kapai*y(1) + kapar*y(2)
      y(5) = z1
    end if 
    
    az(1) = z1
    arw(1) = radius(z1) ! store radius
    anrho(1) = rho(z1)/1.72E-8 ! store the resistivity(normalized)
    afamp(1) = sqrt(y(1)**2+y(2)**2) ! store the rf field
    afphse(1) = atan2(y(2), y(1)) ! store the transit angle
    
    ! integrate deq for field array y from left boundary to right boundary.
    Do i = 1, 5
      q(i) = 0.0
    End Do
    l = zmark(izmark) - z1
    delz = l/float(izstep)
    Do iz = 1, izstep ! advance rf field array y by one step.
      Call rkint(difeq, y, dy, q, 1, 5, delz) 
      If (idiag ==1) then
        z = y(5)
        az(iz+1) = z
        arw(iz+1) = radius(z) ! store radius
        anrho(iz+1) = rho(z)/1.72E-8 ! store the resistivity(normalized)
        afamp(iz+1) = sqrt(y(1)**2+y(2)**2) ! store the rf field
        afphse(iz+1) = atan2(y(2), y(1)) ! store the transit angle
      end if 
    End Do
  
    ! calculate cbc(cw)
    fr = y(1)
    fi = y(2)
    fpr = y(3)
    fpi = y(4)
    cf = cmplx(fr, fi)
    cfp = cmplx(fpr, fpi)
    zf = y(5)
    rwf = radius(zf)
    If ((wr/realc)**2 >= xmn**2/rwf**2) then ! above(equal) cutoff
      ckz2 = (cw/realc)**2 - xmn**2*closs(zf)/rwf**2
      ckz = csqrt(ckz2)
      cbc = cfp - ci*ckz*cf  !check the D(w)=0
    else ! below cutoof
      ckapa2 = xmn**2*closs(zf)/rwf**2 - (cw/realc)**2
      ckapa = csqrt(ckapa2)
      cbc = cfp + ckapa*cf
    end if 
    icont = icont + 1
    Return
  End Function cbc