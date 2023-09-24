! rho(z) is the wall resistivity at position z.
Function rho(z)
    Implicit Real (A, B, D-H, J-Z), Complex (C)
    Dimension zmark(100), rwl(100), rwr(100), rhol(100), rhor(100)
    Common /ckt/zmark, rwl, rwr, rhol, rhor, izmark, izstep
    imax = izmark - 1
    Do i = 1, imax
      If (z>=zmark(i) .And. z<=zmark(i+1)) then
        rho = rhor(i) + (rhol(i+1)-rhor(i))*(z-zmark(i))/(zmark(i+1)-zmark(i))
      end if 
    End Do
  
    If (z<zmark(1)) rho = rhol(1)
    If (z>zmark(izmark)) rho = rhor(izmark)
    Return
  End Function rho
  