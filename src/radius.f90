! radius(z) is the wall radius at position z.
Function radius(z)
    Implicit Real (A, B, D-H, J-Z), Complex (C)
    Dimension zmark(100), rwl(100), rwr(100), rhol(100), rhor(100)
    Common /ckt/zmark, rwl, rwr, rhol, rhor, izmark, izstep
    imax = izmark - 1
  
    Do i = 1, imax
      If (z>=zmark(i) .And. z<=zmark(i+1)) then 
      radius = rwr(i) + (rwl(i+1)-rwr(i))*(z-zmark(i))/(zmark(i+1)-zmark(i))
      end if 
    End Do
    If (z<zmark(1)) then
      radius = rwl(1)
    end if 
    if (z>zmark(izmark)) then
      radius = rwr(izmark)
    end if 
    Return
  End Function radius