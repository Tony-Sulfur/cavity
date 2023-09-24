Subroutine rkint(derivy, y, dy, q, neqfst, neqlst, dx)
    Real y(neqlst), dy(neqlst), q(neqlst)
    Real a(4), b(4), c(4)
    Real dx, t
    External derivy
    Data a/0.5E0, 0.29289322E0, 1.7071068E0, 0.16666667E0/
    Data b/2.0E0, 1.0E0, 1.0E0, 2.0E0/
    Data c/0.5E0, 0.29289322E0, 1.7071068E0, 0.5E0/
    
    first: do  j = 1, 4
      Call derivy(y, dy, neqfst, neqlst)
      second: do i = neqfst, neqlst
        t = a(j)*(dy(i)-b(j)*q(i))
        y(i) = y(i) + dx*t
        q(i) = q(i) + 3.0E0*t - c(j)*dy(i)
      End Do second
    End Do first
    Return
  End Subroutine rkint