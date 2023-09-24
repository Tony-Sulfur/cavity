Program cavity

!======================================================================!
!   date of first version: 1985
!   date of this version: 2023/09/12
!======================================================================!


Implicit Real (A, B, D-H, J-Z), Complex (C)
Dimension axmn(9, 8), cw(20)
Dimension zmark(100), rwl(100), rwr(100), rhol(100), rhor(100)
Dimension az(10001), arw(10001), anrho(10001), afamp(10001), afphse(10001)
Common /const/ci, pi, realc
Common /ckt/zmark, rwl, rwr, rhol, rhor, izmark, izstep
Common /mode/wr, wi, fr0, fi0, xmn, im
Common /diag/az, arw, anrho, afamp, afphse, idiag, icont
External cbc
axmn = reshape([3.832, 1.841, 3.054, 4.201, 5.318, 6.416, 7.501, 8.578, 9.647,&
  7.016, 5.331, 6.706, 8.015, 9.282, 10.520, 11.735, 12.932, 14.116,&
  10.174,  8.536,  9.970, 11.346, 12.682, 13.987, 15.268, 16.529, 17.774,&
  13.324, 11.706, 13.170, 14.586, 15.964, 17.313, 18.637, 19.942, 21.229,&
  16.471, 14.864, 16.348, 17.789, 19.196, 20.576, 21.932, 23.268, 24.587,&
  19.616, 18.016, 19.513, 20.973, 22.401, 23.804, 25.184, 26.545, 27.889,&
  22.760, 21.164, 22.672, 24.145, 25.590, 27.010, 28.410, 29.791, 31.155,&
  25.904, 24.311, 25.826, 27.310, 28.768, 30.203, 31.618, 33.015, 34.397], [9, 8])
! the n-th root of jm'(x)=0

! universal constants
  ci = cmplx(0.0, 1.0)
  pi = 3.1415926
  realc = 2.99792E10

! No. of points to be used to mark the positions of boundaries
  izmark = 5

  r = 0.9   !  the mainbody radius    
  r1 = 0.5  !  the cutoff raduis(right)
  r2 = 1.1  !  the output radius(left)
  l = 11.7  !  the mainbody length 
  l1 = 3.0  !  the cutoff length(right) 
  l2 = 15.0  !  the outpput length(right)
  theta = 10.0  !  the taper angle
  
  zmark(1) = 0.0
  zmark(2) = zmark(1) + l1
  zmark(3) = zmark(2) + l
  zmark(4) = zmark(3) + (r2-r)/tan(theta*pi/180.0)
  zmark(5) = zmark(4) + l2
  
  rwl(1) = r1
  rwr(1) = r1
  rwl(2) = r1
  rwr(2) = r
  rwl(3) = r
  rwr(3) = r
  rwl(4) = r2
  rwr(4) = r2
  rwl(5) = r2
  rwr(5) = r2

! Wall resistivity arrays rhol and rhor in ohm-m.
  rhocu = 1.72E-8
  rhomks = 0.0  ! defalut value:0
  rhol(1) = rhomks
  rhor(1) = rhomks
  rhol(2) = rhomks
  rhor(2) = rhomks
  rhol(3) = rhomks
  rhor(3) = rhomks
  rhol(4) = rhomks
  rhor(4) = rhomks
  rhol(5) = rhomks
  rhor(5) = rhomks

! (for a new problem, always check convergence with respect to izstep).
  izstep = 5000

  ! specify m and n of te(m,n,l) mode
  im = 1
  in = 1
  xmn = axmn(iabs(im)+1, in)

! guess resonant frequency and q of the l-th axial mode
    
    ilmode = 1  !  choose the axianl mode you want to find
    wguess = realc*sqrt((ilmode*pi/l)**2+(xmn/r)**2)
    qguess = 500
    cw(1) = wguess*cmplx(1.0, -0.5/qguess)

! specify the field at the left boundary (i.e. at z=zmark(1))
    fr0 = 1.0E-4
    fi0 = 0.0

!  prepare the parameter to calculate w and q by muller 
    ep1 = 1.0E-6 !tolerance(ineration is stopped if rts<ep1)
    ep2 = ep1
    iroot = 1   ! the no. of root to find
    imaxit = 50 ! maximum no. of function evaluations allowed per zero
    icont = 0   ! the no. of searching time
    idiag = 0   ! instruction for diagnostics
    Call muller(0, iroot, cw, imaxit, ep1, ep2, cbc, .False.) 
    idiag = 1   ! if idiag = 1, cbc would store the profile, resistivity
    value = cabs(cbc(cw(1))) ! the accuracy of the root

    wr = cw(1)
    wi = -ci*cw(1)
    fhz = wr/2.0/pi
    q = -wr/(2.0*wi)

  ! calculate the ratio of q to qref
    lamda = realc/10.396e9
    qref = 4.0*pi*(l/lamda)**2/float(ilmode)
    qratio = q/qref

    

!  write results(fHz, q, afamp(amplitude), afphse(phase), value
  effectl1 = ilmode*pi/((2*pi*fhz/realc)**2 - (xmn/r)**2)**0.5


  Write (*, *) wr, wi, icont, value, fhz, q, qratio, zmark(3)+(12.67-11.7)
  
  

  
End Program cavity
