program main
  implicit none
  integer :: ifreq, iTx, iRx, nfreq, nTx, nRx
  real, dimension(:), allocatable :: freqs
  real, dimension(:), allocatable :: y_Tx, z_Tx, y_Rx, z_Rx
  character(256) :: num, filename
  real :: xx,yy,zr,xr,yr
  integer :: ios, isreceiver,nrec
  real :: x1(10000),x2(10000),x3(10000)
  
  filename = 'receivers.txt'
  open(15,file=filename)
  read(15,*)
  ios=0
  nrec=0
  do while (ios==0)   
     read(15,*,iostat=ios,end=999) xr,yr,zr,xx,yy
     nrec = nrec + 1
     x1(nrec)=xr
     x2(nrec)=yr
     x3(nrec)=zr
  end do
999 close(15)

  print *, 'nrec=', nrec

  nfreq = 3
  nTx = 1
  nRx = nrec
  
  allocate(freqs(nfreq))
  allocate(y_Tx(nTx))
  allocate(z_Tx(nTx))
  allocate(y_Rx(nRx))
  allocate(z_Rx(nRx))

  filename='sources.txt'
  open(15,file=filename)
  read(15,*)
  read(15,*) xr,yr,zr,xx,yy
  close(15)
  y_Tx(1) = 0.
  z_Tx(1) = zr !set the source depth
  
  y_Rx = x1(1:nrec)
  z_Rx = x3(1:nrec)

  freqs = [0.25, 0.75, 1.25]
  
  filename = 'test.emdata'
  open(10,file=filename,status='replace')

  write(10,*) 'Format:  EMData_2.3'
  write(10,*) 'Phase Convention: lag'
  write(10,*) '# CSEM Frequencies:', nfreq
  do ifreq = 1,nfreq
     write(10,*) freqs(ifreq)
  enddo
  write(10,*) 'Reciprocity used: no'  
  write(10,*) '# Transmitters:', nTx
  do iTx = 1,nTx
     write(num,'(i3.3)') iTx
     write(10,*) 0., y_Tx(iTx), z_Tx(iTx), 90., 0., 0., 0., 'edipole  ', 'Tx'//trim(adjustl(num))
  enddo
  write(10,*) '# CSEM Receivers:', nRx
  do iRx = 1,nRx
     write(num,'(i3.3)') iRx
     write(10,*) 0., y_Rx(iRx), z_Rx(iRx), 0., 0., 0., 0., 0., 'Rx'//trim(adjustl(num))
  enddo
  write(10,*) '# Data:', 2*nfreq*nTx*nRx
  do ifreq=1,nfreq
     do iTx=1,nTx
        do iRx=1,nRx
           write(10,*) 3, ifreq, iTx, iRx, 0, 0
           write(10,*) 4, ifreq, iTx, iRx, 0, 0
        enddo
     enddo
  enddo

  print *, 'test.emdata generated!'
  
  close(10)
end program main
