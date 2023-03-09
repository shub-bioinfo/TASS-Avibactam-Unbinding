program histogram
implicit none
real*8, allocatable :: p(:)
real*8 :: grid_min, grid_max, grid_width, x, dum, norm, value
integer :: nbin, bin, i, steps

open(2,file='COLVAR',status='old')

! 1. Provide the grid min, grid max and the grid width.
! 2. If we have to give the number of MD steps, mention after the steps count 
!     suroutine, steps=number_of_steps.
! 3. Check for the Column which has the correspoding CV value in line 61.
! 4. Probability distribution generated in "output" file.


grid_width=0.034d0
grid_min=13.0d0
grid_max=47.0d0
!CALL grid_min_max(2,grid_min,grid_max)

nbin = NINT((grid_max-grid_min)/grid_width)+1

print*, 'grid_min, grid_max, nbin and grid width=', grid_min, grid_max, nbin,  grid_width

allocate(p(nbin))

!rewind(2)

call step_count(2,steps)
print *, 'steps=', steps

rewind(2)

!steps=1600

norm=0.d0
p(1:nbin)=0.d0
do i=1,steps
   read(2,*) dum,value!,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum
   bin = nint((value-grid_min)/grid_width)+1 
!   print*, i, steps
   if(bin.gt.0.and.bin.le.nbin) then
      p(bin) = p(bin)+ 1.d0
   end if
   norm = norm + 1.d0
end do

print*, 'norm =', norm
norm = norm*grid_width  
print*, 'norm*area =', norm

open(10,file='output',status='replace')
do i=1,nbin
 x= grid_min + float(i-1)*grid_width
write(10,'(2f18.9)') x,p(i)/norm
end do

print *, 'done'

end program

subroutine step_count(file_number,steps)
integer :: file_number, steps, ios,i
steps=0
do
 read(file_number,*,iostat=ios)
 if(ios.ne.0) exit
  steps=steps+1
end do
end subroutine

SUBROUTINE grid_min_max(num,grid_min,grid_max)
implicit none
integer :: num, ios
real*8 :: grid_min, grid_max
real*8 :: cv, dumm

rewind(num)
read(num,*,iostat=ios) dumm,cv
if (ios.ne.0) stop 'error reading colvar file'
grid_min=cv
grid_max=cv

rloop : DO
     read(num,*,iostat=ios)dumm,cv 
        if(ios.ne.0) exit rloop
     grid_min=MIN(cv,grid_min)
     grid_max=MAX(cv,grid_max)
end do rloop
end subroutine

