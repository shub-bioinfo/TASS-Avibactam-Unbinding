program histogram
implicit none
real*8 :: grid_min1, grid_max1, grid_width1, &
          grid_min2, grid_max2, grid_width2, &
           s1, s2, x, dum, norm, cv1, cv2, prob        
ALLOCATABLE cv1(:), cv2(:), prob(:,:)
integer :: nbin1, nbin2, bin, index1, index2, i, i_s1, i_s2, steps, step_loop

!To perform the convergence 
!read (*,*) step_loop
open(2,file='COLVAR',status='old')
!open(2,file='colvar',status='old')

! 1. Provide the grid min, grid max and the grid width for both the CVs.
! 2. If we have to give the number of MD steps, mention after the steps count 
!     suroutine, steps=number_of_steps.
! 3. Check for the Column which has the correspoding CV value in line 78.
! 4. Probability distribution generated in "output" file.

grid_width1 = 0.170d0 
grid_min1=13.0d0
grid_max1=47.00d0

!grid_width2 =0.005
!grid_min2=0.0d0
!grid_max2=0.5d0

grid_width2 =1.0
grid_min2=25.0d0
grid_max2=125.0d0
!CALL grid_min_max(2,grid_min,grid_max)

print*, 'CV1', grid_min1, grid_max1
print*, 'CV2', grid_min2, grid_max2

nbin1 = NINT((grid_max1-grid_min1)/grid_width1)+1
nbin2 = NINT((grid_max2-grid_min2)/grid_width2)+1

print *, 'Number of bins for CV1 =', nbin1
print *, 'Number of bins for CV2 =', nbin2

!rewind(2)

call step_count(2,steps)
print *, 'steps=', steps

rewind(2)

allocate(cv1(steps),cv2(steps))
allocate(prob(nbin1,nbin2))

!To give same number of MD steps for the simulation run.
!steps = 2500000
print *, 'steps Used=', steps


!Calculating Probability 
norm=0.d0
prob=0.d0
do i=1,steps
!  write(*,*)i
  read(2,*) dum,cv1(i),dum,dum,cv2(i)
  index1 = nint((cv1(i)-grid_min1)/grid_width1)+1 
  index2 = nint((cv2(i)-grid_min2)/grid_width2)+1 
!   print*, i, steps
    if(index1.gt.0 .and. index2.gt.0 .and. index1.le.nbin1 .and. index2.le.nbin2) then
      prob(index1,index2) = prob(index1,index2)+ 1.d0
      norm=norm+1.d0
    end if
end do

print*, 'norm =', norm
norm = norm*grid_width1*grid_width2
print*, 'norm*area =', norm

!Writting Probability
open(10,file='output',status='replace')
  do i_s1=1,nbin1
    s1= grid_min1 + float(i_s1-1)*grid_width1
    do i_s2=1,nbin2
      s2= grid_min2 + float(i_s2-1)*grid_width2
      write(10,'(3F18.9)') s1, s2, prob(i_s1,i_s2)/norm
    end do 
  write(10,*) 
 end do

print *, 'Unbiased distribution written in output file'

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

