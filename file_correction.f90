program main 
    implicit none

    integer, parameter :: end_file = 3901
    integer i
    real(8) :: t(1:end_file), x(1:end_file)

    open(10,file='xt_6_10_032')
    open(20,file='xt_6_10_032_v2')
    do i = 1,end_file
        read(10, *) t(i), x(i)
    enddo
    x(2668:end_file) = x(2668:end_file) - 0.69

    do i = 1,end_file
        write(20, *) t(i), x(i)
    enddo
end program main