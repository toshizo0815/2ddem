module constant
    implicit none 

    real(8), parameter :: Pi = acos(-1.0d0)
    integer, parameter :: Nstep = 1000000, particle_number = 2, time_step = 1000      !time_stepは描写する時間間隔
    integer, parameter :: graph_step = 500, vtk_step = 1000
    integer, parameter :: step_fixed = 50000           !円筒を止めておくステップ数
    integer, parameter :: array_size = particle_number + 1    !粒子番号の列の大きさ


    !粒子の密度
    real(8), parameter :: rho_particle = 3.6d0    !アルミナ
    ! real(8), parameter :: rho_particle = 2.5d0    !ガラス

    real(8), parameter :: rho_cylinder = 1.2d0   !円筒の密度

    !real(8), parameter :: d_particle = 1.0d0
    real(8), parameter :: d_particle = 2.0d0
    !real(8), parameter :: d_particle = 6.0d0

    real(8), parameter :: d_cylinder = 0.06d0

    real(8), parameter :: g_original = 9.8d0, d_original = 0.001d0


    !粒子のパラメータ
    real(8), parameter :: time = 50.0d0, d = 1.0d0, m = 1.0d0, inertia = (2.0d0 / 5.0d0) * m * ((d / 2.0d0) ** 2)    !m,g,dで規格化
    real(8), parameter :: step_width = time / Nstep               !時間刻み幅
    real(8), parameter :: mu = 0.9d0, mu_cyliner = 0.9d0             !粒子, 円筒の密度
    real(8), parameter :: kn = 2.0d0 * 10 ** 5, kt = (2.0d0 / 7.0d0) * kn, e = 0.1d0
    real(8), parameter :: etan = -2 * log(e) * sqrt((m/2.0d0) * kn / ((acos(-1.0d0)) ** 2 + (log(e)) ** 2)), etat = etan
    !real(8), parameter :: etan = 0.0d0, etat = etan

    !円筒のパラメータ
    real(8), parameter :: r_cylinder = 30.0d0 / d_particle, deltad = 3.0d0 / d_particle, height = d     !deltadは円筒の厚さ, heightは円筒の長さ
    real(8), parameter :: m_cylinder = 6 * (rho_cylinder / rho_particle) * (2 * r_cylinder * deltad + deltad ** 2)
    real(8), parameter :: inertia_cylinder = (m_cylinder / 2.0d0) * ((r_cylinder + deltad) ** 2.0d0 + r_cylinder ** 2.0d0)

    !斜面のパラメータ
    real(8), parameter :: alpha_degree = 10.0d0
    real(8), parameter :: h = 0.0d0, alpha = alpha_degree * Pi / 180.0d0                         !hは初期位置の斜面の高さ, alphaは傾斜角度
    real(8), parameter :: a(0:2) = (/0.0d0, 1.0d0, -h/)

    !接触判定格子のパラメータ
    integer, parameter :: num_cell_x = 4 * r_cylinder, num_cell_y = 2 * r_cylinder, num_cell = num_cell_x * num_cell_y
    integer, parameter :: surplus = 20       !計算範囲外に余分にとる格子数(上下で2*surplus)
    real(8), parameter :: cell_size = 1.2d0 * d    !格子の大きさ
    real(8), parameter :: range_x = num_cell_x * cell_size, range_y = num_cell_y * cell_size      !計算領域のサイズ

    ! real(8), parameter :: g(0:1) = (/sin(alpha), -cos(alpha)/)
    real(8), parameter :: g(0:1) = (/0.0d0, 0.0d0/)

    real(8), parameter :: zero(0:1) = 0
    integer, parameter :: FillingRate = 125 * particle_number / (2 * r_cylinder) ** 2

    integer, parameter :: v_cell = 2 * r_cylinder

end module constant




!----------------------------------------------------------------------------------------------------------------





module print     !ファイル作成用モジュール
    use constant
    implicit none
contains
    subroutine makedirs(outdir)                                    !outdirというフォルダがなければ作る
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine makedirs

    subroutine print_profile_condition                             !python可視化用. 初期状態の情報をファイルにまとめる
        character(len = *), parameter :: base = "./initial_condition.d"
        character(len = 40) filename

        write(filename, '(a,i5.5,a)') base
        open(1,file = filename)    

        write(1,*) particle_number
        write(1,*) Nstep
        write(1,*) a(0)
        write(1,*) a(1)
        write(1,*) a(2)
        write(1,*) range_x
        write(1,*) range_y
        write(1,*) time_step

        close(1)
    end subroutine print_profile_condition

    subroutine print_profile_py(step, x, x_center, theta, theta_center)           !python可視化用. 粒子と円筒の位置と回転角度をファイルにまとめる
        real(8), intent(in) :: x(0:1,1:array_size), theta(1:array_size)
        real(8), intent(in) :: x_center(0:1), theta_center
        integer, intent(in) :: step
        integer :: i = 1
        character(len = *), parameter :: base = "./data_py/xy."
        character(len = 7) serial_num
        character(len = 40) filename

        if (mod(step, time_step) == 0) then
            write(serial_num, '(i7.7)') step
            write(filename, '(a,i7.7,a)') base//serial_num
            open(1,file = filename)

            if(particle_number > 0) then 
                do i = 1, particle_number
                    write(1,'(4(1x,f8.4))') x(0:1,i), theta(i), d / 2
                enddo
            endif

            write(1,'(4(1x,f8.4))') x_center(0:1), theta_center, r_cylinder
            write(1,'(4(1x,f8.4))') x_center(0:1), theta_center, r_cylinder + deltad


            if(range_x - r_cylinder <= x_center(0) .and. x_center(0) <= range_x) then 
                write(1,'(4(1x,f8.4))') x_center(0) - range_x, x_center(1), theta_center, r_cylinder
                write(1,'(4(1x,f8.4))') x_center(0) - range_x, x_center(1), theta_center, r_cylinder + deltad
            elseif(0 <= x_center(0) .and. x_center(0) <= r_cylinder) then 
                write(1,'(4(1x,f8.4))') x_center(0) + range_x, x_center(1), theta_center, r_cylinder
                write(1,'(4(1x,f8.4))') x_center(0) + range_x, x_center(1), theta_center, r_cylinder + deltad
            else 
                write(1,'(4(1x,f8.4))') x_center(0:1), theta_center, r_cylinder
                write(1,'(4(1x,f8.4))') x_center(0:1), theta_center, r_cylinder + deltad
            endif

            close(1)
        endif
    end subroutine print_profile_py


    subroutine mkvtk(step, dirname, x, theta)    !dirnameというフォルダ内にstep目のvtkファイルを作成
        integer, intent(in) :: step
        real(8), intent(in) :: x(0:1,1:array_size), theta(1:array_size)
        character(len=*), intent(in) :: dirname
        integer i
        character(len=40) str_step
        character(len=40) str_N

        write(str_step, '(i7.7)') step
        write(str_N, '(i4)') particle_number
        open(10,file = './'//trim(dirname)//'/data'//trim(str_step)//'.vtk')

        write(10,'(a)') '# vtk DataFile Version 5.1'
        write(10,'(a)') 'Data.vtk'
        write(10,'(a)') 'ASCII'
        write(10,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(10,'(a)') ' '

        !点の位置情報
        write(10,'(a)') 'POINTS '//trim(str_N)//' float'
        do i = 1, particle_number
            write(10,'(f8.2,f8.2,f8.2)') x(:,i), 0
        enddo
        write(10,'(a)') ' '

        !半径の情報
        write(10,'(a)') 'POINT_DATA '//trim(str_N)
        write(10,'(a)') 'SCALARS radius float'
        write(10,'(a)') 'LOOKUP_TABLE default'
        do i = 1, particle_number
            write(10,*) d / 2.0d0
        enddo
        write(10,'(a)') ' '

        !方向指定ベクトルの向き
        write(10,'(a)') 'VECTORS direction float'
        do i = 1, particle_number
            write(10,'(f6.2,f6.2,f6.2)') (d/2.0d0)*cos(theta(i)), (d/2.0d0)*sin(theta(i)), 0
        enddo

        close(10)
    end subroutine mkvtk

    subroutine mkvtk_cylinder(dirname, x_center)
        real(8), intent(in) :: x_center(0:1,0:Nstep)
        character(len=*), intent(in) :: dirname
        integer step
        character(len=40) str_step
        character(len=40) str_N

        do step = 0, Nstep
            if (mod(step,vtk_step)==0) then
                write(str_step, '(i7.7)') step
                write(str_N, '(i4)') 2
                open(10,file = './'//trim(dirname)//'/data_cylinder'//trim(str_step)//'.vtk')
                write(10,'(a)') '# vtk DataFile Version 5.1'
                write(10,'(a)') 'Data.vtk'
                write(10,'(a)') 'ASCII'
                write(10,'(a)') 'DATASET UNSTRUCTURED_GRID'
                write(10,'(a)') ' '
                write(10,'(a)') 'POINTS '//trim(str_N)//' float'
                write(10,*) x_center(:,step), 0
                write(10,*) x_center(:,step), 0
                write(10,'(a)') ' '
                write(10,'(a)') 'CELL_TYPES '//trim(str_N)
                write(10,*) 1
                write(10,*) 1
                write(10,'(a)') ' '
                write(10,'(a)') 'POINT_DATA '//trim(str_N)
                write(10,'(a)') 'SCALARS radius float'
                write(10,'(a)') 'LOOKUP_TABLE default'
                write(10,*) r_cylinder
                write(10,*) r_cylinder + deltad
                close(10)
            endif
        enddo
    end subroutine mkvtk_cylinder
end module print





!--------------------------------------------------------------------------------------------------------------------





module calculate_particle_collision     !数値計算のサブルーチンをまとめたモジュール
    use constant
    implicit none
contains

    subroutine initial_value(x, x_center)     !粒子の初期位置を設定するサブルーチン
        real(8), intent(in) :: x_center(0:1)
        real(8), intent(inout) :: x(0:1,1:array_size)
        integer :: i = 0
        real(8) :: rnd1, rnd2
        integer(4), allocatable :: seed(:)
        integer(4) :: c, sz

        call random_seed( size = sz ) !コンパイラによってシードの配列サイズが異なるので，まずはシードの配列サイズを取得する．
        allocate( seed(sz) ) !得られたサイズでseedをallocateする．
        call random_seed( get = seed ) !デフォルトでセットされているシードを取得する．
        call system_clock( count = c ) !システム時刻を取得して，
        seed(1) = c                    !シードにその時刻を使用することで，毎回異なる乱数となるようにする．
        call random_seed( put = seed ) !組み直したシードを設定する．
        
        if(particle_number > 0) then 
            do i = 1, particle_number

                call random_number(rnd1)
                call random_number(rnd2)

                x(0,i) = (r_cylinder - (d/2.0d0)) * sqrt(rnd1) * cos(2*Pi*rnd2) + x_center(0)
                x(1,i) = (r_cylinder - (d/2.0d0)) * sqrt(rnd1) * sin(2*Pi*rnd2) + x_center(1)
            enddo
        endif 

    end subroutine initial_value

    function vec_length(p, q) result(length)         !2点間の距離を計算する関数
        real(8), intent(in) :: p(0:1), q(0:1)
        real(8) length

        length = sqrt((p(0) - q(0))**2 + (p(1) - q(1))**2)
    end function vec_length

    function point_line_distance(x, p) result(length)    !点と直線の距離を計算する関数
        real(8), intent(in) :: x(0:1)
        real(8), intent(in) :: p(0:2)
        real(8) length
            
        length = abs(p(0) * x(0) + p(1) * x(1) + p(2)) / sqrt((p(0)) ** 2 + (p(1)) ** 2)
    end function point_line_distance

    function inner_product(p, q) result(product)       !内積を計算するサブルーチン
        real(8), intent(in) :: p(:), q(:)
        real(8) product
        integer i 

        product = 0
        do i = 1, size(p)
            product = product + p(i) * q(i)
        enddo
    end function inner_product

    subroutine periodic_boundary(x, round)        
        real(8), intent(inout) :: x(0:1)
        integer, intent(out) :: round

        if(x(0) >= range_x) then 
            x(0) = x(0) - range_x
            round = round + 1
        elseif(x(0) < 0) then 
            x(0) = x(0) + range_x
        endif 
    end subroutine periodic_boundary

    subroutine cal_cog(x, x_g, round)
        real(8), intent(in) :: x(0:1,1:array_size)
        integer, intent(in) :: round(1:array_size)
        real(8), intent(out) :: x_g(0:1)
        real(8) :: unit(0:1) = (/1,0/)
        real(8)  numerator(0:1)
        integer :: i = 0

        numerator(0:1) = 0
        if(particle_number > 0) then 
            do i = 1, particle_number
                numerator(0:1) = numerator(0:1) + x(0:1,i) + round(i) * range_x * unit(0:1)
            enddo
        
            x_g(0:1) = numerator(0:1) / particle_number
        endif 
        numerator(0:1) = 0
    end subroutine cal_cog

    subroutine rot_x(x, beta, x_rot)
        real(8), intent(in) :: x(0:1)
        real(8), intent(in) :: beta
        real(8), intent(out) :: x_rot(0:1)

        x_rot(0) = cos(beta) * x(0) + sin(beta) * x(1)
        x_rot(1) = cos(beta) * x(1) - sin(beta) * x(0)
    end subroutine rot_x
end module calculate_particle_collision






!------------------------------------------------------------------------------------------------------------






module calculate_contact
    use calculate_particle_collision
    use constant
    implicit none
contains

    subroutine particle_contact(d1, d2, x1, x2, v1, v2, w1, w2, deltat_prev, deltat, fn, ft, torque, Us, slip)
        real(8), intent(in) :: d1, d2
        real(8), intent(in) :: x1(0:1), x2(0:1), v1(0:1), v2(0:1)
        real(8), intent(in) :: w1, w2
        real(8), intent(inout) :: deltat_prev(0:1), deltat(0:1)
        integer, intent(out) :: slip
        real(8), intent(out) :: fn(0:1), ft(0:1)
        real(8), intent(out) :: torque
        real(8), intent(inout) :: Us

        real(8) :: u = 0
        real(8) :: n(0:1) = 0, t(0:1) = 0
        real(8) :: vn(0:1) = 0, vt(0:1) = 0, vt_rot(0:1) = 0
        real(8) :: deltan(0:1) = 0, deltat_tmp(0:1) = 0

        slip = 0

        !接触判定
        if (vec_length(x1(0:1), x2(0:1)) <= (d1 + d2)/2) then
            u = 1
        else
            u = 0
        endif

        if (x1(0) /= x2(0) .or. x1(0) /= x2(0)) then 
            n(0:1) = (x1(0:1) - x2(0:1)) / vec_length(x1(0:1), x2(0:1))
        endif 

        t(0) = n(1)
        t(1) = -n(0)

        deltan(0:1) = (vec_length(x1(0:1), x2(0:1)) - ((d1 + d2)/2)) * n(0:1)
        vn(0:1) = inner_product(v1(0:1)-v2(0:1), n(0:1)) * n(0:1)
        vt(0:1) = inner_product(v1(0:1)-v2(0:1), t(0:1)) * t(0:1)
        vt_rot(0:1) = (inner_product(v1(0:1)-v2(0:1), t(0:1)) + (d1 * w1 + d2 * w2)/2) * t(0:1)
        deltat(0:1) = u * (inner_product(deltat_prev(0:1), t(0:1)) * t(0:1) + vt_rot(0:1) * step_width)
        deltat_tmp(0:1) = deltat_prev(0:1)
        deltat_prev(0:1) = deltat(0:1)

        Us = u * (1.0d0 / 2.0d0) * (kn * vec_length(deltan(0:1),zero(0:1))**2 + kt * vec_length(deltat(0:1),zero(0:1))**2)

        fn(0:1) = u * (-kn * deltan(0:1) -etan * vn(0:1))
        ft(0:1) = u * (-kt * deltat(0:1) -etat * vt_rot(0:1))

        if (vec_length(ft(0:1),zero(0:1)) > mu * vec_length(fn(0:1), zero(0:1))) then 
            slip = 1
            ft(0:1) = mu * vec_length(fn(0:1), zero(0:1)) * ft(0:1) / vec_length(ft(0:1), zero(0:1))
            deltat_prev(0:1) = inner_product(deltat_tmp(0:1), t(0:1)) * t(0:1)
        endif 

        torque =  (d1 / 2) * inner_product(ft(0:1), t(0:1))
    end subroutine particle_contact

    subroutine particle_contact_dissipation(d1, d2, x1, x2, v1, v2, w1, w2, x_center, dis, slip, ft)
        integer, intent(in) :: slip
        real(8), intent(in) :: d1, d2
        real(8), intent(in) :: x1(0:1), x2(0:1), v1(0:1), v2(0:1)
        real(8), intent(in) :: w1, w2
        real(8), intent(in) :: x_center(0:1)
        real(8), intent(in) :: ft(0:1)
        real(8), intent(inout) :: dis(0:v_cell-1,0:v_cell-1)

        integer :: nx=0, ny=0
        real(8) :: dissipation = 0.0d0
        real(8) :: u = 0
        real(8) :: n(0:1) = 0, t(0:1) = 0
        real(8) :: vn(0:1) = 0, vt_rot(0:1) = 0

        !接触判定
        if (vec_length(x1(0:1), x2(0:1)) <= (d1 + d2)/2) then
            u = 1
        else
            u = 0
        endif

        if (x1(0) /= x2(0) .or. x1(1) /= x2(1)) then 
            n(0:1) = (x1(0:1) - x2(0:1)) / vec_length(x1(0:1), x2(0:1))
        endif 

        t(0) = n(1)
        t(1) = -n(0)

        vn(0:1) = inner_product(v1(0:1)-v2(0:1), n(0:1)) * n(0:1)     
        vt_rot(0:1) = (inner_product(v1(0:1)-v2(0:1), t(0:1)) + (d1 * w1 + d2 * w2)/2) * t(0:1)
        nx=floor((x1(0)+x2(0))/2.0d0 - x_center(0))+v_cell/2
        ny=floor((x1(1)+x2(1))/2.0d0 - x_center(1))+v_cell/2
        dis(nx,ny) = dis(nx,ny) +(- u * etan * (inner_product(vn(:),vn(:))) &
         + slip*u*inner_product(ft(:),vt_rot(:)))*step_width

        !dis(nx,ny) = dis(nx,ny) +(- u * etan * (inner_product(vn(:),vn(:))))*step_width

        !  dis(nx,ny) = dis(nx,ny) + slip*u*inner_product(ft(:),vt_rot(:))
    end subroutine particle_contact_dissipation



    subroutine particle_cylinder_contact(x_center, v_center, w_wall, x, v, w, deltat_prev, deltat, fwn, fwt, &
    torque1, torque2, u_wall, Us, slip_wall)
        integer, intent(out) :: slip_wall 
        real(8), intent(in) :: w_wall, w
        real(8), intent(in) :: x_center(0:1), v_center(0:1), x(0:1), v(0:1), deltat_prev(0:1)
        integer, intent(out) :: u_wall
        real(8), intent(out) :: torque1, torque2
        real(8), intent(out) :: fwn(0:1), fwt(0:1), deltat(0:1)
        real(8), intent(inout) :: Us
        real(8) :: vn(0:1) = 0, vt(0:1) = 0, n(0:1) = 0, t(0:1) = 0, deltan(0:1) = 0, vt_rot(0:1) = 0

        slip_wall = 0

        if(vec_length(x_center(0:1), x(0:1)) + d / 2 > r_cylinder) then
            u_wall = 1
        else
            u_wall = 0
        endif 
     
        if (vec_length(x_center(0:1), x(0:1)) == 0) then
            n(0:1) = 0
        else 
            n(0:1) = (x_center(0:1) - x(0:1)) / vec_length(x_center(0:1), x(0:1))
        endif 

        t(0) = n(1)
        t(1) = -n(0)

        vn(0:1) = inner_product(v(0:1) - v_center(0:1), n(0:1)) * n(0:1)
        vt(0:1) = inner_product(v(0:1) - v_center(0:1), t(0:1)) * t(0:1)
        deltan(0:1) = u_wall * (r_cylinder - (vec_length(x_center(0:1), x(0:1)) + d / 2)) * n(0:1)
        vt_rot(0:1) = vt(0:1) + (d / 2 * w - r_cylinder * w_wall) * t(0:1)
        deltat(0:1) = (inner_product(deltat_prev(0:1), t(0:1)) * t(0:1) + vt_rot(0:1) * step_width) * u_wall

        Us = (1.0d0 / 2.0d0) * (kn * vec_length(deltan(0:1),zero(0:1))**2 + kt * vec_length(deltat(0:1),zero(0:1))**2)

        fwn(0:1) = -u_wall * (kn * deltan(0:1) + etan * vn(0:1))
        fwt(0:1) = -u_wall * (kt * deltat(0:1) + etat * vt_rot(0:1))

        if (vec_length(fwt(0:1), zero(0:1)) > mu_cyliner * vec_length(fwn(0:1), zero(0:1))) then 
            slip_wall = 1
            fwt(0:1) = mu_cyliner * vec_length(fwn(0:1), zero(0:1)) * t(0:1) * sign(1.0d0, inner_product(fwt(0:1), t(0:1)))

            deltat(0:1) = deltat(0:1) - vt_rot(0:1) * step_width
        endif 

        torque1 = inner_product(fwt(0:1), t(0:1)) * (d / 2)
        torque2 = inner_product(-fwt(0:1), t(0:1)) * r_cylinder

    end subroutine particle_cylinder_contact

    subroutine particle_cylinder_contact_dissipation(x_center, v_center, w_wall, x, v, w, u_wall, Us, dis, slip_wall, fwt)
            integer, intent(in) :: slip_wall
            real(8), intent(in) :: w_wall, w
            real(8), intent(in) :: x_center(0:1), v_center(0:1), x(0:1), v(0:1)
            real(8), intent(in) :: fwt(0:1)
            integer, intent(out) :: u_wall
            real(8), intent(inout) :: Us
            real(8), intent(inout) :: dis(0:v_cell-1,0:v_cell-1)
            real(8) :: vn(0:1) = 0, vt(0:1) = 0, n(0:1) = 0, t(0:1) = 0, vt_rot(0:1) = 0
            integer :: nx = 0, ny = 0
    
            if(vec_length(x_center(0:1), x(0:1)) + d / 2 > r_cylinder) then
                u_wall = 1
            else
                u_wall = 0
            endif 
         
            if (vec_length(x_center(0:1), x(0:1)) == 0) then
                n(0:1) = 0
            else 
                n(0:1) = (x_center(0:1) - x(0:1)) / vec_length(x_center(0:1), x(0:1))
            endif 
    
            t(0) = n(1)
            t(1) = -n(0)
    
            vn(0:1) = inner_product(v(0:1) - v_center(0:1), n(0:1)) * n(0:1)
            vt(0:1) = inner_product(v(0:1) - v_center(0:1), t(0:1)) * t(0:1)
            vt_rot(0:1) = vt(0:1) + (d / 2 * w - r_cylinder * w_wall) * t(0:1)

            nx=floor(x(0) - x_center(0))+v_cell/2
            ny=floor(x(1) - x_center(1))+v_cell/2
            dis(nx,ny) = dis(nx,ny) - u_wall * etan * (inner_product(vn(:),vn(:))) !+ slip_wall*u_wall*inner_product(fwt(:),vt_rot(:))
    
        end subroutine particle_cylinder_contact_dissipation

    subroutine cylinder_slope_contact(x, v, w, deltat, fn, ft, torque)
        real(8), intent(in) ::  w
        real(8), intent(in) :: x(0:1), v(0:1)
        real(8), intent(out) :: fn(0:1), ft(0:1)
        real(8), intent(out) :: torque
        real(8), intent(inout) :: deltat(0:1)
        integer :: u_slope = 0
        real(8) :: L = 0, r_outside = r_cylinder + deltad
        real(8) :: n(0:1) = 0, t(0:1) = 0, vn(0:1) = 0, vt(0:1) = 0, v_contact(0:1) = 0, deltan(0:1) = 0


        L = point_line_distance(x(0:1), a(0:2))

        if(L <= r_outside) then
            u_slope = 1
        else 
            u_slope = 0
        endif 

        n(0) = a(0) / sqrt((a(0)) ** 2 + (a(1)) ** 2)
        n(1) = a(1) / sqrt((a(0)) ** 2 + (a(1)) ** 2)
        t(0) = n(1)
        t(1) = -n(0)

        vn(0:1) = inner_product(v(0:1), n(0:1)) * n(0:1)
        vt(0:1) = v(0:1) - vn(0:1)
        v_contact(0:1) = vt(0:1) + r_outside * w * t(0:1)

        deltan(0:1) = (L - r_outside) * n(0:1)
        deltat(0:1) = u_slope * (deltat(0:1) + v_contact(0:1) * step_width)
        
        fn(0:1) = u_slope * (-kn * deltan(0:1) - etan * vn(0:1))
        ft(0:1) = u_slope * (-kt * deltat(0:1) - etat * v_contact(0:1))

        ! fn(0:1) = u_slope * (-kn * deltan(0:1))    !斜面と円筒の接触による減衰なし
        ! ft(0:1) = u_slope * (-kt * deltat(0:1))

        torque = inner_product(ft(0:1), t(0:1)) * r_outside


    end subroutine cylinder_slope_contact

end module calculate_contact

!--------------------------------------------------------------------------------------------------------------------------






module neighbor_paritcle_search
    use constant
    implicit none
contains

    subroutine get_cell_index(x, idx_x, idx_y, cell_index) 
        real(8), intent(in) :: x(0:1)
        integer, intent(out) :: idx_x, idx_y, cell_index

        idx_x = floor(x(0) / cell_size) + 1
        idx_y = floor(x(1) / cell_size) + 1 + surplus
        cell_index = num_cell_x * (idx_y - 1) + idx_x          !num_cellは偶数限定, セル番号は0から
    end subroutine get_cell_index

end module neighbor_paritcle_search







!-------------------------------------------------------------------------------------------------------------------





module theory
    use constant
    use calculate_particle_collision

contains
    subroutine theory100(step, x_theory, inertia_sub)   !inertia_subは粉体の慣性モーメント
        integer, intent(in) :: step
        real(8), intent(in) :: inertia_sub
        real(8), intent(inout) :: x_theory
        real(8) :: nume = 0, deno = 0, inertia_all = 0

        inertia_all = inertia_cylinder + inertia_sub

        nume = (m_cylinder + (2880.0d0/(d_particle**2)) * m) * g(0)
        deno = (m_cylinder + (2880.0d0/(d_particle**2)) * m) + inertia_all / (r_cylinder + deltad)**2

        x_theory = (1.0d0/2.0d0) * (nume / deno) * (step * step_width)**2
    end subroutine

    subroutine theory0(step, x_theory)
        integer, intent(in) :: step
        real(8), intent(inout) :: x_theory
        real(8) :: nume = 0, deno = 0

        nume = m_cylinder * g(0)
        deno = m_cylinder + inertia_cylinder / (r_cylinder + deltad)**2

        x_theory = (1.0d0/2.0d0) * (nume / deno) * (step * step_width)**2
    end subroutine
end module theory


!-------------------------------------------------------------------------------------------------------------------




program dem_2d
    use calculate_particle_collision
    use print
    use constant
    use calculate_contact
    use neighbor_paritcle_search
    use theory
    implicit none

    real(8) :: x(0:1,1:array_size) = 0, v(0:1,1:array_size) = 0, p(0:1,1:array_size) = 0
    real(8) :: x_prev(0:1,1:array_size) = 0, v_prev(0:1,1:array_size) = 0
    real(8) :: x_tmp(0:1) = 0, x_g(0:1) = 0, x_g_rot(0:1) = 0
    real(8) :: deltan(0:1) = 0, vt_rot(0:1) = 0, vt_rot_wall(0:1) = 0
    real(8) :: deltat_prev(0:1,1:array_size,1:array_size) = 0, deltat(0:1,1:array_size,1:array_size) = 0
    real(8) :: vn(0:1) = 0, vt(0:1) = 0
    real(8) :: torque = 0
    real(8) :: torque_all(1:array_size) = 0, torque_reaction(1:array_size) = 0
    real(8) :: torque_wall_all = 0
    real(8) :: w(1:array_size) = 0, theta(1:array_size) = 0
    real(8) :: w_prev(1:array_size) = 0, theta_prev(1:array_size) = 0
    real(8) :: n(0:1) = 0, t(0:1) = 0
    real(8) :: fcn(0:1) = 0, fct(0:1) = 0, fc_reaction(0:1,1:array_size) = 0
    real(8) :: fwn(0:1,1:array_size) = 0, fwt(0:1,1:array_size) = 0, fw_all(0:1) = 0
    real(8) :: fc_all(0:1,1:array_size) = 0
    real(8) :: kinetic_energy(1:array_size) = 0, potential_energy(1:array_size) = 0
    real(8) :: energy_sum = 0
    real(8) :: time1 = 0, time2 = 0
    integer :: i = 0, j = 0, step = 0, u = 0, round = 0, progress = 0
    integer :: round_particle(1:array_size) = 0


    integer :: u_wall = 0, u_slope = 0
    real(8) :: torque_wall(1:array_size) = 0, torque_wall2(1:array_size) = 0
    real(8) :: x_center(0:1,0:Nstep) = 0, v_center(0:1,0:Nstep), p_center(0:1,0:Nstep) = 0
    real(8) :: x_center_tmp(0:1) = 0, x_contact(0:1) = 0, x_contact_rot(0:1) = 0
    real(8) :: x_original(0:1,0:Nstep) = 0, x_original_rot(0:1) = 0
    real(8) :: x_theory = 0
    real(8) :: theta_center(0:Nstep) = 0, w_wall(0:Nstep) = 0
    real(8) :: deltat_wall_prev(0:1) = 0, deltat_wall(0:1) = 0
    real(8) :: deltat_cylinder(0:1) = 0
    real(8) :: fn_cylinder(0:1,0:Nstep), ft_cylinder(0:1,0:Nstep)
    real(8) :: torque_cylinder(0:Nstep) = 0
    real(8) :: L = 0

    real(8) :: inertia_theory = 0, inertia_simulation = 0

    real(8) :: Uall = 0, Us = 0

    real(8) :: dis(0:v_cell-1,0:v_cell-1) = 0.0d0, dis_ave(0:v_cell-1,0:v_cell-1) = 0.0d0
    integer :: counter = 0, slip = 0, slip_wall = 0


    integer :: i_cell = 0, j_cell = 0
    integer :: particle_index = 0, idx_x = 0, idx_y = 0, cell_index = 0, last_preview = 0, k = 0
    integer :: c_idx_x = 0, c_idx_y = 0, tmp = 0, c_idx = 0
    integer :: nextOf(1:particle_number) = -1
    integer :: first(1:num_cell_x * (num_cell_y + 2 * surplus)) = -1
    integer :: last(1:num_cell_x * (num_cell_y + 2 * surplus)) = -1

    character(len=40) filename_xt, filename_vt, filename_fr
    character(len = *), parameter :: base_xt = "./xt_", base_vt = "./vt_data/vt_n=", base_fr = './fr_'
    character(len = *), parameter :: base_data = './data_'
    character(len = 3) particle_num_char
    character(len = 2) alpha_char 
    character(len = 1) d_char 
    character(len = 7) step_char
    character(len = 3) :: chmyrank = "ave"

    write(particle_num_char, '(i3.3)') FillingRate
    write(alpha_char, '(i2.2)') int(alpha_degree) 
    write(d_char, '(i1.1)') int(d_particle)
    write(filename_xt, '(a)') base_xt//d_char//'_'//alpha_char//'_'//particle_num_char
    write(filename_vt, '(a,i6.6,a)') base_vt//particle_num_char
    write(filename_fr, '(a)') base_fr//d_char//'_'//alpha_char//'_'//particle_num_char



    open(30,file='dissipation_t.d')

    call makedirs('./data_py')
    call makedirs('./vtk_data')
    call print_profile_condition

!初期状態を設定 {
    do step = 0, step_fixed
        x_center(0:1,step) = (/r_cylinder + deltad, r_cylinder + deltad/)
    enddo

    x_prev(0,1) = 10.0d0
    x_prev(1,1) = 5.0d0
    x_prev(0,2) = 0.0d0
    x_prev(1,2) = 5.0d0
    v_prev(0,1) = -1.0d0
    w_prev(1) = -5.0d0
    v_prev(0,2) = 0.0d0

    ! call initial_value(x_prev(0:1,1:array_size), x_center(0:1,0))
    ! call mkvtk(0, './vtk_data', x_prev(:,:), theta_prev(:))

    do step = 1, Nstep 
        do i = 1, particle_number
            Uall = Uall + 0.5d0 * m * (vec_length(v(:,i),zero(:))**2) + 0.5d0 * inertia * (w(i))**2
        enddo 

        if(particle_number > 0) then 

            !各粒子の格子住所を配列nextOfに登録 {
            do c_idx = 1, num_cell_x * (num_cell_y + 2 * surplus)
                first(c_idx) = -1
                last(c_idx) = -1
            enddo

            do particle_index = 1, particle_number
                nextOf(particle_index) = -1
            enddo

            do particle_index = 1, particle_number
                call get_cell_index(x_prev(0:1,particle_index), idx_x, idx_y, cell_index)

                last_preview = last(cell_index)
                
                last(cell_index) = particle_index
        
                if (last_preview == -1) then 
                    first(cell_index) = particle_index
                else 
                    nextOf(last_preview) = particle_index
                endif 

                cell_index = 0
            enddo 
            !}


            do i = 1, particle_number   !注目する粒子
                call get_cell_index(x_prev(0:1,i), idx_x, idx_y, cell_index)

                do i_cell = -1, 1
                    do j_cell = -num_cell_x, num_cell_x, num_cell_x

                        c_idx_y = j_cell
                        c_idx_x = i_cell
                        if(i_cell + idx_x > num_cell_x) then      !(判定粒子の格子番号) + (注目粒子のx方向の格子番号)
                            tmp = first(cell_index + c_idx_x + c_idx_y - num_cell_x)
                        elseif(i_cell + idx_x < 1) then 
                            tmp = first(cell_index + c_idx_x + c_idx_y + num_cell_x)
                        else 
                            tmp = first(cell_index + c_idx_x + c_idx_y)
                        endif 


                        do 
                            if (tmp == -1) then 
                                exit
                            elseif (tmp /= i) then
                                if(i < tmp) then
                                    if(range_x - x_prev(0,i) + x_prev(0,tmp) < d) then 
                                        x_tmp(0) = range_x
                                    elseif(range_x - x_prev(0,tmp) + x_prev(0,i) < d) then 
                                        x_tmp(0) = -range_x 
                                    else 
                                        x_tmp(0) = 0
                                    endif 

                                    call particle_contact(d, d, x_prev(0:1,i), x_prev(0:1,tmp) + x_tmp(0:1), v_prev(0:1,i), &
                                    v_prev(0:1,tmp), w_prev(i), w_prev(tmp), deltat_prev(0:1,i,tmp), deltat(0:1,i,tmp), &
                                    fcn(0:1), fct(0:1), torque, Us, slip)
                                    fc_reaction(0:1,tmp) = fc_reaction(0:1,tmp) - fcn(0:1) - fct(0:1)
                                    torque_reaction(tmp) = torque_reaction(tmp) + torque
                                    Uall = Uall + Us
                                    Us = 0.0d0

                                    call particle_contact_dissipation(d, d, x_prev(0:1,i), x_prev(0:1,tmp) + x_tmp(0:1), v_prev(0:1,i), &
                                    v_prev(0:1,tmp), w_prev(i), w_prev(tmp), x_center(0:1,step-1), dis(:,:),slip, fct(0:1))
                                endif

                                fc_all(0:1,i) = fc_all(0:1,i) + fcn(0:1) + fct(0:1)
                                torque_all(i) = torque_all(i) + torque
                                fcn(0:1) = 0
                                fct(0:1) = 0
                                torque = 0
                            endif
                            tmp = nextOf(tmp)
                        enddo
                        tmp = -1
                    enddo
                enddo


                fc_all(0:1,i) = fc_all(0:1,i) + fc_reaction(0:1,i)
                torque_all(i) = torque_all(i) + torque_reaction(i)

                torque_all(i) = torque_all(i) + torque_wall(i)    !粒子iに作用するトルク
                fc_all(0:1,i) = fc_all(0:1,i) + fwn(0:1,i) + fwt(0:1,i)    
                fw_all(0:1) = fw_all(0:1) - fwn(0:1,i) - fwt(0:1,i)

                torque_wall(i) = 0
                torque_wall2(i) = 0
                fwn(0:1,i) = 0
                fwt(0:1,i) = 0

                v(0:1,i) = v_prev(0:1,i) + (fc_all(0:1,i) / m + g(0:1)) * step_width
                x(0:1,i) = x_prev(0:1,i) + v(0:1,i)* step_width
                w(i) = w_prev(i) + torque_all(i) * step_width / inertia
                theta(i) = theta_prev(i) + w(i) * step_width

                v_prev(0:1,i) = v(0:1,i)
                x_prev(0:1,i) = x(0:1,i)
                w_prev(i) = w(i)
                theta_prev(i) = theta(i)
            enddo

            ! ステップ数vtk_stepごとにvtkファイルを書き出す
            if (mod(step, vtk_step)==0) then
                call mkvtk(step, './vtk_data', x(:,:), theta(:))
            endif
        endif 

        ! if (mod(step,10000)==0 .and. step >= step_fixed) then
        !     write(step_char,'(i7.7)') step - step_fixed
        !     open(30,file='dissipation_data_'//step_char)

        !     do i = 0, v_cell - 1
        !         do j = 0, v_cell - 1
        !             write(30, *) j, i, -dis(j,i)
        !         enddo
        !     enddo
        !     close(30)
        ! endif
        dis_ave(:,:) = dis_ave(:,:) + dis(:,:)
        write(30,*) step, sum(-dis(:,:)), Uall, sum(-dis_ave(:,:)), Uall+sum(-dis_ave(:,:))
        Uall = 0.0d0
        dis(:,:) = 0.0d0


        fc_reaction(0:1,1:array_size) = 0
        fw_all(0:1) = 0
        fc_all(0:1,1:array_size) = 0
        torque_reaction(1:array_size) = 0
        torque_all(1:array_size) = 0
        torque_wall_all = 0
        deltat(0:1,1:array_size,1:array_size) = 0
        deltat_prev(0:1,1:array_size,1:array_size) = 0

        ! !円筒の変位の時間変化
        ! if (step >= step_fixed .and. mod(step, graph_step) == 0) then
        !     write(15,*) (step - step_fixed)*step_width*sqrt(dble(d_particle)/(g_original*1000.0d0)), &
        !     (x_original(0,step) -x_center(0,0))/(2*Pi*(r_cylinder+deltad))
        ! endif

        if(step == (Nstep / 10) * (progress + 1)) then 
            progress = progress + 1
            print *, '進捗',progress * 10,'%'
        endif 
    enddo
    !stepのdoループが終了

    write(*,*) counter
    write(*,*) sum(-dis_ave(:,:))
    write(*,*) (m_cylinder+particle_number*m)*sin(alpha)*v_center(0,400000)

    call mkvtk_cylinder('./vtk_data', x_center(:,:))

    close(15)
    close(20)
    close(30)


end program dem_2d