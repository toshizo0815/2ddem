module constant
    implicit none 

    real(8), parameter :: Pi = acos(-1.0d0)

    !粒子の密度
    real(8), parameter :: rho_particle = 3.6d0    !アルミナ
    ! real(8), parameter :: rho_particle = 2.5d0    !ガラス

    real(8), parameter :: rho_cylinder = 1.2d0   !円筒の密度

    real(8), parameter :: d_particle = 1.0d0
    !real(8), parameter :: d_particle = 2.0d0
    !real(8), parameter :: d_particle = 6.0d0

    real(8), parameter :: g_original = 9.8d0


    !粒子のパラメータ
    real(8), parameter :: d = 1.0d0, m = 1.0d0    !m,g,dで規格化

    !円筒のパラメータ
    real(8), parameter :: R = 30.0d0 / d_particle, deltad = 3.0d0 / d_particle, height = d     !deltadは円筒の厚さ, heightは円筒の長さ
    real(8), parameter :: Mc = 6.0d0 * (rho_cylinder / rho_particle) * (2.0d0 * R * deltad + deltad ** 2.0d0)

    !斜面のパラメータ
    real(8), parameter :: alpha_degree = 2.0d0
    real(8), parameter :: h = 0.0d0, alpha = alpha_degree * Pi / 180.0d0                         !hは初期位置の斜面の高さ, alphaは傾斜角度

    ! real(8), parameter :: g(0:1) = (/sin(alpha), -cos(alpha)/)
    real(8), parameter :: g = 1
end module constant


program cylinder 
    use constant
    implicit none 

    integer :: N = 100, i = 0, step = 0
    real(8) :: Mp = 0.0d0
    real(8) :: R_all = R + deltad, l_all = 0
    real(8) :: Ip = 0, Iall = 0, lg = 0
    real(8) :: A = 0, B = 0, C = 0
    real(8) :: judge = 0

    real(8) :: alpha_eq1 = 0, alpha_eq2 = 0, phi_eq = 0, Sxy1 = 0.50d0, Sxy2 = 0
    integer :: xi = 0, Nxy = 500
    real(8) :: theta_r = 0, ar = 0, br = 0

    !ニュートン法
    real(8) :: error = 0.0001, phi = 0.50d0, phi2 = 0
    real(8) :: theta_1 = Pi/2.0d0, theta_2 = 0
    real(8) :: func = 0, func_dash = 0

    character(len = 1) d_char 

    write(d_char, '(i1.1)') int(d_particle)

    open(40, file = 'phase_data_1_'//d_char)    !寺井ライン
    open(50, file = 'phase_data_2')


    theta_r = 20.0d0 * Pi / 180.0d0    !安息角
    ar = theta_r / (1.0d0 - cos(theta_r))    !寺井ライン
    br = 1.0d0 / sin(theta_r)    !小野ライン

    do xi = 1, Nxy
        phi_eq = dble(xi)/dble(Nxy)    !充填率を変化させる
        Mp = 2880.0d0 * phi_eq * m / d_particle**2
        func = 0
        func_dash = 0

        if (phi_eq == 0.50d0) then    !充填率が0.5のとき，特異点
            Sxy1 = 1.0d0
        else 
            !傾きのとり方によってはSxyが1を超えてしまうので矯正する
            if (phi_eq <= 0.5d0) then 
                Sxy1 = 0.9999d0
            else 
                Sxy1 = 0.001d0
            endif 

            Sxy2 = 0.0d0
            do i =0, 100
                func = Sxy1 - sin(Sxy1*sqrt(1.0d0-Sxy1**2.0d0) + Pi*phi_eq)
                func_dash = 1.0d0 - ((1.0d0 - 2.0d0*Sxy1**2.0d0)/sqrt(1.0d0 - Sxy1**2.0d0))&
                *cos(Sxy1*sqrt(1.0d0-Sxy1**2.0d0) + Pi*phi_eq)
                Sxy2 = Sxy1 - func/func_dash
                if (abs(Sxy2 - Sxy1) < error) exit
                Sxy1 = Sxy2
                if (1.0d0 - Sxy1**2.0d0 < 0) exit
            enddo 
        endif

        ! alpha_eq1 = asin((80.0d0*phi_eq/(Mc + 80.0d0*phi_eq))*(Sxy1**3.0d0)/(3.0d0*Pi*ar*phi_eq))
        ! alpha_eq2 = asin((80.0d0*phi_eq/(Mc + 80.0d0*phi_eq))*(Sxy1**3.0d0)/(3.0d0*Pi*br*phi_eq))
        alpha_eq1 = asin((Mp*Sxy1**3)/(1.50d0*Pi*ar*(1+deltad/R)*(Mc+Mp)*phi_eq))

        if (1.0d0 - Sxy1**2.0d0 >= 0) then
            write(40,*) alpha_eq1 * 180.0d0 / Pi, phi_eq
            write(50,*) alpha_eq2 * 180.0d0 / Pi, phi_eq
        endif
    enddo

    close(40)
    close(50)
end program cylinder
