subroutine icrs2cirs(t, ra_icrs, dec_icrs, ra_cirs, dec_cirs)
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp),intent(in) :: t
    real(dp),intent(in) :: ra_icrs, dec_icrs
    real(dp),intent(out) :: ra_cirs, dec_cirs
    real(dp),parameter :: pi = acos(-1.0_dp)
    real(dp) :: t_jd_j2000, EO
    real(dp) :: ra_icrs_rad, dec_icrs_rad
    t_jd_j2000 = (t-2000.0_dp)*365.25_dp-0.5_dp
    ra_icrs_rad = ra_icrs/12.0_dp*pi
    dec_icrs_rad = dec_icrs/180.0_dp*pi
    call iau_ATCI13(ra_icrs_rad, dec_icrs_rad, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                    2451545.0_dp, t_jd_j2000, ra_cirs, dec_cirs, EO)
end subroutine icrs2cirs

subroutine ecli2cirs(t, ra_ecli, dec_ecli, ra_cirs, dec_cirs)
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp),intent(in) :: t
    real(dp),intent(in) :: ra_ecli, dec_ecli
    real(dp),intent(out) :: ra_cirs, dec_cirs
    real(dp),parameter :: pi = acos(-1.0_dp)
    real(dp) :: t_jd_j2000, EO
    real(dp) :: ra_ecli_rad, dec_ecli_rad
    real(dp) :: ra_icrs_rad, dec_icrs_rad
    t_jd_j2000 = (t-2000.0_dp)*365.25_dp-0.5_dp
    ra_ecli_rad = ra_ecli/12.0_dp*pi
    dec_ecli_rad = dec_ecli/180.0_dp*pi
    call iau_ECEQ06(2451545.0_dp, t_jd_j2000, ra_ecli_rad, dec_ecli_rad, &
                    ra_icrs_rad, dec_icrs_rad)
    call iau_ATCI13(ra_icrs_rad, dec_icrs_rad, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                    2451545.0_dp, t_jd_j2000, ra_cirs, dec_cirs, EO)
end subroutine ecli2cirs

subroutine ecli2cirs_lt(t, ra_ecli, dec_ecli, ra_cirs, dec_cirs)
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp),intent(in) :: t
    real(dp),intent(in) :: ra_ecli, dec_ecli
    real(dp),intent(out) :: ra_cirs, dec_cirs
    double precision :: iau_EPJ
    real(dp),parameter :: pi = acos(-1.0_dp)
    real(dp) :: t_jd_j2000, t_jep, EO
    real(dp) :: ra_ecli_rad, dec_ecli_rad
    real(dp) :: ra_icrs_rad, dec_icrs_rad
    t_jd_j2000 = (t-2000.0_dp)*365.25_dp-0.5_dp
    t_jep = iau_EPJ(2451545.0_dp, t_jd_j2000)
    ra_ecli_rad = ra_ecli/12.0_dp*pi
    dec_ecli_rad = dec_ecli/180.0_dp*pi
    call iau_LTECEQ(t_jep, ra_ecli_rad, dec_ecli_rad, &
                    ra_icrs_rad, dec_icrs_rad)
    call iau_ATCI13(ra_icrs_rad, dec_icrs_rad, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                    2451545.0_dp, t_jd_j2000, ra_cirs, dec_cirs, EO)
end subroutine ecli2cirs_lt

subroutine altaz(theta_2, theta_3, sa_ra, sa_dec, a_ra, a_dec, &
                 alt_a, az_a)
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp),intent(in) :: theta_2, theta_3
    real(dp),intent(in) :: sa_ra, sa_dec
    real(dp),intent(in) :: a_ra, a_dec
    real(dp),intent(out) :: alt_a, az_a
    real(dp),parameter :: pi = acos(-1.0_dp)
    real(dp) :: theta_1
    real(dp) :: x, y, z
    theta_1 = sa_ra &
             -acos(sin(sa_dec)*sin(theta_2)-(sin(theta_3)) &
                  /(cos(sa_dec)*cos(theta_2)))
    x = sin(a_dec)*cos(theta_2) &
       +cos(a_dec)*cos(a_ra-theta_1)*sin(theta_2)
    y = cos(a_dec)*sin(a_ra-theta_1)
    z = sin(a_dec)*sin(theta_2) &
       -cos(a_dec)*cos(a_ra-theta_1)*cos(theta_2)
    alt_a = atan2(y, x)/pi*180.0_dp
    az_a = asin(z)/pi*180.0_dp
end subroutine altaz

program yao
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp),parameter :: pi = acos(-1.0_dp)
    integer :: i
    real(dp) :: t
    real(dp) :: theta_2, theta_3
    real(dp) :: a_icrs_ra, b_icrs_ra, c_icrs_ra, d_icrs_ra
    real(dp) :: a_icrs_dec, b_icrs_dec, c_icrs_dec, d_icrs_dec
    real(dp) :: a_cirs_ra, b_cirs_ra, c_cirs_ra, d_cirs_ra
    real(dp) :: a_cirs_dec, b_cirs_dec, c_cirs_dec, d_cirs_dec
    real(dp) :: sa_ra, sb_ra, sc_ra, sd_ra
    real(dp) :: sa_dec, sb_dec, sc_dec, sd_dec
    real(dp) :: sa_cirs_ra, sb_cirs_ra, sc_cirs_ra, sd_cirs_ra
    real(dp) :: sa_cirs_dec, sb_cirs_dec, sc_cirs_dec, sd_cirs_dec
    real(dp) :: alt_a, alt_b, alt_c, alt_d
    real(dp) :: az_a, az_b, az_c, az_d
    do i = 2000, -3000, -100
    t = i
    a_icrs_ra = 09.0_dp+27.0_dp/60.0_dp+35.24270/3600.0_dp
    b_icrs_ra = 16.0_dp+29.0_dp/60.0_dp+24.45970/3600.0_dp
    c_icrs_ra = 21.0_dp+31.0_dp/60.0_dp+33.5317148/3600.0_dp
    d_icrs_ra = 03.0_dp+46.0_dp/60.0_dp+24.2/3600.0_dp
    a_icrs_dec = -(08.0_dp+39.0_dp/60.0_dp+30.9583/3600.0_dp)
    b_icrs_dec = -(26.0_dp+25.0_dp/60.0_dp+55.2094/3600.0_dp)
    c_icrs_dec = -(05.0_dp+34.0_dp/60.0_dp+16.232006/3600.0_dp)
    d_icrs_dec = +(24.0_dp+06.0_dp/60.0_dp+50.0/3600.0_dp)
    call icrs2cirs(t, a_icrs_ra, a_icrs_dec, a_cirs_ra, a_cirs_dec)
    call icrs2cirs(t, b_icrs_ra, b_icrs_dec, b_cirs_ra, b_cirs_dec)
    call icrs2cirs(t, c_icrs_ra, c_icrs_dec, c_cirs_ra, c_cirs_dec)
    call icrs2cirs(t, d_icrs_ra, d_icrs_dec, d_cirs_ra, d_cirs_dec)
    !print *, a_cirs_ra, b_cirs_ra, c_cirs_ra, d_cirs_ra
    !print *, a_cirs_dec, b_cirs_dec, c_cirs_dec, d_cirs_dec
    sa_ra = 00.0_dp
    sb_ra = 06.0_dp
    sc_ra = 12.0_dp
    sd_ra = 18.0_dp
    sa_dec = 0.0_dp
    sb_dec = 0.0_dp
    sc_dec = 0.0_dp
    sd_dec = 0.0_dp
    call ecli2cirs_lt(t, sa_ra, sa_dec, sa_cirs_ra, sa_cirs_dec)
    call ecli2cirs_lt(t, sb_ra, sb_dec, sb_cirs_ra, sb_cirs_dec)
    call ecli2cirs_lt(t, sc_ra, sc_dec, sc_cirs_ra, sc_cirs_dec)
    call ecli2cirs_lt(t, sd_ra, sd_dec, sd_cirs_ra, sd_cirs_dec)
    !print *, sa_cirs_ra, sb_cirs_ra, sc_cirs_ra, sd_cirs_ra
    !print *, sa_cirs_dec, sb_cirs_dec, sc_cirs_dec, sd_cirs_dec
    call ecli2cirs(t, sa_ra, sa_dec, sa_cirs_ra, sa_cirs_dec)
    call ecli2cirs(t, sb_ra, sb_dec, sb_cirs_ra, sb_cirs_dec)
    call ecli2cirs(t, sc_ra, sc_dec, sc_cirs_ra, sc_cirs_dec)
    call ecli2cirs(t, sd_ra, sd_dec, sd_cirs_ra, sd_cirs_dec)
    !print *, sa_cirs_ra, sb_cirs_ra, sc_cirs_ra, sd_cirs_ra
    !print *, sa_cirs_dec, sb_cirs_dec, sc_cirs_dec, sd_cirs_dec
    theta_2 = +35.0_dp/180.0_dp*pi
    theta_3 = -15.0_dp/180.0_dp*pi
    call altaz(theta_2, theta_3, sa_cirs_ra, sa_cirs_dec, &
               a_cirs_ra, a_cirs_dec, alt_a, az_a)
    call altaz(theta_2, theta_3, sb_cirs_ra, sb_cirs_dec, &
               b_cirs_ra, b_cirs_dec, alt_b, az_b)
    call altaz(theta_2, theta_3, sc_cirs_ra, sc_cirs_dec, &
               c_cirs_ra, c_cirs_dec, alt_c, az_c)
    call altaz(theta_2, theta_3, sd_cirs_ra, sd_cirs_dec, &
               d_cirs_ra, d_cirs_dec, alt_d, az_d)
    print *, i
    print *, alt_a, alt_b, alt_c, alt_d
    print *, az_a, az_b, az_c, az_d
    end do
end program yao
