program yao
    use iso_fortran_env, only: dp => real64
    implicit none
    real(dp),parameter :: pi = acos(-1.0_dp)
    real(dp) :: t, t_jd_j2000, EO
    real(dp) :: a_icrs_ra, b_icrs_ra, c_icrs_ra, d_icrs_ra
    real(dp) :: a_icrs_dec, b_icrs_dec, c_icrs_dec, d_icrs_dec
    real(dp) :: a_cirs_ra, b_cirs_ra, c_cirs_ra, d_cirs_ra
    real(dp) :: a_cirs_dec, b_cirs_dec, c_cirs_dec, d_cirs_dec
    a_icrs_ra = 09.0_dp+27.0_dp/60.0_dp+35.24270/3600.0_dp
    b_icrs_ra = 16.0_dp+29.0_dp/60.0_dp+24.45970/3600.0_dp
    c_icrs_ra = 21.0_dp+31.0_dp/60.0_dp+33.5317148/3600.0_dp
    d_icrs_ra = 03.0_dp+46.0_dp/60.0_dp+24.2/3600.0_dp
    a_icrs_dec = -(08.0_dp+39.0_dp/60.0_dp+30.9583/3600.0_dp)
    b_icrs_dec = -(26.0_dp+25.0_dp/60.0_dp+55.2094/3600.0_dp)
    c_icrs_dec = -(05.0_dp+34.0_dp/60.0_dp+16.232006/3600.0_dp)
    d_icrs_dec = +(24.0_dp+06.0_dp/60.0_dp+50.0/3600.0_dp)
    a_icrs_ra = a_icrs_ra*(pi/12.0_dp)
    b_icrs_ra = b_icrs_ra*(pi/12.0_dp)
    c_icrs_ra = c_icrs_ra*(pi/12.0_dp)
    d_icrs_ra = d_icrs_ra*(pi/12.0_dp)
    a_icrs_dec = a_icrs_dec*(pi/180.0_dp)
    b_icrs_dec = b_icrs_dec*(pi/180.0_dp)
    c_icrs_dec = c_icrs_dec*(pi/180.0_dp)
    d_icrs_dec = d_icrs_dec*(pi/180.0_dp)
    t = 2000.0_dp
    t_jd_j2000 = (t-2000.0_dp)*365.25_dp-0.5_dp
    call iau_ATCI13(a_icrs_ra, a_icrs_dec, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & 
                    2451545.0_dp, t_jd_j2000, a_cirs_ra, a_cirs_dec, EO)
    call iau_ATCI13(b_icrs_ra, b_icrs_dec, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & 
                    2451545.0_dp, t_jd_j2000, b_cirs_ra, b_cirs_dec, EO)
    call iau_ATCI13(c_icrs_ra, c_icrs_dec, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & 
                    2451545.0_dp, t_jd_j2000, c_cirs_ra, c_cirs_dec, EO)
    call iau_ATCI13(d_icrs_ra, d_icrs_dec, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & 
                    2451545.0_dp, t_jd_j2000, d_cirs_ra, d_cirs_dec, EO)
end program yao
