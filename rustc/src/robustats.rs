#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#![register_tool(c2rust)]
#![feature(register_tool)]
extern "C" {
    fn fabs(_: libc::c_double) -> libc::c_double;
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(__ptr: *mut libc::c_void);
}
use std::iter::zip;
mod my_base;

pub type int64_t = libc::c_long;
#[no_mangle]
pub unsafe extern "C" fn weighted_median(
    mut x: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut begin: int64_t,
    mut end: int64_t,
) -> libc::c_double {
    let mut xw_n: int64_t = 0;
    let mut n: int64_t = 0;
    let mut i: int64_t = 0;
    let mut median_index: int64_t = 0;
    let mut median: libc::c_double = 0.;
    let mut w_lower_sum: libc::c_double = 0.;
    let mut w_lower_sum_norm: libc::c_double = 0.;
    let mut w_higher_sum: libc::c_double = 0.;
    let mut w_higher_sum_norm: libc::c_double = 0.;
    xw_n = end - begin + 1 as libc::c_int as libc::c_long;
    let mut xw: *mut *mut libc::c_double = my_base::zip(x, w, xw_n) as *mut *mut libc::c_double;
    let mut w_sum: libc::c_double = my_base::sum_double(w, xw_n) as libc::c_double;
    loop {
        n = end - begin + 1 as libc::c_int as libc::c_long;
        if n == 1 as libc::c_int as libc::c_long {
            my_base::free_zip_memory(xw, xw_n);
            return *x.offset(begin as isize);
        } else {
            if n == 2 as libc::c_int as libc::c_long {
                my_base::free_zip_memory(xw, xw_n);
                if *w.offset(begin as isize) >= *w.offset(end as isize) {
                    return *x.offset(begin as isize)
                } else {
                    return *x.offset(end as isize)
                }
            } else {
                median_index = begin
                    + (n - 1 as libc::c_int as libc::c_long)
                        / 2 as libc::c_int as libc::c_long;
                median = my_base::partition_on_kth_smallest_2d(
                    xw,
                    begin,
                    end,
                    2 as libc::c_int,
                    0 as libc::c_int,
                    median_index,
                ) as libc::c_double;
                w_lower_sum = 0.0f64;
                i = begin;
                while i < median_index {
                    w_lower_sum
                        += *(*xw.offset(i as isize)).offset(1 as libc::c_int as isize);
                    i += 1;
                }
                w_lower_sum_norm = w_lower_sum / w_sum;
                w_higher_sum = 0.0f64;
                i = median_index + 1 as libc::c_int as libc::c_long;
                while i <= end {
                    w_higher_sum
                        += *(*xw.offset(i as isize)).offset(1 as libc::c_int as isize);
                    i += 1;
                }
                w_higher_sum_norm = w_higher_sum / w_sum;
                if w_lower_sum_norm < 0.5f64 && w_higher_sum_norm < 0.5f64 {
                    my_base::free_zip_memory(xw, xw_n);
                    return median;
                } else {
                    if w_lower_sum_norm > 0.5f64 {
                        *(*xw.offset(median_index as isize))
                            .offset(
                                1 as libc::c_int as isize,
                            ) = *(*xw.offset(median_index as isize))
                            .offset(1 as libc::c_int as isize) + w_higher_sum;
                        end = median_index;
                    } else {
                        *(*xw.offset(median_index as isize))
                            .offset(
                                1 as libc::c_int as isize,
                            ) = *(*xw.offset(median_index as isize))
                            .offset(1 as libc::c_int as isize) + w_lower_sum;
                        begin = median_index;
                    }
                }
            }
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn h_kernel(
    mut i: int64_t,
    mut j: int64_t,
    mut z_plus: *mut libc::c_double,
    mut n_plus: int64_t,
    mut z_minus: *mut libc::c_double,
    mut n_minus: int64_t,
    mut epsilon: libc::c_double,
) -> libc::c_double {
    let mut a: libc::c_double = *z_plus.offset(i as isize);
    let mut b: libc::c_double = *z_minus.offset(j as isize);
    if fabs(a - b) <= 2 as libc::c_int as libc::c_double * epsilon {
        return my_base::sign(n_plus - i - j - 1 as libc::c_int as libc::c_long) as libc::c_double
    } else {
        return (a + b) / (a - b)
    };
}
#[no_mangle]
pub unsafe extern "C" fn where_h_greater_than_u(
    mut p: *mut int64_t,
    mut n_p: int64_t,
    mut z_plus: *mut libc::c_double,
    mut n_plus: int64_t,
    mut z_minus: *mut libc::c_double,
    mut n_minus: int64_t,
    mut u: libc::c_double,
    mut epsilon: libc::c_double,
    mut k_epsilon: libc::c_double,
) {
    my_base::fill_array_int(p, n_p, 0 as libc::c_int);
    let mut h: libc::c_double = 0.;
    let mut j: int64_t = 0 as libc::c_int as int64_t;
    let mut i: int64_t = n_plus - 1 as libc::c_int as libc::c_long;
    while i >= 0 as libc::c_int as libc::c_long {
        h = h_kernel(i, j, z_plus, n_plus, z_minus, n_minus, k_epsilon);
        while j < n_minus && h - u > epsilon {
            j += 1;
            h = h_kernel(i, j, z_plus, n_plus, z_minus, n_minus, k_epsilon);
        }
        *p.offset(i as isize) = j - 1 as libc::c_int as libc::c_long;
        i -= 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn where_h_less_than_u(
    mut q: *mut int64_t,
    mut n_q: int64_t,
    mut z_plus: *mut libc::c_double,
    mut n_plus: int64_t,
    mut z_minus: *mut libc::c_double,
    mut n_minus: int64_t,
    mut u: libc::c_double,
    mut epsilon: libc::c_double,
    mut k_epsilon: libc::c_double,
) {
    my_base::fill_array_int(q, n_q, 0 as libc::c_int);
    let mut h: libc::c_double = 0.;
    let mut j: int64_t = n_minus - 1 as libc::c_int as libc::c_long;
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n_plus {
        h = h_kernel(i, j, z_plus, n_plus, z_minus, n_minus, k_epsilon);
        while j >= 0 as libc::c_int as libc::c_long && h - u < -epsilon {
            j -= 1;
            h = h_kernel(i, j, z_plus, n_plus, z_minus, n_minus, k_epsilon);
        }
        *q.offset(i as isize) = j + 1 as libc::c_int as libc::c_long;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn medcouple(
    mut x: *mut libc::c_double,
    mut n: int64_t,
    mut epsilon1: libc::c_double,
    mut epsilon2: libc::c_double,
) -> libc::c_double {
    let mut i: int64_t = 0;
    let mut j: int64_t = 0;
    if n < 3 as libc::c_int as libc::c_long {
        return 0.0f64;
    }
    let mut median_index: int64_t = n / 2 as libc::c_int as libc::c_long;
    let mut median: libc::c_double = *x.offset(median_index as isize);
    if fabs(*x.offset(0 as libc::c_int as isize) - median)
        < epsilon1 * (epsilon1 + fabs(median))
    {
        return -1.0f64;
    }
    if fabs(*x.offset((n - 1 as libc::c_int as libc::c_long) as isize) - median)
        < epsilon1 * (epsilon1 + fabs(median))
    {
        return 1.0f64;
    }
    let mut scale_factor: libc::c_double = (2 as libc::c_int
        * my_base::max_(
            *x.offset(0 as libc::c_int as isize) - median,
            median - *x.offset((n - 1 as libc::c_int as libc::c_long) as isize),
        )) as libc::c_double;
    let mut lowest_median_index: int64_t = median_index;
    let mut lowest_median: libc::c_double = median;
    while lowest_median == median {
        lowest_median_index += 1;
        lowest_median = *x.offset(lowest_median_index as isize);
    }
    lowest_median_index -= 1;
    let mut n_plus: int64_t = lowest_median_index + 1 as libc::c_int as libc::c_long;
    let mut z_plus: *mut libc::c_double = malloc(
        (n_plus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    ) as *mut libc::c_double;
    i = 0 as libc::c_int as int64_t;
    while i < n_plus {
        *z_plus.offset(i as isize) = (*x.offset(i as isize) - median) / scale_factor;
        i += 1;
    }
    let mut highest_median_index: int64_t = median_index;
    let mut highest_median: libc::c_double = median;
    while highest_median == median {
        highest_median_index -= 1;
        highest_median = *x.offset(highest_median_index as isize);
    }
    highest_median_index += 1;
    let mut n_minus: int64_t = n - highest_median_index;
    let mut z_minus: *mut libc::c_double = malloc(
        (n_minus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    ) as *mut libc::c_double;
    i = 0 as libc::c_int as int64_t;
    while i < n_minus {
        *z_minus
            .offset(
                i as isize,
            ) = (*x.offset((highest_median_index + i) as isize) - median) / scale_factor;
        i += 1;
    }
    let mut left_border: *mut int64_t = malloc(
        (n_plus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<int64_t>() as libc::c_ulong),
    ) as *mut int64_t;
    my_base::fill_array_int(left_border, n_plus, 0 as libc::c_int);
    let mut right_border: *mut int64_t = malloc(
        (n_plus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<int64_t>() as libc::c_ulong),
    ) as *mut int64_t;
    my_base::fill_array_int(right_border, n_plus, n_minus - 1 as libc::c_int as libc::c_long);
    let mut left_total: int64_t = 0 as libc::c_int as int64_t;
    let mut right_total: int64_t = n_minus * n_plus;
    let mut medcouple_index: int64_t = right_total / 2 as libc::c_int as libc::c_long;
    let mut mid_border: int64_t = 0;
    let mut right_tent_total: int64_t = 0;
    let mut left_tent_total: int64_t = 0;
    let mut row_medians: *mut libc::c_double = malloc(
        (n_plus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    ) as *mut libc::c_double;
    let mut weights: *mut libc::c_double = malloc(
        (n_plus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    ) as *mut libc::c_double;
    let mut w_median: libc::c_double = 0.;
    let mut wm_epsilon: libc::c_double = 0.;
    let mut left_border_tent: *mut int64_t = malloc(
        (n_plus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<int64_t>() as libc::c_ulong),
    ) as *mut int64_t;
    let mut right_border_tent: *mut int64_t = malloc(
        (n_plus as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<int64_t>() as libc::c_ulong),
    ) as *mut int64_t;
    while right_total - left_total > n_plus {
        let mut n_middle_indices: int64_t = 0 as libc::c_int as int64_t;
        i = 0 as libc::c_int as int64_t;
        while i < n_plus {
            if *left_border.offset(i as isize) <= *right_border.offset(i as isize) {
                n_middle_indices += 1;
            }
            i += 1;
        }
        j = 0 as libc::c_int as int64_t;
        i = 0 as libc::c_int as int64_t;
        while i < n_plus {
            if *left_border.offset(i as isize) <= *right_border.offset(i as isize) {
                mid_border = (*left_border.offset(i as isize)
                    + *right_border.offset(i as isize))
                    / 2 as libc::c_int as libc::c_long;
                *row_medians
                    .offset(
                        j as isize,
                    ) = h_kernel(
                    i,
                    mid_border,
                    z_plus,
                    n_plus,
                    z_minus,
                    n_minus,
                    epsilon2,
                );
                *weights
                    .offset(
                        j as isize,
                    ) = (*right_border.offset(i as isize)
                    - *left_border.offset(i as isize) + 1 as libc::c_int as libc::c_long)
                    as libc::c_double;
                j += 1;
            }
            i += 1;
        }
        w_median = weighted_median(
            row_medians,
            weights,
            0 as libc::c_int as int64_t,
            n_middle_indices - 1 as libc::c_int as libc::c_long,
        );
        wm_epsilon = epsilon1 * (epsilon1 + fabs(w_median));
        where_h_greater_than_u(
            right_border_tent,
            n_plus,
            z_plus,
            n_plus,
            z_minus,
            n_minus,
            w_median,
            wm_epsilon,
            epsilon2,
        );
        where_h_less_than_u(
            left_border_tent,
            n_plus,
            z_plus,
            n_plus,
            z_minus,
            n_minus,
            w_median,
            wm_epsilon,
            epsilon2,
        );
        right_tent_total = my_base::sum_int(right_border_tent, n_plus) as libc::c_long + n_plus;
        left_tent_total = my_base::sum_int(left_border_tent, n_plus) as int64_t;
        if medcouple_index <= right_tent_total - 1 as libc::c_int as libc::c_long {
            my_base::copy_array_int(right_border_tent, right_border, n_plus);
            right_total = right_tent_total;
        } else if medcouple_index > left_tent_total - 1 as libc::c_int as libc::c_long {
            my_base::copy_array_int(left_border_tent, left_border, n_plus);
            left_total = left_tent_total;
        } else {
            free(row_medians as *mut libc::c_void);
            free(weights as *mut libc::c_void);
            free(right_border_tent as *mut libc::c_void);
            free(left_border_tent as *mut libc::c_void);
            free(z_minus as *mut libc::c_void);
            free(z_plus as *mut libc::c_void);
            free(left_border as *mut libc::c_void);
            free(right_border as *mut libc::c_void);
            return w_median;
        }
    }
    free(row_medians as *mut libc::c_void);
    free(weights as *mut libc::c_void);
    free(right_border_tent as *mut libc::c_void);
    free(left_border_tent as *mut libc::c_void);
    let mut n_remaining: int64_t = 0 as libc::c_int as int64_t;
    i = 0 as libc::c_int as int64_t;
    while i < n_plus {
        j = *left_border.offset(i as isize);
        while j <= *right_border.offset(i as isize) {
            n_remaining += 1;
            j += 1;
        }
        i += 1;
    }
    let mut remaining: *mut libc::c_double = malloc(
        (n_remaining as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    ) as *mut libc::c_double;
    let mut k: int64_t = 0 as libc::c_int as int64_t;
    i = 0 as libc::c_int as int64_t;
    while i < n_plus {
        j = *left_border.offset(i as isize);
        while j <= *right_border.offset(i as isize) {
            *remaining
                .offset(
                    k as isize,
                ) = -h_kernel(i, j, z_plus, n_plus, z_minus, n_minus, epsilon2);
            k += 1;
            j += 1;
        }
        i += 1;
    }
    let mut medcouple_: libc::c_double = -my_base::select_kth_smallest(
        remaining,
        n_remaining,
        medcouple_index - left_total,
    ) as libc::c_double;
    free(z_minus as *mut libc::c_void);
    free(z_plus as *mut libc::c_void);
    free(left_border as *mut libc::c_void);
    free(right_border as *mut libc::c_void);
    free(remaining as *mut libc::c_void);
    return medcouple_;
}
#[no_mangle]
pub unsafe extern "C" fn mode(
    mut x: *mut libc::c_double,
    mut n: int64_t,
) -> libc::c_double {
    let mut m: int64_t = 0;
    let mut m_half: int64_t = 0;
    let mut i: int64_t = 0;
    let mut j: int64_t = 0;
    let mut width: libc::c_double = 0.;
    let mut min_width: libc::c_double = 0.;
    let mut begin: int64_t = 0 as libc::c_int as int64_t;
    let mut end: int64_t = n - 1 as libc::c_int as libc::c_long;
    loop {
        m = end - begin + 1 as libc::c_int as libc::c_long;
        if m == 1 as libc::c_int as libc::c_long {
            return *x.offset(begin as isize)
        } else {
            if m == 2 as libc::c_int as libc::c_long {
                return (*x.offset(begin as isize) + *x.offset(end as isize)) / 2.0f64
            } else {
                if m == 3 as libc::c_int as libc::c_long {
                    if *x.offset((begin + 1 as libc::c_int as libc::c_long) as isize)
                        - *x.offset(begin as isize)
                        < *x.offset(end as isize)
                            - *x
                                .offset((begin + 1 as libc::c_int as libc::c_long) as isize)
                    {
                        return (*x.offset(begin as isize)
                            + *x
                                .offset(
                                    (begin + 1 as libc::c_int as libc::c_long) as isize,
                                )) / 2.0f64
                    } else if *x
                        .offset((begin + 1 as libc::c_int as libc::c_long) as isize)
                        - *x.offset(begin as isize)
                        > *x.offset(end as isize)
                            - *x
                                .offset((begin + 1 as libc::c_int as libc::c_long) as isize)
                    {
                        return (*x
                            .offset((begin + 1 as libc::c_int as libc::c_long) as isize)
                            + *x.offset(end as isize)) / 2.0f64
                    } else {
                        return *x
                            .offset((begin + 1 as libc::c_int as libc::c_long) as isize)
                    }
                } else {
                    min_width = *x.offset(end as isize) - *x.offset(begin as isize);
                    m_half = (m + 1 as libc::c_int as libc::c_long)
                        / 2 as libc::c_int as libc::c_long;
                    j = begin;
                    i = begin;
                    while i <= begin + m - m_half {
                        width = *x
                            .offset(
                                (i + m_half - 1 as libc::c_int as libc::c_long) as isize,
                            ) - *x.offset(i as isize);
                        if width < min_width {
                            min_width = width;
                            j = i;
                        }
                        i += 1;
                    }
                    begin = j;
                    end = j + m_half - 1 as libc::c_int as libc::c_long;
                }
            }
        }
    };
}
