#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#![register_tool(c2rust)]
#![feature(register_tool)]
extern "C" {
    fn rand() -> libc::c_int;
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(__ptr: *mut libc::c_void);
}
pub type int64_t = libc::c_long;
#[no_mangle]
pub unsafe extern "C" fn max_(
    mut a: libc::c_double,
    mut b: libc::c_double,
) -> libc::c_double {
    if a >= b { return a } else { return b };
}
#[no_mangle]
pub unsafe extern "C" fn sign(mut x: libc::c_double) -> libc::c_double {
    if x > 0.0f64 {
        return 1.0f64
    } else if x < 0.0f64 {
        return -1.0f64
    } else {
        return 0.0f64
    };
}
#[no_mangle]
pub unsafe extern "C" fn sum_int(mut x: *mut int64_t, mut n: int64_t) -> int64_t {
    let mut sum: int64_t = 0 as libc::c_int as int64_t;
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n {
        sum += *x.offset(i as isize);
        i += 1;
    }
    return sum;
}
#[no_mangle]
pub unsafe extern "C" fn sum_double(
    mut x: *mut libc::c_double,
    mut n: int64_t,
) -> libc::c_double {
    let mut sum: libc::c_double = 0 as libc::c_int as libc::c_double;
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n {
        sum += *x.offset(i as isize);
        i += 1;
    }
    return sum;
}
#[no_mangle]
pub unsafe extern "C" fn swap(
    mut x: *mut libc::c_double,
    mut i: int64_t,
    mut j: int64_t,
) {
    let mut temp: libc::c_double = *x.offset(i as isize);
    *x.offset(i as isize) = *x.offset(j as isize);
    *x.offset(j as isize) = temp;
}
#[no_mangle]
pub unsafe extern "C" fn swap_2d(
    mut x: *mut *mut libc::c_double,
    mut n2: int64_t,
    mut i: int64_t,
    mut j: int64_t,
) {
    let mut temp: libc::c_double = 0.;
    let mut k: int64_t = 0 as libc::c_int as int64_t;
    while k < n2 {
        temp = *(*x.offset(i as isize)).offset(k as isize);
        *(*x.offset(i as isize))
            .offset(k as isize) = *(*x.offset(j as isize)).offset(k as isize);
        *(*x.offset(j as isize)).offset(k as isize) = temp;
        k += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn compare_ascending(
    mut i: *const libc::c_void,
    mut j: *const libc::c_void,
) -> libc::c_int {
    let mut a: libc::c_double = *(i as *mut libc::c_double);
    let mut b: libc::c_double = *(j as *mut libc::c_double);
    if a > b {
        return 1 as libc::c_int
    } else if a < b {
        return -(1 as libc::c_int)
    } else {
        return 0 as libc::c_int
    };
}
#[no_mangle]
pub unsafe extern "C" fn compare_descending(
    mut i: *const libc::c_void,
    mut j: *const libc::c_void,
) -> libc::c_int {
    let mut a: libc::c_double = *(i as *mut libc::c_double);
    let mut b: libc::c_double = *(j as *mut libc::c_double);
    if a > b {
        return -(1 as libc::c_int)
    } else if a < b {
        return 1 as libc::c_int
    } else {
        return 0 as libc::c_int
    };
}
#[no_mangle]
pub unsafe extern "C" fn compare_ascending_2d(
    mut i: *const libc::c_void,
    mut j: *const libc::c_void,
    mut k: int64_t,
) -> libc::c_int {
    let mut a: *mut libc::c_double = *(i as *mut *mut libc::c_double);
    let mut b: *mut libc::c_double = *(j as *mut *mut libc::c_double);
    if *a.offset(k as isize) > *b.offset(k as isize) {
        return 1 as libc::c_int
    } else if *a.offset(k as isize) < *b.offset(k as isize) {
        return -(1 as libc::c_int)
    } else {
        return 0 as libc::c_int
    };
}
#[no_mangle]
pub unsafe extern "C" fn compare_descending_2d(
    mut i: *const libc::c_void,
    mut j: *const libc::c_void,
    mut k: int64_t,
) -> libc::c_int {
    let mut a: *mut libc::c_double = *(i as *mut *mut libc::c_double);
    let mut b: *mut libc::c_double = *(j as *mut *mut libc::c_double);
    if *a.offset(k as isize) > *b.offset(k as isize) {
        return -(1 as libc::c_int)
    } else if *a.offset(k as isize) < *b.offset(k as isize) {
        return 1 as libc::c_int
    } else {
        return 0 as libc::c_int
    };
}
#[no_mangle]
pub unsafe extern "C" fn random_0_to_1() -> libc::c_double {
    return rand() as libc::c_double / 2147483647 as libc::c_int as libc::c_double;
}
#[no_mangle]
pub unsafe extern "C" fn random_range(mut begin: int64_t, mut end: int64_t) -> int64_t {
    return begin
        + (random_0_to_1()
            * (end - begin + 1 as libc::c_int as libc::c_long) as libc::c_double)
            as int64_t;
}
#[no_mangle]
pub unsafe extern "C" fn fill_array_int(
    mut x: *mut int64_t,
    mut n: int64_t,
    mut value: int64_t,
) {
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n {
        *x.offset(i as isize) = value;
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn copy_array_int(
    mut x: *mut int64_t,
    mut y: *mut int64_t,
    mut n: int64_t,
) {
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n {
        *y.offset(i as isize) = *x.offset(i as isize);
        i += 1;
    }
}
#[no_mangle]
pub unsafe extern "C" fn zip(
    mut array_a: *mut libc::c_double,
    mut array_b: *mut libc::c_double,
    mut n: int64_t,
) -> *mut *mut libc::c_double {
    let mut zip_arr: *mut *mut libc::c_double = malloc(
        (n as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<*mut libc::c_double>() as libc::c_ulong),
    ) as *mut *mut libc::c_double;
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n {
        let ref mut fresh0 = *zip_arr.offset(i as isize);
        *fresh0 = malloc(
            (2 as libc::c_int as libc::c_ulong)
                .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
        ) as *mut libc::c_double;
        *(*zip_arr.offset(i as isize))
            .offset(0 as libc::c_int as isize) = *array_a.offset(i as isize);
        *(*zip_arr.offset(i as isize))
            .offset(1 as libc::c_int as isize) = *array_b.offset(i as isize);
        i += 1;
    }
    return zip_arr;
}
#[no_mangle]
pub unsafe extern "C" fn free_zip_memory(
    mut zip_array: *mut *mut libc::c_double,
    mut n: int64_t,
) {
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n {
        free(*zip_array.offset(i as isize) as *mut libc::c_void);
        i += 1;
    }
    free(zip_array as *mut libc::c_void);
}
#[no_mangle]
pub unsafe extern "C" fn partition_on_value(
    mut x: *mut libc::c_double,
    mut begin: int64_t,
    mut end: int64_t,
    mut value: libc::c_double,
) -> int64_t {
    let mut i: int64_t = begin;
    let mut j: int64_t = begin;
    while i <= end {
        if *x.offset(i as isize) < value {
            swap(x, j, i);
            j += 1;
        }
        i += 1;
    }
    return j;
}
#[no_mangle]
pub unsafe extern "C" fn partition_on_kth_element(
    mut x: *mut libc::c_double,
    mut begin: int64_t,
    mut end: int64_t,
    mut k: int64_t,
) -> int64_t {
    let mut value: libc::c_double = *x.offset(k as isize);
    swap(x, k, end);
    let mut i: int64_t = begin;
    let mut j: int64_t = begin;
    while j < end {
        if *x.offset(j as isize) < value {
            swap(x, i, j);
            i += 1;
        }
        j += 1;
    }
    swap(x, i, end);
    return i;
}
#[no_mangle]
pub unsafe extern "C" fn partition_on_kth_element_2d(
    mut x: *mut *mut libc::c_double,
    mut begin: int64_t,
    mut end: int64_t,
    mut n2: int64_t,
    mut m: int64_t,
    mut k: int64_t,
) -> int64_t {
    let mut value: libc::c_double = *(*x.offset(k as isize)).offset(m as isize);
    swap_2d(x, n2, k, end);
    let mut i: int64_t = begin;
    let mut j: int64_t = begin;
    while j < end {
        if *(*x.offset(j as isize)).offset(m as isize) < value {
            swap_2d(x, n2, i, j);
            i += 1;
        }
        j += 1;
    }
    swap_2d(x, n2, i, end);
    return i;
}
#[no_mangle]
pub unsafe extern "C" fn partition_on_kth_smallest(
    mut x: *mut libc::c_double,
    mut begin: int64_t,
    mut end: int64_t,
    mut k: int64_t,
) -> libc::c_double {
    loop {
        if begin == end {
            return *x.offset(begin as isize);
        }
        let mut pivot_index: int64_t = random_range(begin, end);
        pivot_index = partition_on_kth_element(x, begin, end, pivot_index);
        if k == pivot_index {
            return *x.offset(k as isize)
        } else {
            if k < pivot_index {
                end = pivot_index - 1 as libc::c_int as libc::c_long;
            } else {
                begin = pivot_index + 1 as libc::c_int as libc::c_long;
            }
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn partition_on_kth_smallest_2d(
    mut x: *mut *mut libc::c_double,
    mut begin: int64_t,
    mut end: int64_t,
    mut n2: int64_t,
    mut m: int64_t,
    mut k: int64_t,
) -> libc::c_double {
    loop {
        if begin == end {
            return *(*x.offset(begin as isize)).offset(m as isize);
        }
        let mut pivot_index: int64_t = random_range(begin, end);
        pivot_index = partition_on_kth_element_2d(x, begin, end, n2, m, pivot_index);
        if k == pivot_index {
            return *(*x.offset(k as isize)).offset(m as isize)
        } else {
            if k < pivot_index {
                end = pivot_index - 1 as libc::c_int as libc::c_long;
            } else {
                begin = pivot_index + 1 as libc::c_int as libc::c_long;
            }
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn select_kth_smallest(
    mut x: *mut libc::c_double,
    mut n: int64_t,
    mut k: int64_t,
) -> libc::c_double {
    let mut x_copy: *mut libc::c_double = malloc(
        (n as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    ) as *mut libc::c_double;
    let mut i: int64_t = 0 as libc::c_int as int64_t;
    while i < n {
        *x_copy.offset(i as isize) = *x.offset(i as isize);
        i += 1;
    }
    let mut kth_smallest: libc::c_double = partition_on_kth_smallest(
        x_copy,
        0 as libc::c_int as int64_t,
        n - 1 as libc::c_int as libc::c_long,
        k,
    );
    free(x_copy as *mut libc::c_void);
    return kth_smallest;
}
