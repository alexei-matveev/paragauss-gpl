#define UPPERCASE	1
#define UPPERCASE_	2
#define lowercase	3
#define lowercase_	4

#ifndef F77_EXT_NAMES
# define F77_EXT_NAMES	lowercase_
#endif

#if F77_EXT_NAMES == lowercase
#define F77_pkbyte_scalar	pvm_pkbyte_scalar
#define F77_pkstring		pvm_pkstring
#define F77_pkcplx_scalar	pvm_pkcplx_scalar
#define F77_pkdcplx_scalar	pvm_pkdcplx_scalar
#define F77_pkdouble_scalar	pvm_pkdouble_scalar
#define F77_pkfloat_scalar	pvm_pkfloat_scalar
#define F77_pkint_scalar	pvm_pkint_scalar
#define F77_pkshort_scalar	pvm_pkshort_scalar
#define F77_upkbyte_scalar	pvm_upkbyte_scalar
#define F77_upkstring		pvm_upkstring
#define F77_upkcplx_scalar	pvm_upkcplx_scalar
#define F77_upkdcplx_scalar	pvm_upkdcplx_scalar
#define F77_upkdouble_scalar	pvm_upkdouble_scalar
#define F77_upkfloat_scalar	pvm_upkfloat_scalar
#define F77_upkint_scalar	pvm_upkint_scalar
#define F77_upkshort_scalar	pvm_upkshort_scalar
#define F77_pkbyte_vec		pvm_pkbyte_vec
#define F77_pkcplx_vec		pvm_pkcplx_vec
#define F77_pkdcplx_vec		pvm_pkdcplx_vec
#define F77_pkdouble_vec	pvm_pkdouble_vec
#define F77_pkfloat_vec		pvm_pkfloat_vec
#define F77_pkint_vec		pvm_pkint_vec
#define F77_pkshort_vec		pvm_pkshort_vec
#define F77_upkbyte_vec		pvm_upkbyte_vec
#define F77_upkcplx_vec		pvm_upkcplx_vec
#define F77_upkdcplx_vec	pvm_upkdcplx_vec
#define F77_upkdouble_vec	pvm_upkdouble_vec
#define F77_upkfloat_vec	pvm_upkfloat_vec
#define F77_upkint_vec		pvm_upkint_vec
#define F77_upkshort_vec	pvm_upkshort_vec
#define F77_pkbyte_vecsc	pvm_pkbyte_vecsc
#define F77_pkcplx_vecsc	pvm_pkcplx_vecsc
#define F77_pkdcplx_vecsc	pvm_pkdcplx_vecsc
#define F77_pkdouble_vecsc	pvm_pkdouble_vecsc
#define F77_pkfloat_vecsc	pvm_pkfloat_vecsc
#define F77_pkint_vecsc		pvm_pkint_vecsc
#define F77_pkshort_vecsc	pvm_pkshort_vecsc
#define F77_upkbyte_vecsc	pvm_upkbyte_vecsc
#define F77_upkcplx_vecsc	pvm_upkcplx_vecsc
#define F77_upkdcplx_vecsc	pvm_upkdcplx_vecsc
#define F77_upkdouble_vecsc	pvm_upkdouble_vecsc
#define F77_upkfloat_vecsc	pvm_upkfloat_vecsc
#define F77_upkint_vecsc	pvm_upkint_vecsc
#define F77_upkshort_vecsc	pvm_upkshort_vecsc

#elif F77_EXT_NAMES == lowercase_
#define F77_pkbyte_scalar	pvm_pkbyte_scalar_
#define F77_pkstring		pvm_pkstring_
#define F77_pkcplx_scalar	pvm_pkcplx_scalar_
#define F77_pkdcplx_scalar	pvm_pkdcplx_scalar_
#define F77_pkdouble_scalar	pvm_pkdouble_scalar_
#define F77_pkfloat_scalar	pvm_pkfloat_scalar_
#define F77_pkint_scalar	pvm_pkint_scalar_
#define F77_pkshort_scalar	pvm_pkshort_scalar_
#define F77_upkbyte_scalar	pvm_upkbyte_scalar_
#define F77_upkstring		pvm_upkstring_
#define F77_upkcplx_scalar	pvm_upkcplx_scalar_
#define F77_upkdcplx_scalar	pvm_upkdcplx_scalar_
#define F77_upkdouble_scalar	pvm_upkdouble_scalar_
#define F77_upkfloat_scalar	pvm_upkfloat_scalar_
#define F77_upkint_scalar	pvm_upkint_scalar_
#define F77_upkshort_scalar	pvm_upkshort_scalar_
#define F77_pkbyte_vec		pvm_pkbyte_vec_
#define F77_pkcplx_vec		pvm_pkcplx_vec_
#define F77_pkdcplx_vec		pvm_pkdcplx_vec_
#define F77_pkdouble_vec	pvm_pkdouble_vec_
#define F77_pkfloat_vec		pvm_pkfloat_vec_
#define F77_pkint_vec		pvm_pkint_vec_
#define F77_pkshort_vec		pvm_pkshort_vec_
#define F77_upkbyte_vec		pvm_upkbyte_vec_
#define F77_upkcplx_vec		pvm_upkcplx_vec_
#define F77_upkdcplx_vec	pvm_upkdcplx_vec_
#define F77_upkdouble_vec	pvm_upkdouble_vec_
#define F77_upkfloat_vec	pvm_upkfloat_vec_
#define F77_upkint_vec		pvm_upkint_vec_
#define F77_upkshort_vec	pvm_upkshort_vec_
#define F77_pkbyte_vecsc	pvm_pkbyte_vecsc_
#define F77_pkcplx_vecsc	pvm_pkcplx_vecsc_
#define F77_pkdcplx_vecsc	pvm_pkdcplx_vecsc_
#define F77_pkdouble_vecsc	pvm_pkdouble_vecsc_
#define F77_pkfloat_vecsc	pvm_pkfloat_vecsc_
#define F77_pkint_vecsc		pvm_pkint_vecsc_
#define F77_pkshort_vecsc	pvm_pkshort_vecsc_
#define F77_upkbyte_vecsc	pvm_upkbyte_vecsc_
#define F77_upkcplx_vecsc	pvm_upkcplx_vecsc_
#define F77_upkdcplx_vecsc	pvm_upkdcplx_vecsc_
#define F77_upkdouble_vecsc	pvm_upkdouble_vecsc_
#define F77_upkfloat_vecsc	pvm_upkfloat_vecsc_
#define F77_upkint_vecsc	pvm_upkint_vecsc_
#define F77_upkshort_vecsc	pvm_upkshort_vecsc_
#elif 1
#error "need extension for UPPERCASE"
#endif
