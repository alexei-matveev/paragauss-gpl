#define UPPERCASE       1
#define UPPERCASE_      2
#define lowercase       3
#define lowercase_      4

#ifndef F77_EXT_NAMES
# ifdef	__SR8000
#  define F77_EXT_NAMES		UPPERCASE
# else
#  define F77_EXT_NAMES		lowercase_
# endif
#endif

#if F77_EXT_NAMES == UPPERCASE
#define	comm_send		COMM_SEND
#define	comm_msgtag_buf		COMM_MSGTAG_BUF
#define	mpix_pkint_scalar	MPIX_PKINT_SCALAR
#define	mpix_pkint_vecsc		MPIX_PKINT_VECSC
#define	mpix_upkint_vecsc	MPIX_UPKINT_VECSC
#define	mpix_upkdouble_vecsc	MPIX_UPKDOUBLE_VECSC
#define	mpix_upkdouble_scalar	MPIX_UPKDOUBLE_SCALAR
#define	comm_sending_host_buf	COMM_SENDING_HOST_BUF
#define	mpix_pkdouble_scalar	MPIX_PKDOUBLE_SCALAR
#define	mpix_upkint_scalar	MPIX_UPKINT_SCALAR
#define	mpix_pkint_vec		MPIX_PKINT_VEC
#define	mpix_upkint_vec		MPIX_UPKINT_VEC
#define	mpix_pkdouble_vec	MPIX_PKDOUBLE_VEC
#define	mpix_upkdouble_vec	MPIX_UPKDOUBLE_VEC
#define	mpix_pkdouble_vecsc	MPIX_PKDOUBLE_VECSC
#define	comm_init_buffer_data	COMM_INIT_BUFFER_DATA
#define	comm_init_send_buf	COMM_INIT_SEND_BUF
#define	comm_save_recv_c	COMM_SAVE_RECV_C
#define	comm_save_recv_nonblocking_buf	COMM_SAVE_RECV_NONBLOCKING_BUF
#define	mpix_pkshort_scalar	MPIX_PKSHORT_SCALAR
#define	mpix_upkshort_scalar	MPIX_UPKSHORT_SCALAR
#define	mpix_pkshort_vec		MPIX_PKSHORT_VEC
#define	mpix_upkshort_vec	MPIX_UPKSHORT_VEC
#define	mpix_pkstring		MPIX_PKSTRING
#define	mpix_upkstring		MPIX_UPKSTRING
#elif F77_EXT_NAMES == lowercase
#define	comm_send		comm_send
#define	comm_msgtag_buf		comm_msgtag_buf
#define	mpix_pkint_scalar	mpix_pkint_scalar
#define	mpix_pkint_vecsc		mpix_pkint_vecsc
#define	mpix_upkint_vecsc	mpix_upkint_vecsc
#define	mpix_upkdouble_vecsc	mpix_upkdouble_vecsc
#define	mpix_upkdouble_scalar	mpix_upkdouble_scalar
#define	comm_sending_host_buf	comm_sending_host_buf
#define	mpix_pkdouble_scalar	mpix_pkdouble_scalar
#define	mpix_upkint_scalar	mpix_upkint_scalar
#define	mpix_pkint_vec		mpix_pkint_vec
#define	mpix_upkint_vec		mpix_upkint_vec
#define	mpix_pkdouble_vec	mpix_pkdouble_vec
#define	mpix_upkdouble_vec	mpix_upkdouble_vec
#define	mpix_pkdouble_vecsc	mpix_pkdouble_vecsc
#define	comm_init_buffer_data	comm_init_buffer_data
#define	comm_init_send_buf	comm_init_send_buf
#define	comm_save_recv_c	comm_save_recv_c
#define	comm_save_recv_nonblocking_buf	comm_save_recv_nonblocking_buf
#define	mpix_pkshort_scalar	mpix_pkshort_scalar
#define	mpix_upkshort_scalar	mpix_upkshort_scalar
#define	mpix_pkshort_vec		mpix_pkshort_vec
#define	mpix_upkshort_vec	mpix_upkshort_vec
#define	mpix_pkstring		mpix_pkstring
#define	mpix_upkstring		mpix_upkstring
#elif F77_EXT_NAMES == lowercase_
#define	comm_send		comm_send_
#define	comm_msgtag_buf		comm_msgtag_buf_
#define	mpix_pkint_scalar	mpix_pkint_scalar_
#define	mpix_pkint_vecsc		mpix_pkint_vecsc_
#define	mpix_upkint_vecsc	mpix_upkint_vecsc_
#define	mpix_upkdouble_vecsc	mpix_upkdouble_vecsc_
#define	mpix_upkdouble_scalar	mpix_upkdouble_scalar_
#define	comm_sending_host_buf	comm_sending_host_buf_
#define	mpix_pkdouble_scalar	mpix_pkdouble_scalar_
#define	mpix_upkint_scalar	mpix_upkint_scalar_
#define	mpix_pkint_vec		mpix_pkint_vec_
#define	mpix_upkint_vec		mpix_upkint_vec_
#define	mpix_pkdouble_vec	mpix_pkdouble_vec_
#define	mpix_upkdouble_vec	mpix_upkdouble_vec_
#define	mpix_pkdouble_vecsc	mpix_pkdouble_vecsc_
#define	comm_init_buffer_data	comm_init_buffer_data_
#define	comm_init_send_buf	comm_init_send_buf_
#define	comm_save_recv_c	comm_save_recv_c_
#define	comm_save_recv_nonblocking_buf	comm_save_recv_nonblocking_buf_
#define	mpix_pkshort_scalar	mpix_pkshort_scalar_
#define	mpix_upkshort_scalar	mpix_upkshort_scalar_
#define	mpix_pkshort_vec		mpix_pkshort_vec_
#define	mpix_upkshort_vec	mpix_upkshort_vec_
#define	mpix_pkstring		mpix_pkstring_
#define	mpix_upkstring		mpix_upkstring_
#else
# error 'F77_EXT_NAMES: no such case'
#endif
