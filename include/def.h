
#if defined __STDC__ && __STDC__
# define __CONCAT(x,y) x ## y
# define __STRING(x) #x
#else
# define __CONCAT(x,y) x/**/y
# define __STRING(x) "x"
#endif

#define WARN(X)      print *, "WARNING: ", X ," at ", __FILE__ , __LINE__
#define ABORT(X)     print *,"ABORT: at ",__FILE__,__LINE__;call error_handler(X)
#define TRACE(X) !!$FPP$ call trace_hook(X, __FILE__, __LINE__)

#if FPP_NODEBUG
# define ASSERT(expr) !!$FPP$ check for expr here!
#else
# define ASSERT(expr) if(.not.(expr))call ASSERT_FAILED(__STRING(expr),__FILE__,__LINE__)
# define FPP_ASSERT
#endif

#ifdef USE_MPI_MODULE
# define USE_MPI use mpi /* some build systems provide MPI module */
#else
# define USE_MPI use commparameter_module /* our wrapper around mpif.h */
#endif

#ifdef FPP_DEBUG
# define _DPRINT 1
# define DPRINT print *,   "dbg:",
# define DWRITE(X,Y) write(X,Y)
# define DPRINTF(X) print X,
# define DLPRINT(L) if( L .LE. FPP_DEBUG ) print *,"dbg(",L,"):",
# define USE_DEBUG use debug
# define DCALL call
# define DTRACE(x) print*,x,"dbg:TRACE: passed ",__FILE__,":",__LINE__
# define DSHOW(x)  print*,"dbg:",__STRING(x=),x
#else
# undef _DPRINT
# define DPRINT !!$FPP$ print *,
# define DWRITE !!$FPP$ write(*,*)
# define DPRINTF(X) !!$FPP$ print X,
# define DLPRINT !!$FPP$
# define USE_DEBUG !!$FPP$ dont use debug
# define DCALL !!$FPP$ dont call
# define DTRACE(x) !!$FPP$ dont print*,'dbg: TRACE: passed ',__FILE__,':',__LINE__
# define DSHOW(x) !!$FPP$ dont print*,__STRING(x=),x
#endif

! define ARRSTAT(X) sum(X),sum(abs(X)),maxval(X),minval(X)
# define ARRSTAT(X) sum(X),sum(abs(X))

#ifdef WITH_MEMLOG
# define USE_MEMLOG use memlog_module
# define MEMLOG(a) call meminc(a)
# define MEMSET(a) call memset(a)
# define MEMUSAGE(a) memusage()
#else
# define USE_MEMLOG !!$dont use memlog_module
# define MEMLOG(a) !!$dont call meminc(a)
# define MEMSET(a) !!$dont call memset(a)
# define MEMUSAGE(a) 0
#endif

#if FPP_FAST_COMPILE
# define FPP_NOMDA
#endif

! Machine specific settings:
#if _VPP
# define _VECTOR
#endif

#if FPP_TIMERS
#    define FPP_TIMER_VARS(t)   __CONCAT(t,c_),__CONCAT(t,o_),__CONCAT(t,__)
#    define FPP_TIMER_ZERO(t)   __CONCAT(t,c_)=0; __CONCAT(t,o_)=0; __CONCAT(t,__)=0
#    define FPP_USR_TIMER_DECL(t)   real::__CONCAT(t,c_),__CONCAT(t,o_),__CONCAT(t,__)=0.0
#    define FPP_USR_TIMER_VALUE(t)  __CONCAT(t,__)
#    define FPP_USR_TIMER_SLICE(t)  (__CONCAT(t,c_)-__CONCAT(t,o_))
#    define FPP_USR_TIMER_START(t)  call CPU_TIME(__CONCAT(t,o_))
#    define FPP_USR_TIMER_STOP(t)   call CPU_TIME(__CONCAT(t,c_));\
__CONCAT(t,__)=__CONCAT(t,__)+FPP_TIMER_SLICE(t)

#    define FPP_CLK_TIMER_DECL(t)   integer::__CONCAT(t,c_),__CONCAT(t,r_),__CONCAT(t,o_),__CONCAT(t,m_);\
real::__CONCAT(t,__)=0.0
#    define FPP_CLK_TIMER_VALUE(t)  (__CONCAT(t,__)/__CONCAT(t,r_))
#    define FPP_CLK_TIMER_SLICE(t)  ((real(__CONCAT(t,c_))-__CONCAT(t,o_))/__CONCAT(t,r_))
#    define FPP_CLK_TIMER_START(t)  call SYSTEM_CLOCK(__CONCAT(t,o_))
#    define FPP_CLK_TIMER_STOP(t)   call SYSTEM_CLOCK(__CONCAT(t,c_),__CONCAT(t,r_),__CONCAT(t,m_));\
if(__CONCAT(t,c_)<__CONCAT(t,o_))__CONCAT(t,o_)=__CONCAT(t,o_)-__CONCAT(t,m_);\
__CONCAT(t,__)=__CONCAT(t,__)+__CONCAT(t,c_)-__CONCAT(t,o_)

#    define FPP_MPI_TIMER_DECL(t)   double precision::__CONCAT(t,c_),__CONCAT(t,o_),__CONCAT(t,__)=0.0
#    define FPP_MPI_TIMER_VALUE(t)  __CONCAT(t,__)
#    define FPP_MPI_TIMER_SLICE(t)  (__CONCAT(t,c_)-__CONCAT(t,o_))
#    define FPP_MPI_TIMER_START(t)  __CONCAT(t,o_)=MPI_WTIME()
#    define FPP_MPI_TIMER_STOP(t)   __CONCAT(t,c_)=MPI_WTIME();__CONCAT(t,__)=__CONCAT(t,__)+FPP_TIMER_SLICE(t)

#  if FPP_TIMERS == 1
#    define FPP_TIMER_DECL(t)  FPP_CLK_TIMER_DECL(t)
#    define FPP_TIMER_SLICE(t) FPP_CLK_TIMER_SLICE(t)
#    define FPP_TIMER_START(t) FPP_CLK_TIMER_START(t)
#    define FPP_TIMER_STOP(t)  FPP_CLK_TIMER_STOP(t)
#    define FPP_TIMER_VALUE(t) FPP_CLK_TIMER_VALUE(t)
#  endif

#  if FPP_TIMERS == 2
#    define FPP_TIMER_DECL(t)  FPP_USR_TIMER_DECL(t)
#    define FPP_TIMER_SLICE(t) FPP_USR_TIMER_SLICE(t)
#    define FPP_TIMER_START(t) FPP_USR_TIMER_START(t)
#    define FPP_TIMER_STOP(t)  FPP_USR_TIMER_STOP(t)
#    define FPP_TIMER_VALUE(t) FPP_USR_TIMER_VALUE(t)
#  endif

#  if FPP_TIMERS == 3
#    define FPP_TIMER_DECL(t)  FPP_MPI_TIMER_DECL(t)
#    define FPP_TIMER_SLICE(t) FPP_MPI_TIMER_SLICE(t)
#    define FPP_TIMER_START(t) FPP_MPI_TIMER_START(t)
#    define FPP_TIMER_STOP(t)  FPP_MPI_TIMER_STOP(t)
#    define FPP_TIMER_VALUE(t) FPP_MPI_TIMER_VALUE(t)
#  endif

#  define FPP_TIMER_PRINT(t)  DPRINT ' time=', FPP_TIMER_VALUE(t), '(' , __STRING(t),')'
#else
#  define FPP_TIMER_DECL(t)  !!$FPP$ decl: timer t
#  define FPP_TIMER_SLICE(t) 0
#  define FPP_TIMER_START(t) !!$FPP$ start: timer t
#  define FPP_TIMER_STOP(t)  !!$FPP$ stop: timer t
#  define FPP_TIMER_PRINT(t) !!$FPP$ print: timer t
#  define FPP_TIMER_VALUE(t) 0
#endif




