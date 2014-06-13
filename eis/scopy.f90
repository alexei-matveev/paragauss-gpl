        subroutine scopy(n,const,nd1,a,nd2)
        use type_module
        implicit real(kind=r8_kind) (A-H,O-Z)
        real(kind=r8_kind) A(n)
        do i=1,n
        a(i)=const
        end do
        end
