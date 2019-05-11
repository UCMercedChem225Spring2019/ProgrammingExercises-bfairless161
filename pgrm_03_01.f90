      program pgrm_03_01
!
!     This program computes the inverse square-root of a matrix.
!
!
!     At execution time, the program expects 2 command line arguments: (1) nDim;
!     and (2) the matrix being raised to the (1/2) power.
!

      implicit none
      integer,parameter::unitIn=10
      integer::i,iError,nDim,lenSym
      real,dimension(:),allocatable::inputSymMatrix
      real,dimension(:,:),allocatable::inputSqMatrix,invSqrtInputMatrix
      character(len=256)::cmdlineArg,Filename
!
!     Begin by reading the leading dimension of the matrix and the input file
!     name from the command line. Then, open the file and read the input matrix,
!     inputSymMatrix
!


      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nDim
      lenSym = (nDim*(nDim+1))/2
      allocate(inputSymMatrix(lenSym),inputSqMatrix(nDim,nDim),  &
        invSqrtInputMatrix(nDim,nDim))
      call Get_Command_Argument(2,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) inputSymMatrix(i)
      endDo
      close(Unit=unitIn)
!
!     Form the square-root of inputSymMatrix. The result is loaded into square
!     (full storage) matrix invSqrtInputMatrix.
!

      write(*,*)' The matrix loaded (column) upper-triangle packed:'
      call SymmetricPacked2Matrix_UpperPac(nDim,inputSymMatrix,  &
        inputSqMatrix)
      write(*,*)' Input Matrix:'
      call Print_Matrix_Full_Real(inputSqMatrix,nDim,nDim)
      call InvSQRT_SymMatrix(nDim,inputSymMatrix,invSqrtInputMatrix)
      write(*,*)' Inverse SQRT Matrix:'
      call Print_Matrix_Full_Real(invSqrtInputMatrix,nDim,nDim)
      write(*,*)' Matrix product that should be the identity.'
      call Print_Matrix_Full_Real(MatMul(MatMul(invSqrtInputMatrix,  &
        invSqrtInputMatrix),inputSqMatrix),nDim,nDim)
!
      end program pgrm_03_01

      Subroutine Print_Matrix_Full_Real(AMat,M,N)
!
!     This subroutine prints a real matrix that is fully dimension -
!     i.e.,
!     not stored in packed form. AMat is the matrix, which is
!     dimensioned
!     (M,N).
!
!     The output of this routine is sent to unit number 6 (set by the
!     local
!     parameter integer IOut).
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::M,N
      real,dimension(M,N),intent(in)::AMat
!
!     Local variables
      integer,parameter::IOut=6,NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      Return
      End Subroutine Print_Matrix_Full_Real

        Subroutine SymmetricPacked2Matrix_UpperPac(N,ArrayIn,AMatOut)
       
        Implicit None
        Real,Dimension(N*N)::ArrayIn
        Real,Dimension(N,N)::AMatOut
        Integer::i,j,k,N
        k=1
        Do j=1, N
        Do i=1, j
        AMatOut(i,j)=ArrayIn(k)
        AMatOut(j,i)=AMatOut(i,j)
        k=k+1
        EndDo
        EndDo

! Do i=1, N
! Do j=i, N
! AMatOut(i,j)=ArrayIn(j*(j-1)/2+i)
! EndDo
! EndDo
! Do i=1, N
! Do j=i, N
! AMatOut(i,j)=AMatOut(j,i)
! EndDo
! EndDo
!
Return
        End Subroutine SymmetricPacked2Matrix_UpperPac

        Subroutine InvSQRT_SymMatrix(nDim,inputSymMatrix,invSqrtInputMatrix)
      
      IMPLICIT NONE
      Integer, INTENT(IN)  :: nDim
      Real,Dimension((nDim*(nDim+1))/2),INTENT(IN) :: inputSymMatrix
      Real,Dimension(nDim,nDim),  INTENT(OUT) :: invSqrtInputMatrix

      Integer:: i,j,k,IError,s
      Real,Dimension(:,:), Allocatable::EVecs, Temp_Matrix,EValMatx
      Real,Dimension(:), Allocatable::EVals, Temp_Vector,TempDump
      Allocate (EVals(nDim), EVecs(nDim,nDim),
        Temp_Vector(3*nDim),EValMatx(nDim,nDim))
      Allocate (Temp_Matrix(nDim,nDim),TempDump((nDim*(nDim+1))/2))
      Do i=1,(nDim*(nDim+1))/2
        TempDump(i) = inputSymMatrix(i)
      Enddo
      Call SSPEV("V", "U", nDim, TempDump, EVals, EVecs,nDim, & 
         Temp_Vector,IError)

        Do i=1,nDim
          EValMatx(i,i) = 1/(sqrt(EVals(i)))
        Enddo
 
        Mat = Matmul(Matmul(Evecs,Evalmay),transpose(Evecs))
        RETURN

      END SUBROUTINE InvSQRT_SymMatrix

