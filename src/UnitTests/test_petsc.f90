program main
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
  PetscErrorCode ierr
  Mat A
  Vec b
  Vec x
  KSP ksp !求解方法
  PC pc   !矩阵预处理方法
  PetscInt :: ione=1 !长整型的1，设置矩阵和向量元素时需要设置每次插入值的个数
  PetscInt :: m=10 !矩阵规模，本案例，矩阵规模为m行，m列
  PetscInt i, j !矩阵的行和列
  PetscInt rstart,rend !本进程，行的范围
  PetscScalar  v !临时变量，用于设置矩阵元素
  PetscReal :: tol=1.0e-7 !求解误差
  PetscReal :: rdefault=-2.0 !对应C++版Petsc中的PETSC_DEFAULT
  PetscInt :: idefault=-2 !对应C++版Petsc中的PETSC_DEFAULT
  
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRA(ierr)

    ! 创建矩阵和向量
  call MatCreate(PETSC_COMM_WORLD, A, ierr);CHKERRA(ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, m, ierr);CHKERRA(ierr)
  call MatSetType(A, MATMPIAIJ, ierr);CHKERRA(ierr)
  call MatSetUp(A, ierr)
  
  call MatGetOwnershipRange(A, rstart, rend, ierr);CHKERRA(ierr)
  
  ! 设置系数矩阵中元素的数值
  do i = rstart, rend-1
	v=1.0
    call MatSetValues(A, ione,[i],ione, [i], [v], INSERT_VALUES, ierr);CHKERRA(ierr)
    j=i-1;
    if (i > 0) call MatSetValues(A, ione,[i],ione, [j], [v], INSERT_VALUES, ierr);CHKERRA(ierr)
	j=i+1
    if (i < m-1) call MatSetValues(A,ione, [i], ione,[j], [v], INSERT_VALUES, ierr);CHKERRA(ierr)
  end do
  
  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRA(ierr)
  !call MatView(A, PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr) !显示系数矩阵

  ! 创建右端项向量b和未知向量x
  call VecCreate(PETSC_COMM_WORLD, b, ierr);CHKERRA(ierr)
  call VecSetSizes(b, PETSC_DECIDE, m,ierr);CHKERRA(ierr)
  call VecSetFromOptions(b,ierr);CHKERRA(ierr)
  !call VecGetOwnershipRange(b, rstart, rend,ierr);CHKERRA(ierr)
  
  !设置向量b的元素值
  do i = rstart, rend-1
	v=i
    call VecSetValues(b,ione,[i],[v],INSERT_VALUES,ierr);CHKERRA(ierr)
  end do

  call VecAssemblyBegin(b,ierr);CHKERRA(ierr)
  call VecAssemblyEnd(b,ierr);CHKERRA(ierr)  
  
  !call VecView(b, PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr) !显示b向量
  
  call VecDuplicate(b, x, ierr);CHKERRA(ierr)

  ! 配置KSP求解器
  call KSPCreate(PETSC_COMM_WORLD, ksp, ierr);CHKERRA(ierr)
  call KSPSetOperators(ksp, A, A, ierr);CHKERRA(ierr)
  call KSPGetPC(ksp, pc, ierr);CHKERRA(ierr)
  call PCSetType(pc, PCJACOBI, ierr);CHKERRA(ierr)  ! 设置雅可比预处理器
  call KSPSetFromOptions(ksp, ierr);CHKERRA(ierr)
  call KSPSetType(ksp, KSPCG, ierr);CHKERRA(ierr)
  

  call KSPSetTolerances(ksp,tol,rdefault,rdefault,idefault,ierr);CHKERRA(ierr)

  ! 求解
  call KSPSolve(ksp, b, x, ierr);CHKERRA(ierr)
  
  !输出计算结果
  call VecView(x, PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
  
  call KSPDestroy(ksp, ierr);CHKERRA(ierr)
  call MatDestroy(A, ierr);CHKERRA(ierr)
  call VecDestroy(x, ierr);CHKERRA(ierr)
  call VecDestroy(b, ierr);CHKERRA(ierr)
  
  call PetscFinalize(ierr);CHKERRA(ierr)
end program
