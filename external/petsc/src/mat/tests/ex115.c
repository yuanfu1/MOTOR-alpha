static char help[] = "Tests MatHYPRE\n";

#include <petscmathypre.h>

int main(int argc, char **args)
{
  Mat                 A, B, C, D;
  Mat                 pAB, CD;
  hypre_ParCSRMatrix *parcsr;
  PetscReal           err;
  PetscInt            i, j, N = 6, M = 6;
  PetscBool           flg, testptap = PETSC_TRUE, testmatmatmult = PETSC_TRUE;
  PetscReal           norm;
  char                file[256];
  MatType             mtype = MATAIJ;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, NULL, help));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-f", file, sizeof(file), &flg));
#if defined(PETSC_USE_COMPLEX)
  testptap       = PETSC_FALSE;
  testmatmatmult = PETSC_FALSE;
  PetscCall(PetscOptionsInsertString(NULL, "-options_left 0"));
#endif
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-ptap", &testptap, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-matmatmult", &testmatmatmult, NULL));
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
#if PetscDefined(HAVE_HYPRE_DEVICE)
  #if PetscDefined(HAVE_HIP)
  mtype = MATAIJHIPSPARSE;
  #elif PetscDefined(HAVE_CUDA)
  mtype = MATAIJCUSPARSE;
  #endif
#endif

  if (!flg) { /* Create a matrix and test MatSetValues */
    PetscMPIInt size;

    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-M", &M, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-N", &N, NULL));
    PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, M, N));
    PetscCall(MatSetType(A, mtype));
    PetscCall(MatSeqAIJSetPreallocation(A, 9, NULL));
    PetscCall(MatMPIAIJSetPreallocation(A, 9, NULL, 9, NULL));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &B));
    PetscCall(MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, M, N));
#if PetscDefined(HAVE_HYPRE_DEVICE)
    PetscCall(MatSetType(B, mtype));
#else
    PetscCall(MatSetType(B, MATHYPRE));
#endif
    PetscCall(MatSeqAIJSetPreallocation(B, 9, NULL));
    PetscCall(MatMPIAIJSetPreallocation(B, 9, NULL, 9, NULL));
    if (M == N) {
      PetscCall(MatHYPRESetPreallocation(B, 9, NULL, 9, NULL));
    } else {
      PetscCall(MatHYPRESetPreallocation(B, 6, NULL, 6, NULL));
    }
    if (M == N) {
      for (i = 0; i < M; i++) {
        PetscInt    cols[] = {0, 1, 2, 3, 4, 5};
        PetscScalar vals[] = {0, 1. / size, 2. / size, 3. / size, 4. / size, 5. / size};
        for (j = i - 2; j < i + 1; j++) {
          if (j >= N) {
            PetscCall(MatSetValue(A, i, N - 1, (1. * j * N + i) / (3. * N * size), ADD_VALUES));
            PetscCall(MatSetValue(B, i, N - 1, (1. * j * N + i) / (3. * N * size), ADD_VALUES));
          } else if (i > j) {
            PetscCall(MatSetValue(A, i, PetscMin(j, N - 1), (1. * j * N + i) / (2. * N * size), ADD_VALUES));
            PetscCall(MatSetValue(B, i, PetscMin(j, N - 1), (1. * j * N + i) / (2. * N * size), ADD_VALUES));
          } else {
            PetscCall(MatSetValue(A, i, PetscMin(j, N - 1), -1. - (1. * j * N + i) / (4. * N * size), ADD_VALUES));
            PetscCall(MatSetValue(B, i, PetscMin(j, N - 1), -1. - (1. * j * N + i) / (4. * N * size), ADD_VALUES));
          }
        }
        PetscCall(MatSetValues(A, 1, &i, PetscMin(6, N), cols, vals, ADD_VALUES));
        PetscCall(MatSetValues(B, 1, &i, PetscMin(6, N), cols, vals, ADD_VALUES));
      }
    } else {
      PetscInt  rows[2];
      PetscBool test_offproc = PETSC_FALSE;

      PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_offproc", &test_offproc, NULL));
      if (test_offproc) {
        const PetscInt *ranges;
        PetscMPIInt     rank;

        PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
        PetscCall(MatGetOwnershipRanges(A, &ranges));
        rows[0] = ranges[(rank + 1) % size];
        rows[1] = ranges[(rank + 1) % size + 1];
      } else {
        PetscCall(MatGetOwnershipRange(A, &rows[0], &rows[1]));
      }
      for (i = rows[0]; i < rows[1]; i++) {
        PetscInt    cols[] = {0, 1, 2, 3, 4, 5};
        PetscScalar vals[] = {-1, 1, -2, 2, -3, 3};

        PetscCall(MatSetValues(A, 1, &i, PetscMin(6, N), cols, vals, INSERT_VALUES));
        PetscCall(MatSetValues(B, 1, &i, PetscMin(6, N), cols, vals, INSERT_VALUES));
      }
    }
    /* MAT_FLUSH_ASSEMBLY currently not supported */
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
#if PetscDefined(HAVE_HYPRE_DEVICE)
    PetscCall(MatConvert(B, MATHYPRE, MAT_INPLACE_MATRIX, &B));
#endif

#if defined(PETSC_USE_COMPLEX)
    /* make the matrix imaginary */
    PetscCall(MatScale(A, PETSC_i));
    PetscCall(MatScale(B, PETSC_i));
#endif

#if !PetscDefined(HAVE_HYPRE_DEVICE)
    /* MatAXPY further exercises MatSetValues_HYPRE */
    PetscCall(MatAXPY(B, -1., A, DIFFERENT_NONZERO_PATTERN));
    PetscCall(MatConvert(B, MATMPIAIJ, MAT_INITIAL_MATRIX, &C));
    PetscCall(MatNorm(C, NORM_INFINITY, &err));
    PetscCheck(err <= PETSC_SMALL, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatSetValues %g", err);
    PetscCall(MatDestroy(&C));
#endif
    PetscCall(MatDestroy(&B));
  } else {
    PetscViewer viewer;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatLoad(A, viewer));
    PetscCall(MatSetType(A, mtype));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(MatGetSize(A, &M, &N));
  }

  /* check conversion routines */
  PetscCall(MatViewFromOptions(A, NULL, "-view_A"));
  PetscCall(MatConvert(A, MATHYPRE, MAT_INITIAL_MATRIX, &B));
  PetscCall(MatViewFromOptions(B, NULL, "-view_convert"));
  PetscCall(MatMultEqual(B, A, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error Mat HYPRE init");
  PetscCall(MatConvert(A, MATHYPRE, MAT_REUSE_MATRIX, &B));
  PetscCall(MatViewFromOptions(B, NULL, "-view_convert"));
  PetscCall(MatMultEqual(B, A, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error Mat HYPRE reuse");
  PetscCall(MatConvert(B, MATIS, MAT_INITIAL_MATRIX, &D));
  PetscCall(MatConvert(B, MATIS, MAT_REUSE_MATRIX, &D));
  PetscCall(MatMultEqual(D, A, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error Mat IS");
  PetscCall(MatConvert(B, MATAIJ, MAT_INITIAL_MATRIX, &C));
  PetscCall(MatConvert(B, MATAIJ, MAT_REUSE_MATRIX, &C));
  PetscCall(MatMultEqual(C, A, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error Mat AIJ");
  PetscCall(MatAXPY(C, -1., A, SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(C, NORM_INFINITY, &err));
  PetscCheck(err <= PETSC_SMALL, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error Mat AIJ %g", err);
  PetscCall(MatDestroy(&C));
  PetscCall(MatConvert(D, MATAIJ, MAT_INITIAL_MATRIX, &C));
  PetscCall(MatAXPY(C, -1., A, SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(C, NORM_INFINITY, &err));
  PetscCheck(err <= PETSC_SMALL, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error Mat IS %g", err);
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&D));

  /* check MatCreateFromParCSR */
  PetscCall(MatHYPREGetParCSR(B, &parcsr));
  PetscCall(MatCreateFromParCSR(parcsr, MATAIJ, PETSC_COPY_VALUES, &D));
  PetscCall(MatDestroy(&D));
  PetscCall(MatCreateFromParCSR(parcsr, MATHYPRE, PETSC_USE_POINTER, &C));

  /* check MatMult operations */
  PetscCall(MatMultEqual(A, B, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMult B");
  PetscCall(MatMultEqual(A, C, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMult C");
  PetscCall(MatMultAddEqual(A, B, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMultAdd B");
  PetscCall(MatMultAddEqual(A, C, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMultAdd C");
  PetscCall(MatMultTransposeEqual(A, B, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMultTranspose B");
  PetscCall(MatMultTransposeEqual(A, C, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMultTranspose C");
  PetscCall(MatMultTransposeAddEqual(A, B, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMultTransposeAdd B");
  PetscCall(MatMultTransposeAddEqual(A, C, 4, &flg));
  PetscCheck(flg, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMultTransposeAdd C");

  /* check PtAP */
  if (testptap && M == N) {
    Mat pP, hP;

    /* PETSc MatPtAP -> output is a MatAIJ
       It uses HYPRE functions when -matptap_via hypre is specified at command line */
    PetscCall(MatPtAP(A, A, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &pP));
    PetscCall(MatPtAP(A, A, MAT_REUSE_MATRIX, PETSC_DETERMINE, &pP));
    PetscCall(MatNorm(pP, NORM_INFINITY, &norm));
    PetscCall(MatPtAPMultEqual(A, A, pP, 10, &flg));
    PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "Error in MatPtAP_MatAIJ");

    /* MatPtAP_HYPRE_HYPRE -> output is a MatHYPRE */
    PetscCall(MatPtAP(C, B, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &hP));
    PetscCall(MatPtAP(C, B, MAT_REUSE_MATRIX, PETSC_DETERMINE, &hP));
    PetscCall(MatPtAPMultEqual(C, B, hP, 10, &flg));
    PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "Error in MatPtAP_HYPRE_HYPRE");

    /* Test MatAXPY_Basic() */
#if !PetscDefined(HAVE_HYPRE_DEVICE)
    PetscCall(MatAXPY(hP, -1., pP, DIFFERENT_NONZERO_PATTERN));
    PetscCall(MatHasOperation(hP, MATOP_NORM, &flg));
    if (!flg) { /* TODO add MatNorm_HYPRE */
      PetscCall(MatConvert(hP, MATAIJ, MAT_INPLACE_MATRIX, &hP));
    }
    PetscCall(MatNorm(hP, NORM_INFINITY, &err));
    PetscCheck(err / norm <= PETSC_SMALL, PetscObjectComm((PetscObject)hP), PETSC_ERR_PLIB, "Error MatPtAP %g %g", err, norm);
#endif
    PetscCall(MatDestroy(&pP));
    PetscCall(MatDestroy(&hP));

    /* MatPtAP_AIJ_HYPRE -> output can be decided at runtime with -matptap_hypre_outtype */
#if !PetscDefined(HAVE_HYPRE_DEVICE)
    PetscCall(MatPtAP(A, B, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &hP));
    PetscCall(MatPtAP(A, B, MAT_REUSE_MATRIX, PETSC_DETERMINE, &hP));
    PetscCall(MatPtAPMultEqual(A, B, hP, 10, &flg));
    PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "Error in MatPtAP_AIJ_HYPRE");
    PetscCall(MatDestroy(&hP));
#endif
  }
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&B));

  /* check MatMatMult */
  if (testmatmatmult) {
    PetscCall(MatTranspose(A, MAT_INITIAL_MATRIX, &B));
    PetscCall(MatConvert(A, MATHYPRE, MAT_INITIAL_MATRIX, &C));
    PetscCall(MatConvert(B, MATHYPRE, MAT_INITIAL_MATRIX, &D));

    /* PETSc MatMatMult -> output is a MatAIJ
       It uses HYPRE functions when -matmatmult_via hypre is specified at command line */
    PetscCall(MatMatMult(A, B, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &pAB));
    PetscCall(MatMatMult(A, B, MAT_REUSE_MATRIX, PETSC_DETERMINE, &pAB));
    PetscCall(MatNorm(pAB, NORM_INFINITY, &norm));
    PetscCall(MatMatMultEqual(A, B, pAB, 10, &flg));
    PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "Error in MatMatMult_AIJ_AIJ");

    /* MatMatMult_HYPRE_HYPRE -> output is a MatHYPRE */
    PetscCall(MatMatMult(C, D, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &CD));
    PetscCall(MatMatMult(C, D, MAT_REUSE_MATRIX, PETSC_DETERMINE, &CD));
    PetscCall(MatMatMultEqual(C, D, CD, 10, &flg));
    PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_PLIB, "Error in MatMatMult_HYPRE_HYPRE");

    /* Test MatAXPY_Basic() */
#if !PetscDefined(HAVE_HYPRE_DEVICE)
    PetscCall(MatAXPY(CD, -1., pAB, DIFFERENT_NONZERO_PATTERN));

    PetscCall(MatHasOperation(CD, MATOP_NORM, &flg));
    if (!flg) { /* TODO add MatNorm_HYPRE */
      PetscCall(MatConvert(CD, MATAIJ, MAT_INPLACE_MATRIX, &CD));
    }
    PetscCall(MatNorm(CD, NORM_INFINITY, &err));
    PetscCheck((err / norm) <= PETSC_SMALL, PetscObjectComm((PetscObject)CD), PETSC_ERR_PLIB, "Error MatMatMult %g %g", err, norm);
#endif

    PetscCall(MatDestroy(&C));
    PetscCall(MatDestroy(&D));
    PetscCall(MatDestroy(&pAB));
    PetscCall(MatDestroy(&CD));

    /* When configured with HYPRE, MatMatMatMult is available for the triplet transpose(aij)-aij-aij */
#if !PetscDefined(HAVE_HYPRE_DEVICE)
    Mat CAB;
    PetscCall(MatCreateTranspose(A, &C));
    PetscCall(MatMatMatMult(C, A, B, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &CAB));
    PetscCall(MatDestroy(&C));
    PetscCall(MatTranspose(A, MAT_INITIAL_MATRIX, &C));
    PetscCall(MatMatMult(C, A, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &D));
    PetscCall(MatDestroy(&C));
    PetscCall(MatMatMult(D, B, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &C));
    PetscCall(MatNorm(C, NORM_INFINITY, &norm));
    PetscCall(MatAXPY(C, -1., CAB, DIFFERENT_NONZERO_PATTERN));
    PetscCall(MatNorm(C, NORM_INFINITY, &err));
    PetscCheck((err / norm) <= PETSC_SMALL, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatMatMatMult %g %g", err, norm);
    PetscCall(MatDestroy(&C));
    PetscCall(MatDestroy(&D));
    PetscCall(MatDestroy(&CAB));
#endif
    PetscCall(MatDestroy(&B));
  }

  /* Check MatView */
  PetscCall(MatConvert(A, MATHYPRE, MAT_INITIAL_MATRIX, &B));
  PetscCall(MatViewFromOptions(B, NULL, "-view_B"));

  /* Check MatDuplicate/MatCopy */
  for (j = 0; j < 3; j++) {
    MatDuplicateOption dop;

    dop = MAT_COPY_VALUES;
    if (j == 1) dop = MAT_DO_NOT_COPY_VALUES;
    if (j == 2) dop = MAT_SHARE_NONZERO_PATTERN;

    for (i = 0; i < 3; i++) {
      MatStructure str;

      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Dup/Copy tests: %" PetscInt_FMT " %" PetscInt_FMT "\n", j, i));

      str = DIFFERENT_NONZERO_PATTERN;
      if (i == 1) str = SAME_NONZERO_PATTERN;
      if (i == 2) str = SUBSET_NONZERO_PATTERN;

      PetscCall(MatDuplicate(A, dop, &C));
      PetscCall(MatDuplicate(B, dop, &D));
      if (dop != MAT_COPY_VALUES) {
        PetscCall(MatCopy(A, C, str));
        PetscCall(MatCopy(B, D, str));
      }
      /* AXPY with AIJ and HYPRE */
#if !PetscDefined(HAVE_HYPRE_DEVICE)
      PetscCall(MatAXPY(C, -1.0, D, str));
      PetscCall(MatNorm(C, NORM_INFINITY, &err));
#else
      err = 0.0;
#endif
      if (err > PETSC_SMALL) {
        PetscCall(MatViewFromOptions(A, NULL, "-view_duplicate_diff"));
        PetscCall(MatViewFromOptions(B, NULL, "-view_duplicate_diff"));
        PetscCall(MatViewFromOptions(C, NULL, "-view_duplicate_diff"));
        PetscCall(MatViewFromOptions(D, NULL, "-view_duplicate_diff"));
        SETERRQ(PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error test 1 MatDuplicate/MatCopy %g (%" PetscInt_FMT ",%" PetscInt_FMT ")", err, j, i);
      }
      /* AXPY with HYPRE and HYPRE */
      PetscCall(MatAXPY(D, -1.0, B, str));
      if (err > PETSC_SMALL) {
        PetscCall(MatViewFromOptions(A, NULL, "-view_duplicate_diff"));
        PetscCall(MatViewFromOptions(B, NULL, "-view_duplicate_diff"));
        PetscCall(MatViewFromOptions(D, NULL, "-view_duplicate_diff"));
        SETERRQ(PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error test 2 MatDuplicate/MatCopy %g (%" PetscInt_FMT ",%" PetscInt_FMT ")", err, j, i);
      }
#if !PetscDefined(HAVE_HYPRE_DEVICE)
      /* Copy from HYPRE to AIJ */
      PetscCall(MatCopy(B, C, str));
      /* Copy from AIJ to HYPRE */
      PetscCall(MatCopy(A, D, str));
      /* AXPY with HYPRE and AIJ */
      PetscCall(MatAXPY(D, -1.0, C, str));
      PetscCall(MatHasOperation(D, MATOP_NORM, &flg));
      if (!flg) { /* TODO add MatNorm_HYPRE */
        PetscCall(MatConvert(D, MATAIJ, MAT_INPLACE_MATRIX, &D));
      }
      PetscCall(MatNorm(D, NORM_INFINITY, &err));
      if (err > PETSC_SMALL) {
        PetscCall(MatViewFromOptions(A, NULL, "-view_duplicate_diff"));
        PetscCall(MatViewFromOptions(C, NULL, "-view_duplicate_diff"));
        PetscCall(MatViewFromOptions(D, NULL, "-view_duplicate_diff"));
        SETERRQ(PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error test 3 MatDuplicate/MatCopy %g (%" PetscInt_FMT ",%" PetscInt_FMT ")", err, j, i);
      }
#endif
      PetscCall(MatDestroy(&C));
      PetscCall(MatDestroy(&D));
    }
  }
  PetscCall(MatDestroy(&B));

  PetscCall(MatHasCongruentLayouts(A, &flg));
  if (flg) {
    Vec y, y2;

    PetscCall(MatConvert(A, MATHYPRE, MAT_INITIAL_MATRIX, &B));
    PetscCall(MatCreateVecs(A, NULL, &y));
    PetscCall(MatCreateVecs(B, NULL, &y2));
    PetscCall(MatGetDiagonal(A, y));
    PetscCall(MatGetDiagonal(B, y2));
    PetscCall(VecAXPY(y2, -1.0, y));
    PetscCall(VecNorm(y2, NORM_INFINITY, &err));
    if (err > PETSC_SMALL) {
      PetscCall(VecViewFromOptions(y, NULL, "-view_diagonal_diff"));
      PetscCall(VecViewFromOptions(y2, NULL, "-view_diagonal_diff"));
      SETERRQ(PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Error MatGetDiagonal %g", err);
    }
    PetscCall(MatDestroy(&B));
    PetscCall(VecDestroy(&y));
    PetscCall(VecDestroy(&y2));
  }

  PetscCall(MatDestroy(&A));

  PetscCall(PetscFinalize());
  return 0;
}

/*TEST

   build:
      requires: hypre

   test:
      suffix: 1
      args: -N 11 -M 11
      output_file: output/ex115_1.out

   test:
      suffix: 2
      nsize: 3
      args: -N 13 -M 13 -matmatmult_via hypre -options_left 0
      output_file: output/ex115_1.out

   test:
      suffix: 3
      nsize: 4
      args: -M 13 -N 7 -matmatmult_via hypre -options_left 0
      output_file: output/ex115_1.out

   test:
      suffix: 4
      nsize: 2
      args: -M 12 -N 19
      output_file: output/ex115_1.out

   test:
      suffix: 5
      nsize: 3
      args: -M 13 -N 13 -options_left 0 -matptap_via hypre -matptap_hypre_outtype hypre
      output_file: output/ex115_1.out

   test:
      requires: !defined(PETSC_HAVE_HYPRE_DEVICE)
      suffix: 6
      nsize: 3
      args: -M 12 -N 19 -test_offproc
      output_file: output/ex115_1.out

   test:
      requires: !defined(PETSC_HAVE_HYPRE_DEVICE)
      suffix: 7
      nsize: 3
      args: -M 19 -N 12 -test_offproc -view_B ::ascii_info_detail
      output_file: output/ex115_7.out

   test:
      requires: !defined(PETSC_HAVE_HYPRE_DEVICE)
      suffix: 8
      nsize: 3
      args: -M 1 -N 12 -test_offproc
      output_file: output/ex115_1.out

   test:
      requires: !defined(PETSC_HAVE_HYPRE_DEVICE)
      suffix: 9
      nsize: 3
      args: -M 1 -N 2 -test_offproc
      output_file: output/ex115_1.out

TEST*/
