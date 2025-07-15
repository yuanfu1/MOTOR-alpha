#include <../src/vec/is/sf/impls/basic/sfbasic.h>
#include <../src/vec/is/sf/impls/basic/sfpack.h>
#include <petsc/private/viewerimpl.h>

// Init persistent MPI send/recv requests
static PetscErrorCode PetscSFLinkInitMPIRequests_Persistent_Basic(PetscSF sf, PetscSFLink link, PetscSFDirection direction)
{
  PetscSF_Basic     *bas = (PetscSF_Basic *)sf->data;
  PetscInt           cnt;
  PetscMPIInt        nrootranks, ndrootranks, nleafranks, ndleafranks;
  const PetscInt    *rootoffset, *leafoffset;
  MPI_Aint           disp;
  MPI_Comm           comm          = PetscObjectComm((PetscObject)sf);
  MPI_Datatype       unit          = link->unit;
  const PetscMemType rootmtype_mpi = link->rootmtype_mpi, leafmtype_mpi = link->leafmtype_mpi; /* Used to select buffers passed to MPI */
  const PetscInt     rootdirect_mpi = link->rootdirect_mpi, leafdirect_mpi = link->leafdirect_mpi;

  PetscFunctionBegin;
  if (bas->rootbuflen[PETSCSF_REMOTE] && !link->rootreqsinited[direction][rootmtype_mpi][rootdirect_mpi]) {
    PetscCall(PetscSFGetRootInfo_Basic(sf, &nrootranks, &ndrootranks, NULL, &rootoffset, NULL));
    if (direction == PETSCSF_LEAF2ROOT) {
      for (PetscMPIInt i = ndrootranks, j = 0; i < nrootranks; i++, j++) {
        disp = (rootoffset[i] - rootoffset[ndrootranks]) * link->unitbytes;
        cnt  = rootoffset[i + 1] - rootoffset[i];
        PetscCallMPI(MPIU_Recv_init(link->rootbuf[PETSCSF_REMOTE][rootmtype_mpi] + disp, cnt, unit, bas->iranks[i], link->tag, comm, link->rootreqs[direction][rootmtype_mpi][rootdirect_mpi] + j));
      }
    } else { /* PETSCSF_ROOT2LEAF */
      for (PetscMPIInt i = ndrootranks, j = 0; i < nrootranks; i++, j++) {
        disp = (rootoffset[i] - rootoffset[ndrootranks]) * link->unitbytes;
        cnt  = rootoffset[i + 1] - rootoffset[i];
        PetscCallMPI(MPIU_Send_init(link->rootbuf[PETSCSF_REMOTE][rootmtype_mpi] + disp, cnt, unit, bas->iranks[i], link->tag, comm, link->rootreqs[direction][rootmtype_mpi][rootdirect_mpi] + j));
      }
    }
    link->rootreqsinited[direction][rootmtype_mpi][rootdirect_mpi] = PETSC_TRUE;
  }

  if (sf->leafbuflen[PETSCSF_REMOTE] && !link->leafreqsinited[direction][leafmtype_mpi][leafdirect_mpi]) {
    PetscCall(PetscSFGetLeafInfo_Basic(sf, &nleafranks, &ndleafranks, NULL, &leafoffset, NULL, NULL));
    if (direction == PETSCSF_LEAF2ROOT) {
      for (PetscMPIInt i = ndleafranks, j = 0; i < nleafranks; i++, j++) {
        disp = (leafoffset[i] - leafoffset[ndleafranks]) * link->unitbytes;
        cnt  = leafoffset[i + 1] - leafoffset[i];
        PetscCallMPI(MPIU_Send_init(link->leafbuf[PETSCSF_REMOTE][leafmtype_mpi] + disp, cnt, unit, sf->ranks[i], link->tag, comm, link->leafreqs[direction][leafmtype_mpi][leafdirect_mpi] + j));
      }
    } else { /* PETSCSF_ROOT2LEAF */
      for (PetscMPIInt i = ndleafranks, j = 0; i < nleafranks; i++, j++) {
        disp = (leafoffset[i] - leafoffset[ndleafranks]) * link->unitbytes;
        cnt  = leafoffset[i + 1] - leafoffset[i];
        PetscCallMPI(MPIU_Recv_init(link->leafbuf[PETSCSF_REMOTE][leafmtype_mpi] + disp, cnt, unit, sf->ranks[i], link->tag, comm, link->leafreqs[direction][leafmtype_mpi][leafdirect_mpi] + j));
      }
    }
    link->leafreqsinited[direction][leafmtype_mpi][leafdirect_mpi] = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Start MPI requests. If use non-GPU aware MPI, we might need to copy data from device buf to host buf
static PetscErrorCode PetscSFLinkStartCommunication_Persistent_Basic(PetscSF sf, PetscSFLink link, PetscSFDirection direction)
{
  PetscMPIInt    nsreqs = 0, nrreqs = 0;
  MPI_Request   *sreqs = NULL, *rreqs = NULL;
  PetscSF_Basic *bas = (PetscSF_Basic *)sf->data;
  PetscInt       sbuflen, rbuflen;

  PetscFunctionBegin;
  rbuflen = (direction == PETSCSF_ROOT2LEAF) ? sf->leafbuflen[PETSCSF_REMOTE] : bas->rootbuflen[PETSCSF_REMOTE];
  if (rbuflen) {
    if (direction == PETSCSF_ROOT2LEAF) {
      nrreqs = sf->nleafreqs;
      PetscCall(PetscSFLinkGetMPIBuffersAndRequests(sf, link, direction, NULL, NULL, NULL, &rreqs));
    } else { /* leaf to root */
      nrreqs = bas->nrootreqs;
      PetscCall(PetscSFLinkGetMPIBuffersAndRequests(sf, link, direction, NULL, NULL, &rreqs, NULL));
    }
  }

  sbuflen = (direction == PETSCSF_ROOT2LEAF) ? bas->rootbuflen[PETSCSF_REMOTE] : sf->leafbuflen[PETSCSF_REMOTE];
  if (sbuflen) {
    if (direction == PETSCSF_ROOT2LEAF) {
      nsreqs = bas->nrootreqs;
      PetscCall(PetscSFLinkCopyRootBufferInCaseNotUseGpuAwareMPI(sf, link, PETSC_TRUE /*device2host before sending */));
      PetscCall(PetscSFLinkGetMPIBuffersAndRequests(sf, link, direction, NULL, NULL, &sreqs, NULL));
    } else { /* leaf to root */
      nsreqs = sf->nleafreqs;
      PetscCall(PetscSFLinkCopyLeafBufferInCaseNotUseGpuAwareMPI(sf, link, PETSC_TRUE));
      PetscCall(PetscSFLinkGetMPIBuffersAndRequests(sf, link, direction, NULL, NULL, NULL, &sreqs));
    }
  }
  PetscCall(PetscSFLinkSyncStreamBeforeCallMPI(sf, link)); // need to sync the stream to make BOTH sendbuf and recvbuf ready
  if (rbuflen) PetscCallMPI(MPI_Startall_irecv(rbuflen, link->unit, nrreqs, rreqs));
  if (sbuflen) PetscCallMPI(MPI_Startall_isend(sbuflen, link->unit, nsreqs, sreqs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_MPIX_STREAM)
// issue MPIX_Isend/Irecv_enqueue()
static PetscErrorCode PetscSFLinkStartCommunication_MPIX_Stream(PetscSF sf, PetscSFLink link, PetscSFDirection direction)
{
  PetscSF_Basic     *bas = (PetscSF_Basic *)sf->data;
  PetscInt           i, j;
  PetscMPIInt        nrootranks, ndrootranks, nleafranks, ndleafranks, cnt;
  const PetscInt    *rootoffset, *leafoffset;
  MPI_Aint           disp;
  MPI_Comm           stream_comm   = sf->stream_comm;
  MPI_Datatype       unit          = link->unit;
  const PetscMemType rootmtype_mpi = link->rootmtype_mpi, leafmtype_mpi = link->leafmtype_mpi; /* Used to select buffers passed to MPI */
  const PetscInt     rootdirect_mpi = link->rootdirect_mpi, leafdirect_mpi = link->leafdirect_mpi;

  PetscFunctionBegin;
  if (bas->rootbuflen[PETSCSF_REMOTE]) {
    PetscCall(PetscSFGetRootInfo_Basic(sf, &nrootranks, &ndrootranks, NULL, &rootoffset, NULL));
    if (direction == PETSCSF_LEAF2ROOT) {
      for (i = ndrootranks, j = 0; i < nrootranks; i++, j++) {
        disp = (rootoffset[i] - rootoffset[ndrootranks]) * link->unitbytes;
        cnt  = (PetscMPIInt)(rootoffset[i + 1] - rootoffset[i]);
        PetscCallMPI(MPIX_Irecv_enqueue(link->rootbuf[PETSCSF_REMOTE][rootmtype_mpi] + disp, cnt, unit, bas->iranks[i], link->tag, stream_comm, link->rootreqs[direction][rootmtype_mpi][rootdirect_mpi] + j));
      }
    } else { // PETSCSF_ROOT2LEAF
      for (i = ndrootranks, j = 0; i < nrootranks; i++, j++) {
        disp = (rootoffset[i] - rootoffset[ndrootranks]) * link->unitbytes;
        cnt  = (PetscMPIInt)(rootoffset[i + 1] - rootoffset[i]);
        // no need to sync the gpu stream!
        PetscCallMPI(MPIX_Isend_enqueue(link->rootbuf[PETSCSF_REMOTE][rootmtype_mpi] + disp, cnt, unit, bas->iranks[i], link->tag, stream_comm, link->rootreqs[direction][rootmtype_mpi][rootdirect_mpi] + j));
      }
    }
  }

  if (sf->leafbuflen[PETSCSF_REMOTE]) {
    PetscCall(PetscSFGetLeafInfo_Basic(sf, &nleafranks, &ndleafranks, NULL, &leafoffset, NULL, NULL));
    if (direction == PETSCSF_LEAF2ROOT) {
      for (i = ndleafranks, j = 0; i < nleafranks; i++, j++) {
        disp = (leafoffset[i] - leafoffset[ndleafranks]) * link->unitbytes;
        cnt  = (PetscMPIInt)(leafoffset[i + 1] - leafoffset[i]);
        // no need to sync the gpu stream!
        PetscCallMPI(MPIX_Isend_enqueue(link->leafbuf[PETSCSF_REMOTE][leafmtype_mpi] + disp, cnt, unit, sf->ranks[i], link->tag, stream_comm, link->leafreqs[direction][leafmtype_mpi][leafdirect_mpi] + j));
      }
    } else { // PETSCSF_ROOT2LEAF
      for (i = ndleafranks, j = 0; i < nleafranks; i++, j++) {
        disp = (leafoffset[i] - leafoffset[ndleafranks]) * link->unitbytes;
        cnt  = (PetscMPIInt)(leafoffset[i + 1] - leafoffset[i]);
        PetscCallMPI(MPIX_Irecv_enqueue(link->leafbuf[PETSCSF_REMOTE][leafmtype_mpi] + disp, cnt, unit, sf->ranks[i], link->tag, stream_comm, link->leafreqs[direction][leafmtype_mpi][leafdirect_mpi] + j));
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PetscSFLinkFinishCommunication_MPIX_Stream(PetscSF sf, PetscSFLink link, PetscSFDirection direction)
{
  PetscSF_Basic     *bas           = (PetscSF_Basic *)sf->data;
  const PetscMemType rootmtype_mpi = link->rootmtype_mpi, leafmtype_mpi = link->leafmtype_mpi;
  const PetscInt     rootdirect_mpi = link->rootdirect_mpi, leafdirect_mpi = link->leafdirect_mpi;

  PetscFunctionBegin;
  PetscCallMPI(MPIX_Waitall_enqueue(bas->nrootreqs, link->rootreqs[direction][rootmtype_mpi][rootdirect_mpi], MPI_STATUSES_IGNORE));
  PetscCallMPI(MPIX_Waitall_enqueue(sf->nleafreqs, link->leafreqs[direction][leafmtype_mpi][leafdirect_mpi], MPI_STATUSES_IGNORE));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode PetscSFSetCommunicationOps_Basic(PetscSF sf, PetscSFLink link)
{
  PetscFunctionBegin;
  link->InitMPIRequests    = PetscSFLinkInitMPIRequests_Persistent_Basic;
  link->StartCommunication = PetscSFLinkStartCommunication_Persistent_Basic;
#if defined(PETSC_HAVE_MPIX_STREAM)
  const PetscMemType rootmtype_mpi = link->rootmtype_mpi, leafmtype_mpi = link->leafmtype_mpi;
  if (sf->use_stream_aware_mpi && (PetscMemTypeDevice(rootmtype_mpi) || PetscMemTypeDevice(leafmtype_mpi))) {
    link->StartCommunication  = PetscSFLinkStartCommunication_MPIX_Stream;
    link->FinishCommunication = PetscSFLinkFinishCommunication_MPIX_Stream;
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*===================================================================================*/
/*              SF public interface implementations                                  */
/*===================================================================================*/
PETSC_INTERN PetscErrorCode PetscSFSetUp_Basic(PetscSF sf)
{
  PetscSF_Basic *bas = (PetscSF_Basic *)sf->data;
  PetscInt      *rlengths, *ilengths;
  PetscMPIInt    nRemoteRootRanks, nRemoteLeafRanks;
  PetscMPIInt    rank, niranks, *iranks, tag;
  MPI_Comm       comm;
  MPI_Group      group;
  MPI_Request   *rootreqs, *leafreqs;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_group(PETSC_COMM_SELF, &group));
  PetscCall(PetscSFSetUpRanks(sf, group));
  PetscCallMPI(MPI_Group_free(&group));
  PetscCall(PetscObjectGetComm((PetscObject)sf, &comm));
  PetscCall(PetscObjectGetNewTag((PetscObject)sf, &tag));
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  /*
   * Inform roots about how many leaves and from which ranks
   */
  PetscCall(PetscMalloc1(sf->nranks, &rlengths));
  /* Determine number, sending ranks and length of incoming */
  for (PetscMPIInt i = 0; i < sf->nranks; i++) { rlengths[i] = sf->roffset[i + 1] - sf->roffset[i]; /* Number of roots referenced by my leaves; for rank sf->ranks[i] */ }
  nRemoteRootRanks = sf->nranks - sf->ndranks;
  PetscCall(PetscCommBuildTwoSided(comm, 1, MPIU_INT, nRemoteRootRanks, PetscSafePointerPlusOffset(sf->ranks, sf->ndranks), PetscSafePointerPlusOffset(rlengths, sf->ndranks), &niranks, &iranks, (void **)&ilengths));

  /* Sort iranks. See use of VecScatterGetRemoteOrdered_Private() in MatGetBrowsOfAoCols_MPIAIJ() on why.
     We could sort ranks there at the price of allocating extra working arrays. Presumably, niranks is
     small and the sorting is cheap.
   */
  PetscCall(PetscSortMPIIntWithIntArray(niranks, iranks, ilengths));

  /* Partition into distinguished and non-distinguished incoming ranks */
  bas->ndiranks = sf->ndranks;
  bas->niranks  = bas->ndiranks + niranks;
  PetscCall(PetscMalloc2(bas->niranks, &bas->iranks, bas->niranks + 1, &bas->ioffset));
  bas->ioffset[0] = 0;
  for (PetscMPIInt i = 0; i < bas->ndiranks; i++) {
    bas->iranks[i]      = sf->ranks[i];
    bas->ioffset[i + 1] = bas->ioffset[i] + rlengths[i];
  }
  PetscCheck(bas->ndiranks <= 1 && (bas->ndiranks != 1 || bas->iranks[0] == rank), PETSC_COMM_SELF, PETSC_ERR_PLIB, "Broken setup for shared ranks");
  for (PetscMPIInt i = bas->ndiranks; i < bas->niranks; i++) {
    bas->iranks[i]      = iranks[i - bas->ndiranks];
    bas->ioffset[i + 1] = bas->ioffset[i] + ilengths[i - bas->ndiranks];
  }
  bas->itotal = bas->ioffset[bas->niranks];
  PetscCall(PetscFree(rlengths));
  PetscCall(PetscFree(iranks));
  PetscCall(PetscFree(ilengths));

  /* Send leaf identities to roots */
  nRemoteLeafRanks = bas->niranks - bas->ndiranks;
  PetscCall(PetscMalloc1(bas->itotal, &bas->irootloc));
  PetscCall(PetscMalloc2(nRemoteLeafRanks, &rootreqs, nRemoteRootRanks, &leafreqs));
  for (PetscMPIInt i = bas->ndiranks; i < bas->niranks; i++) PetscCallMPI(MPIU_Irecv(bas->irootloc + bas->ioffset[i], bas->ioffset[i + 1] - bas->ioffset[i], MPIU_INT, bas->iranks[i], tag, comm, &rootreqs[i - bas->ndiranks]));
  for (PetscMPIInt i = 0; i < sf->nranks; i++) {
    PetscInt npoints = sf->roffset[i + 1] - sf->roffset[i];
    if (i < sf->ndranks) {
      PetscCheck(sf->ranks[i] == rank, PETSC_COMM_SELF, PETSC_ERR_PLIB, "Cannot interpret distinguished leaf rank");
      PetscCheck(bas->iranks[0] == rank, PETSC_COMM_SELF, PETSC_ERR_PLIB, "Cannot interpret distinguished root rank");
      PetscCheck(npoints == bas->ioffset[1] - bas->ioffset[0], PETSC_COMM_SELF, PETSC_ERR_PLIB, "Distinguished rank exchange has mismatched lengths");
      PetscCall(PetscArraycpy(bas->irootloc + bas->ioffset[0], sf->rremote + sf->roffset[i], npoints));
      continue;
    }
    PetscCallMPI(MPIU_Isend(sf->rremote + sf->roffset[i], npoints, MPIU_INT, sf->ranks[i], tag, comm, &leafreqs[i - sf->ndranks]));
  }
  PetscCallMPI(MPI_Waitall(nRemoteLeafRanks, rootreqs, MPI_STATUSES_IGNORE));
  PetscCallMPI(MPI_Waitall(nRemoteRootRanks, leafreqs, MPI_STATUSES_IGNORE));

  sf->nleafreqs  = nRemoteRootRanks;
  bas->nrootreqs = nRemoteLeafRanks;

  /* Setup fields related to packing, such as rootbuflen[] */
  PetscCall(PetscSFSetUpPackFields(sf));
  PetscCall(PetscFree2(rootreqs, leafreqs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFReset_Basic(PetscSF sf)
{
  PetscSF_Basic *bas  = (PetscSF_Basic *)sf->data;
  PetscSFLink    link = bas->avail, next;

  PetscFunctionBegin;
  PetscCheck(!bas->inuse, PetscObjectComm((PetscObject)sf), PETSC_ERR_ARG_WRONGSTATE, "Outstanding operation has not been completed");
  PetscCall(PetscFree2(bas->iranks, bas->ioffset));
  PetscCall(PetscFree(bas->irootloc));

#if defined(PETSC_HAVE_DEVICE)
  for (int i = 0; i < 2; i++) PetscCall(PetscSFFree(sf, PETSC_MEMTYPE_DEVICE, bas->irootloc_d[i]));
#endif

#if defined(PETSC_HAVE_NVSHMEM)
  PetscCall(PetscSFReset_Basic_NVSHMEM(sf));
#endif

  for (; link; link = next) {
    next = link->next;
    PetscCall(PetscSFLinkDestroy(sf, link));
  }
  bas->avail = NULL;
  PetscCall(PetscSFResetPackFields(sf));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFDestroy_Basic(PetscSF sf)
{
  PetscFunctionBegin;
  PetscCall(PetscSFReset_Basic(sf));
  PetscCall(PetscFree(sf->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_SINGLE_LIBRARY)
  #include <petscmat.h>

PETSC_INTERN PetscErrorCode PetscSFView_Basic_PatternAndSizes(PetscSF sf, PetscViewer viewer)
{
  PetscSF_Basic     *bas = (PetscSF_Basic *)sf->data;
  PetscMPIInt        nrootranks, ndrootranks;
  const PetscInt    *rootoffset;
  PetscMPIInt        rank, size;
  const PetscMPIInt *rootranks;
  MPI_Comm           comm = PetscObjectComm((PetscObject)sf);
  PetscScalar        unitbytes;
  Mat                A;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  /* PetscSFView is most useful for the SF used in VecScatterBegin/End in MatMult etc, where we do
    PetscSFBcast, i.e., roots send data to leaves.  We dump the communication pattern into a matrix
    in senders' view point: how many bytes I will send to my neighbors.

    Looking at a column of the matrix, one can also know how many bytes the rank will receive from others.

    If PetscSFLink bas->inuse is available, we can use that to get tree vertex size. But that would give
    different interpretations for the same SF for different data types. Since we most care about VecScatter,
    we uniformly treat each vertex as a PetscScalar.
  */
  unitbytes = (PetscScalar)sizeof(PetscScalar);

  PetscCall(PetscSFGetRootInfo_Basic(sf, &nrootranks, &ndrootranks, &rootranks, &rootoffset, NULL));
  PetscCall(MatCreateAIJ(comm, 1, 1, size, size, 1, NULL, nrootranks - ndrootranks, NULL, &A));
  PetscCall(MatSetOptionsPrefix(A, "__petsc_internal__")); /* To prevent the internal A from taking any command line options */
  for (PetscMPIInt i = 0; i < nrootranks; i++) PetscCall(MatSetValue(A, rank, bas->iranks[i], (rootoffset[i + 1] - rootoffset[i]) * unitbytes, INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatView(A, viewer));
  PetscCall(MatDestroy(&A));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

PETSC_INTERN PetscErrorCode PetscSFView_Basic(PetscSF sf, PetscViewer viewer)
{
  PetscBool isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &isascii));
  if (isascii && viewer->format != PETSC_VIEWER_ASCII_MATLAB) PetscCall(PetscViewerASCIIPrintf(viewer, "  MultiSF sort=%s\n", sf->rankorder ? "rank-order" : "unordered"));
#if defined(PETSC_USE_SINGLE_LIBRARY)
  else {
    PetscBool isdraw, isbinary;
    PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERDRAW, &isdraw));
    PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERBINARY, &isbinary));
    if ((isascii && viewer->format == PETSC_VIEWER_ASCII_MATLAB) || isdraw || isbinary) PetscCall(PetscSFView_Basic_PatternAndSizes(sf, viewer));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFBcastBegin_Basic(PetscSF sf, MPI_Datatype unit, PetscMemType rootmtype, const void *rootdata, PetscMemType leafmtype, void *leafdata, MPI_Op op)
{
  PetscSFLink link = NULL;

  PetscFunctionBegin;
  /* Create a communication link, which provides buffers, MPI requests etc (if MPI is used) */
  PetscCall(PetscSFLinkCreate(sf, unit, rootmtype, rootdata, leafmtype, leafdata, op, PETSCSF_BCAST, &link));
  /* Pack rootdata to rootbuf for remote communication */
  PetscCall(PetscSFLinkPackRootData(sf, link, PETSCSF_REMOTE, rootdata));
  /* Start communication, e.g., post MPIU_Isend */
  PetscCall(PetscSFLinkStartCommunication(sf, link, PETSCSF_ROOT2LEAF));
  /* Do local scatter (i.e., self to self communication), which overlaps with the remote communication above */
  PetscCall(PetscSFLinkScatterLocal(sf, link, PETSCSF_ROOT2LEAF, (void *)rootdata, leafdata, op));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFBcastEnd_Basic(PetscSF sf, MPI_Datatype unit, const void *rootdata, void *leafdata, MPI_Op op)
{
  PetscSFLink link = NULL;

  PetscFunctionBegin;
  /* Retrieve the link used in XxxBegin() with root/leafdata as key */
  PetscCall(PetscSFLinkGetInUse(sf, unit, rootdata, leafdata, PETSC_OWN_POINTER, &link));
  /* Finish remote communication, e.g., post MPI_Waitall */
  PetscCall(PetscSFLinkFinishCommunication(sf, link, PETSCSF_ROOT2LEAF));
  /* Unpack data in leafbuf to leafdata for remote communication */
  PetscCall(PetscSFLinkUnpackLeafData(sf, link, PETSCSF_REMOTE, leafdata, op));
  /* Recycle the link */
  PetscCall(PetscSFLinkReclaim(sf, &link));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Shared by ReduceBegin and FetchAndOpBegin */
static inline PetscErrorCode PetscSFLeafToRootBegin_Basic(PetscSF sf, MPI_Datatype unit, PetscMemType leafmtype, const void *leafdata, PetscMemType rootmtype, void *rootdata, MPI_Op op, PetscSFOperation sfop, PetscSFLink *out)
{
  PetscSFLink link = NULL;

  PetscFunctionBegin;
  PetscCall(PetscSFLinkCreate(sf, unit, rootmtype, rootdata, leafmtype, leafdata, op, sfop, &link));
  PetscCall(PetscSFLinkPackLeafData(sf, link, PETSCSF_REMOTE, leafdata));
  PetscCall(PetscSFLinkStartCommunication(sf, link, PETSCSF_LEAF2ROOT));
  *out = link;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* leaf -> root with reduction */
PETSC_INTERN PetscErrorCode PetscSFReduceBegin_Basic(PetscSF sf, MPI_Datatype unit, PetscMemType leafmtype, const void *leafdata, PetscMemType rootmtype, void *rootdata, MPI_Op op)
{
  PetscSFLink link = NULL;

  PetscFunctionBegin;
  PetscCall(PetscSFLeafToRootBegin_Basic(sf, unit, leafmtype, leafdata, rootmtype, rootdata, op, PETSCSF_REDUCE, &link));
  PetscCall(PetscSFLinkScatterLocal(sf, link, PETSCSF_LEAF2ROOT, rootdata, (void *)leafdata, op));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFReduceEnd_Basic(PetscSF sf, MPI_Datatype unit, const void *leafdata, void *rootdata, MPI_Op op)
{
  PetscSFLink link = NULL;

  PetscFunctionBegin;
  PetscCall(PetscSFLinkGetInUse(sf, unit, rootdata, leafdata, PETSC_OWN_POINTER, &link));
  PetscCall(PetscSFLinkFinishCommunication(sf, link, PETSCSF_LEAF2ROOT));
  PetscCall(PetscSFLinkUnpackRootData(sf, link, PETSCSF_REMOTE, rootdata, op));
  PetscCall(PetscSFLinkReclaim(sf, &link));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFFetchAndOpBegin_Basic(PetscSF sf, MPI_Datatype unit, PetscMemType rootmtype, void *rootdata, PetscMemType leafmtype, const void *leafdata, void *leafupdate, MPI_Op op)
{
  PetscSFLink link = NULL;

  PetscFunctionBegin;
  PetscCall(PetscSFLeafToRootBegin_Basic(sf, unit, leafmtype, leafdata, rootmtype, rootdata, op, PETSCSF_FETCH, &link));
  PetscCall(PetscSFLinkFetchAndOpLocal(sf, link, rootdata, leafdata, leafupdate, op));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFFetchAndOpEnd_Basic(PetscSF sf, MPI_Datatype unit, void *rootdata, const void *leafdata, void *leafupdate, MPI_Op op)
{
  PetscSFLink link = NULL;

  PetscFunctionBegin;
  PetscCall(PetscSFLinkGetInUse(sf, unit, rootdata, leafdata, PETSC_OWN_POINTER, &link));
  /* This implementation could be changed to unpack as receives arrive, at the cost of non-determinism */
  PetscCall(PetscSFLinkFinishCommunication(sf, link, PETSCSF_LEAF2ROOT));
  /* Do fetch-and-op, the (remote) update results are in rootbuf */
  PetscCall(PetscSFLinkFetchAndOpRemote(sf, link, rootdata, op));
  /* Bcast rootbuf to leafupdate */
  PetscCall(PetscSFLinkStartCommunication(sf, link, PETSCSF_ROOT2LEAF));
  PetscCall(PetscSFLinkFinishCommunication(sf, link, PETSCSF_ROOT2LEAF));
  /* Unpack and insert fetched data into leaves */
  PetscCall(PetscSFLinkUnpackLeafData(sf, link, PETSCSF_REMOTE, leafupdate, MPI_REPLACE));
  PetscCall(PetscSFLinkReclaim(sf, &link));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode PetscSFGetLeafRanks_Basic(PetscSF sf, PetscMPIInt *niranks, const PetscMPIInt **iranks, const PetscInt **ioffset, const PetscInt **irootloc)
{
  PetscSF_Basic *bas = (PetscSF_Basic *)sf->data;

  PetscFunctionBegin;
  if (niranks) *niranks = bas->niranks;
  if (iranks) *iranks = bas->iranks;
  if (ioffset) *ioffset = bas->ioffset;
  if (irootloc) *irootloc = bas->irootloc;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* An optimized PetscSFCreateEmbeddedRootSF. We aggressively make use of the established communication on sf.
   We need one bcast on sf, and no communication anymore to build the embedded sf. Note that selected[]
   was sorted before calling the routine.
 */
PETSC_INTERN PetscErrorCode PetscSFCreateEmbeddedRootSF_Basic(PetscSF sf, PetscInt nselected, const PetscInt *selected, PetscSF *newsf)
{
  PetscSF            esf;
  PetscInt          *esf_roffset, *esf_rmine, *esf_rremote;
  PetscInt           j, p, q, nroots, esf_nleaves, *new_ilocal, minleaf, maxleaf, maxlocal;
  char              *rootdata, *leafdata, *leafmem; /* Only stores 0 or 1, so we can save memory with char */
  PetscMPIInt       *esf_ranks, nranks, ndranks, niranks, esf_nranks, esf_ndranks, ndiranks;
  const PetscMPIInt *ranks, *iranks;
  const PetscInt    *roffset, *rmine, *rremote, *ioffset, *irootloc;
  PetscBool          connected;
  PetscSFNode       *new_iremote;
  PetscSF_Basic     *bas;

  PetscFunctionBegin;
  PetscCall(PetscSFCreate(PetscObjectComm((PetscObject)sf), &esf));
  PetscCall(PetscSFSetFromOptions(esf));
  PetscCall(PetscSFSetType(esf, PETSCSFBASIC)); /* This optimized routine can only create a basic sf */

  /* Find out which leaves are still connected to roots in the embedded sf by doing a Bcast */
  PetscCall(PetscSFGetGraph(sf, &nroots, NULL, NULL, NULL));
  PetscCall(PetscSFGetLeafRange(sf, &minleaf, &maxleaf));
  maxlocal = maxleaf - minleaf + 1;
  PetscCall(PetscCalloc2(nroots, &rootdata, maxlocal, &leafmem));
  leafdata = PetscSafePointerPlusOffset(leafmem, -minleaf);
  /* Tag selected roots */
  for (PetscInt i = 0; i < nselected; ++i) rootdata[selected[i]] = 1;

  PetscCall(PetscSFBcastBegin(sf, MPI_CHAR, rootdata, leafdata, MPI_REPLACE));
  PetscCall(PetscSFBcastEnd(sf, MPI_CHAR, rootdata, leafdata, MPI_REPLACE));
  PetscCall(PetscSFGetLeafInfo_Basic(sf, &nranks, &ndranks, &ranks, &roffset, &rmine, &rremote)); /* Get send info */
  esf_nranks = esf_ndranks = esf_nleaves = 0;
  for (PetscMPIInt i = 0; i < nranks; i++) {
    connected = PETSC_FALSE; /* Is this process still connected to this remote root rank? */
    for (j = roffset[i]; j < roffset[i + 1]; j++) {
      if (leafdata[rmine[j]]) {
        esf_nleaves++;
        connected = PETSC_TRUE;
      }
    }
    if (connected) {
      esf_nranks++;
      if (i < ndranks) esf_ndranks++;
    }
  }

  /* Set graph of esf and also set up its outgoing communication (i.e., send info), which is usually done by PetscSFSetUpRanks */
  PetscCall(PetscMalloc1(esf_nleaves, &new_ilocal));
  PetscCall(PetscMalloc1(esf_nleaves, &new_iremote));
  PetscCall(PetscMalloc4(esf_nranks, &esf_ranks, esf_nranks + 1, &esf_roffset, esf_nleaves, &esf_rmine, esf_nleaves, &esf_rremote));
  p              = 0; /* Counter for connected root ranks */
  q              = 0; /* Counter for connected leaves */
  esf_roffset[0] = 0;
  for (PetscMPIInt i = 0; i < nranks; i++) { /* Scan leaf data again to fill esf arrays */
    connected = PETSC_FALSE;
    for (j = roffset[i]; j < roffset[i + 1]; j++) {
      if (leafdata[rmine[j]]) {
        esf_rmine[q] = new_ilocal[q] = rmine[j];
        esf_rremote[q]               = rremote[j];
        new_iremote[q].index         = rremote[j];
        new_iremote[q].rank          = ranks[i];
        connected                    = PETSC_TRUE;
        q++;
      }
    }
    if (connected) {
      esf_ranks[p]       = ranks[i];
      esf_roffset[p + 1] = q;
      p++;
    }
  }

  /* SetGraph internally resets the SF, so we only set its fields after the call */
  PetscCall(PetscSFSetGraph(esf, nroots, esf_nleaves, new_ilocal, PETSC_OWN_POINTER, new_iremote, PETSC_OWN_POINTER));
  esf->nranks    = esf_nranks;
  esf->ndranks   = esf_ndranks;
  esf->ranks     = esf_ranks;
  esf->roffset   = esf_roffset;
  esf->rmine     = esf_rmine;
  esf->rremote   = esf_rremote;
  esf->nleafreqs = esf_nranks - esf_ndranks;

  /* Set up the incoming communication (i.e., recv info) stored in esf->data, which is usually done by PetscSFSetUp_Basic */
  bas = (PetscSF_Basic *)esf->data;
  PetscCall(PetscSFGetRootInfo_Basic(sf, &niranks, &ndiranks, &iranks, &ioffset, &irootloc)); /* Get recv info */
  /* Embedded sf always has simpler communication than the original one. We might allocate longer arrays than needed here. But we
     we do not care since these arrays are usually short. The benefit is we can fill these arrays by just parsing irootloc once.
   */
  PetscCall(PetscMalloc2(niranks, &bas->iranks, niranks + 1, &bas->ioffset));
  PetscCall(PetscMalloc1(ioffset[niranks], &bas->irootloc));
  bas->niranks = bas->ndiranks = bas->ioffset[0] = 0;
  p                                              = 0; /* Counter for connected leaf ranks */
  q                                              = 0; /* Counter for connected roots */
  for (PetscMPIInt i = 0; i < niranks; i++) {
    connected = PETSC_FALSE; /* Is the current process still connected to this remote leaf rank? */
    for (j = ioffset[i]; j < ioffset[i + 1]; j++) {
      if (rootdata[irootloc[j]]) {
        bas->irootloc[q++] = irootloc[j];
        connected          = PETSC_TRUE;
      }
    }
    if (connected) {
      bas->niranks++;
      if (i < ndiranks) bas->ndiranks++; /* Note that order of ranks (including distinguished ranks) is kept */
      bas->iranks[p]      = iranks[i];
      bas->ioffset[p + 1] = q;
      p++;
    }
  }
  bas->itotal     = q;
  bas->nrootreqs  = bas->niranks - bas->ndiranks;
  esf->persistent = PETSC_TRUE;
  /* Setup packing related fields */
  PetscCall(PetscSFSetUpPackFields(esf));

  /* Copy from PetscSFSetUp(), since this method wants to skip PetscSFSetUp(). */
#if defined(PETSC_HAVE_CUDA)
  if (esf->backend == PETSCSF_BACKEND_CUDA) {
    esf->ops->Malloc = PetscSFMalloc_CUDA;
    esf->ops->Free   = PetscSFFree_CUDA;
  }
#endif

#if defined(PETSC_HAVE_HIP)
  /* TODO: Needs debugging */
  if (esf->backend == PETSCSF_BACKEND_HIP) {
    esf->ops->Malloc = PetscSFMalloc_HIP;
    esf->ops->Free   = PetscSFFree_HIP;
  }
#endif

#if defined(PETSC_HAVE_KOKKOS)
  if (esf->backend == PETSCSF_BACKEND_KOKKOS) {
    esf->ops->Malloc = PetscSFMalloc_Kokkos;
    esf->ops->Free   = PetscSFFree_Kokkos;
  }
#endif
  esf->setupcalled = PETSC_TRUE; /* We have done setup ourselves! */
  PetscCall(PetscFree2(rootdata, leafmem));
  *newsf = esf;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_EXTERN PetscErrorCode PetscSFCreate_Basic(PetscSF sf)
{
  PetscSF_Basic *dat;

  PetscFunctionBegin;
  sf->ops->SetUp                = PetscSFSetUp_Basic;
  sf->ops->Reset                = PetscSFReset_Basic;
  sf->ops->Destroy              = PetscSFDestroy_Basic;
  sf->ops->View                 = PetscSFView_Basic;
  sf->ops->BcastBegin           = PetscSFBcastBegin_Basic;
  sf->ops->BcastEnd             = PetscSFBcastEnd_Basic;
  sf->ops->ReduceBegin          = PetscSFReduceBegin_Basic;
  sf->ops->ReduceEnd            = PetscSFReduceEnd_Basic;
  sf->ops->FetchAndOpBegin      = PetscSFFetchAndOpBegin_Basic;
  sf->ops->FetchAndOpEnd        = PetscSFFetchAndOpEnd_Basic;
  sf->ops->GetLeafRanks         = PetscSFGetLeafRanks_Basic;
  sf->ops->CreateEmbeddedRootSF = PetscSFCreateEmbeddedRootSF_Basic;
  sf->ops->SetCommunicationOps  = PetscSFSetCommunicationOps_Basic;

  sf->persistent = PETSC_TRUE; // currently SFBASIC always uses persistent send/recv
  sf->collective = PETSC_FALSE;

  PetscCall(PetscNew(&dat));
  sf->data = (void *)dat;
  PetscFunctionReturn(PETSC_SUCCESS);
}
