// Declare variables
PetscMPIInt commSize;
PetscViewer petscviewer;
PetscReal abstol, rtol, stol;
// Initialize PETSc
PetscFunctionBeginUser;
PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));
PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &commSize));
PetscCheck(commSize == 1, PETSC_COMM_SELF, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!");
