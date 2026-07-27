# Runtime activation for the `aws` pixi environment (AWS ParallelCluster).
#
# The aws env installs conda-forge's `openmpi external_*` stub, so MPI comes
# from the SYSTEM OpenMPI 5 + EFA under /opt/amazon. Put it on PATH (h5pcc
# resolves bare `mpicc` from PATH at build time) and LD_LIBRARY_PATH (conda
# libhdf5/mpi4py resolve libmpi.so.40 at runtime). Guarded so it is a no-op
# off-cluster; modeled on Thea_Neutronics tools/aws-activation.sh.
if [ -d /opt/amazon/openmpi5 ]; then
    if ! command -v module >/dev/null 2>&1; then
        [ -r /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
    fi
    module use /opt/amazon/modules/modulefiles 2>/dev/null
    module load openmpi5 libfabric-aws 2>/dev/null
    export PATH="/opt/amazon/openmpi5/bin${PATH:+:$PATH}"
    export LD_LIBRARY_PATH="/opt/amazon/openmpi5/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

    # hdf5's h5pcc bakes an absolute CCBASE=<env>/bin/mpicc, which the external
    # openmpi stub never provides; override so h5pcc resolves the system mpicc
    # from the PATH set above. On-cluster only, so off-cluster builds still
    # fail fast instead of silently linking a random local MPI.
    export HDF5_CC=mpicc
    export HDF5_CLINKER=mpicc
fi
