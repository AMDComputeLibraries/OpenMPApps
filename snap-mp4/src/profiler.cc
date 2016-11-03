
#include "profiler.h"

double wtime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1.0E-6;
}

void outer_profiler(struct timers * timers)
{
    if (!profiling)
        return;
#if 0
    cl_int err;

    // Times are in nanoseconds
    cl_ulong tick, tock;

    // Get outer souce update times
    err = clGetEventProfilingInfo(outer_source_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
    check_ocl(err, "Getting outer source start time");
    err = clGetEventProfilingInfo(outer_source_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
    check_ocl(err, "Getting outer source end time");
    timers->outer_source_time += (double)(tock - tick) * 1.0E-9;

    // Get outer parameter times
    // Start is velocity delta start, end is denominator end
    err = clGetEventProfilingInfo(velocity_delta_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
    check_ocl(err, "Getting velocity delta start time");
    err = clGetEventProfilingInfo(denominator_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
    check_ocl(err, "Getting denominator end time");
    timers->outer_params_time += (double)(tock - tick) * 1.0E-9;
#endif
}

void inner_profiler(struct timers * timers, struct problem * problem)
{
    if (!profiling)
        return;
#if 0
    cl_int err;

    // Times are in nanoseconds
    cl_ulong tick, tock;

    // Get inner source update times
    err = clGetEventProfilingInfo(inner_source_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
    check_ocl(err, "Getting inner source start time");
    err = clGetEventProfilingInfo(inner_source_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
    check_ocl(err, "Getting inner source end time");
    timers->inner_source_time += (double)(tock - tick) * 1.0E-9;

    // Get scalar flux reduction times
    err = clGetEventProfilingInfo(scalar_flux_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
    check_ocl(err, "Getting scalar flux start time");
    err = clGetEventProfilingInfo(scalar_flux_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
    check_ocl(err, "Getting scalar flux end time");
    timers->reduction_time += (double)(tock - tick) * 1.0E-9;
    if (problem->cmom-1 > 0)
    {
        err = clGetEventProfilingInfo(scalar_flux_moments_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
        check_ocl(err, "Getting scalar flux moments start time");
        err = clGetEventProfilingInfo(scalar_flux_moments_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
        check_ocl(err, "Getting scalar flux moments end time");
        timers->reduction_time += (double)(tock - tick) * 1.0E-9;
    }
#endif
}


void chunk_profiler(struct timers * timers)
{
    if (!profiling)
        return;
#if 0
    cl_int err;

    // Times are in nanoseconds
    cl_ulong tick, tock;

    // Get recv writes
    if (flux_i_write_event)
    {
        err = clGetEventProfilingInfo(flux_i_write_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
        check_ocl(err, "Getting flux i write start time");
        err = clGetEventProfilingInfo(flux_i_write_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
        check_ocl(err, "Getting flux i write stop time");
        timers->sweep_transfer_time += (double)(tock - tick) * 1.0E-9;
    }

    if (flux_j_write_event)
    {
        err = clGetEventProfilingInfo(flux_j_write_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
        check_ocl(err, "Getting flux j write start time");
        err = clGetEventProfilingInfo(flux_j_write_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
        check_ocl(err, "Getting flux j write stop time");
        timers->sweep_transfer_time += (double)(tock - tick) * 1.0E-9;
    }

    // Get send reads
    if (flux_i_read_event)
    {
        err = clGetEventProfilingInfo(flux_i_read_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
        check_ocl(err, "Getting flux i read start time");
        err = clGetEventProfilingInfo(flux_i_read_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
        check_ocl(err, "Getting flux i read stop time");
        timers->sweep_transfer_time += (double)(tock - tick) * 1.0E-9;
    }

    if (flux_j_read_event)
    {
        err = clGetEventProfilingInfo(flux_j_read_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tick, NULL);
        check_ocl(err, "Getting flux j read start time");
        err = clGetEventProfilingInfo(flux_j_read_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tock, NULL);
        check_ocl(err, "Getting flux j read stop time");
        timers->sweep_transfer_time += (double)(tock - tick) * 1.0E-9;
    }
#endif
}

