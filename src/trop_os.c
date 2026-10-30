/*
 * Operating system functions for the fasttrop C core. This file is separate,
 * so <windows.h> and the R headers are not in the same translation unit.
 *
 * trop_wtime() is a high-resolution monotonic clock. On Windows, the
 * omp_get_wtime() function of MinGW has a resolution of only one millisecond.
 * Thus, trop_wtime() uses QueryPerformanceCounter().
 *
 * On Windows, trop_qos_enter() and trop_qos_leave() stop execution-speed power
 * throttling (EcoQoS) for the thread while the C code runs. Then Windows
 * controls the thread again. The calls can nest. Windows 11 uses power
 * throttling on processes that run in the background for a long time, for
 * example an R session with a long calculation. Power throttling can make the
 * C code slower. On other operating systems, or with throttle = TRUE, these
 * functions have no effect.
 */
#ifdef _WIN32
#if !defined(_WIN32_WINNT) || _WIN32_WINNT < 0x0602
#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0602
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <time.h>
#endif

#include "trop_core.h"

double trop_wtime(void) {
#ifdef _WIN32
  LARGE_INTEGER f, c;
  QueryPerformanceFrequency(&f);
  QueryPerformanceCounter(&c);
  return (double) c.QuadPart / (double) f.QuadPart;
#else
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
#endif
}

#ifdef _WIN32
/* Some MinGW-w64 versions do not declare this type. The layout is from the
 * Windows SDK. */
typedef struct {
  ULONG Version;
  ULONG ControlMask;
  ULONG StateMask;
} trop_thread_throttling_state;
#define TROP_THROTTLING_VERSION 1UL
#define TROP_THROTTLING_EXECUTION_SPEED 0x1UL
#define TROP_THREAD_POWER_THROTTLING ((THREAD_INFORMATION_CLASS) 3) /* ThreadPowerThrottling */

static _Thread_local int qos_depth = 0;

static void set_throttling(ULONG control_mask) {
  trop_thread_throttling_state s;
  s.Version = TROP_THROTTLING_VERSION;
  s.ControlMask = control_mask; /* 0 gives the control to Windows again */
  s.StateMask = 0;              /* with ControlMask: no power throttling */
  SetThreadInformation(GetCurrentThread(), TROP_THREAD_POWER_THROTTLING, &s, sizeof s);
}
#endif

void trop_qos_enter(const trop_opts *o) {
#ifdef _WIN32
  if (o->throttle) return;
  if (qos_depth++ == 0) set_throttling(TROP_THROTTLING_EXECUTION_SPEED);
#else
  (void) o;
#endif
}

void trop_qos_leave(const trop_opts *o) {
#ifdef _WIN32
  if (o->throttle) return;
  if (qos_depth > 0 && --qos_depth == 0) set_throttling(0UL);
#else
  (void) o;
#endif
}
